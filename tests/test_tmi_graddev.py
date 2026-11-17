"""Component and spatial conventions for native and scanner-frame GNC fields."""

import subprocess

import nibabel as nib
import numpy as np
import pytest
from mrtrix3 import MRtrixError, run

from lib.gnc.backends import _write_graddev
from lib.gnc.graddev import load_gnc_field
from lib.io import load_mrtrix


def read_dwi(shape, affine):
    nib.save(nib.Nifti1Image(np.zeros(shape, np.float32), affine), "dwi.nii")
    run.command(["mrconvert", "dwi.nii", "dwi.mif"], show=False)
    return load_mrtrix("dwi.mif")


@pytest.mark.parametrize("handedness", [-1, 1])
@pytest.mark.parametrize("backend", ["tortoise", "gradunwarp"])
def test_native_graddev_maps_coil_to_spatial_direction(tmp_path, monkeypatch, handedness, backend):
    monkeypatch.chdir(tmp_path)
    # Permuted, oblique, anisotropic axes. MRtrix reorders the spatial samples.
    axes = np.array([[0., -0.6, 0.8], [handedness, 0., 0.], [0., 0.8, 0.6]])
    affine = nib.affines.from_matvec(axes @ np.diag([1., 2., 3.]), [3., 4., 5.])
    shape = (3, 4, 5)
    dwi = read_dwi(shape, affine)
    # Coil i acquires a component along voxel axis k. Vary its size by position.
    shear = np.arange(np.prod(shape)).reshape(shape) / 512.
    raw = np.broadcast_to(np.eye(3), shape + (3, 3)).copy()
    raw[..., 0, 2] = shear
    raw[0, 0, 0] = np.nan
    field = raw if backend == "tortoise" else raw.reshape(*shape, 9)
    nib.save(nib.Nifti1Image(field, affine), "native field.nii")
    # Column i of L gives the gradient created by coil i. The first FSL
    # axis points opposite to voxel i on a right-handed grid.
    coil_i = axes[:, 0]
    if backend == "gradunwarp" and handedness == 1:
        coil_i = -coil_i
    expected = np.eye(3) + shear[..., None, None] * np.outer(axes[:, 2], coil_i)
    scanner = nib.Nifti1Image(expected.reshape(*shape, 9), affine)
    nib.save(scanner, "scanner.nii.gz")
    expected_ras = nib.as_closest_canonical(scanner).get_fdata().reshape(*dwi.shape[:3], 3, 3)

    np.testing.assert_allclose(load_gnc_field("native field.nii", dwi, backend), expected_ras, atol=1e-7)
    _write_graddev("native field.nii", "tagged.nii.gz", backend)
    np.testing.assert_allclose(load_gnc_field("tagged.nii.gz", dwi), expected_ras, atol=1e-7)
    np.testing.assert_allclose(load_gnc_field("scanner.nii.gz", dwi, "scanner"), expected_ras)
    run.command(["mrconvert", "scanner.nii.gz", "scanner.mif"], show=False)
    np.testing.assert_allclose(load_gnc_field("scanner.mif", dwi), expected_ras)


@pytest.mark.parametrize("filename,description", [
    ("preproc_graddev.nii", b""),
    ("tortoise_graddev.nii.gz", b"TORTOISE"),
])
@pytest.mark.parametrize("x_spacing", [-2., 2.])
def test_untagged_nifti_default_depends_on_grid(tmp_path, monkeypatch, filename, description, x_spacing):
    monkeypatch.chdir(tmp_path)
    affine = np.diag([x_spacing, 2., 2., 1.])
    dwi = read_dwi((3, 4, 5), affine)
    # Asymmetric entries expose a wrong transpose or frame.
    J = np.eye(3)
    J[0, 1] = 0.25
    J[2, 0] = -0.1
    field = nib.Nifti1Image(np.broadcast_to(J.ravel(), dwi.shape + (9,)), affine)
    field.header["descrip"] = description
    nib.save(field, filename)
    tortoise = load_gnc_field(filename, dwi, "tortoise")
    gradunwarp = load_gnc_field(filename, dwi, "gradunwarp")
    if x_spacing < 0:
        np.testing.assert_array_equal(load_gnc_field(filename, dwi), gradunwarp)
        np.testing.assert_array_equal(tortoise, gradunwarp)
    else:
        with pytest.raises(MRtrixError, match="Supply -gnc_format tortoise"):
            load_gnc_field(filename, dwi)
        assert tortoise[0, 0, 0, 1, 0] == 0.25
        assert gradunwarp[0, 0, 0, 1, 0] == -0.25


def test_explicit_format_overrides_header_tag(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    dwi = read_dwi((3, 4, 5), np.eye(4))
    J = np.eye(3)
    J[0, 1] = 0.25
    field = np.broadcast_to(J.ravel(), dwi.shape + (9,))
    nib.save(nib.Nifti1Image(field, np.eye(4)), "native.nii")
    _write_graddev("native.nii", "renamed.nii", "gradunwarp")
    # On a positive-determinant grid, the FSL x flip changes the sign of the shear.
    assert load_gnc_field("renamed.nii", dwi)[0, 0, 0, 1, 0] == -0.25
    assert load_gnc_field("renamed.nii", dwi, "tortoise")[0, 0, 0, 1, 0] == 0.25


def test_graddev_grid_must_match(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    dwi = read_dwi((3, 4, 5), np.eye(4))
    shifted = nib.affines.from_matvec(np.eye(3), [1., 0., 0.])
    field = np.broadcast_to(np.eye(3).ravel(), dwi.shape + (9,))
    nib.save(nib.Nifti1Image(field, shifted), "shifted.nii")
    with pytest.raises(MRtrixError, match="same spatial grid"):
        load_gnc_field("shifted.nii", dwi, "tortoise")


def test_native_graddev_recovers_fiber_direction(tmp_path, monkeypatch):
    """A known 3-degree coil rotation must recover a fiber along scanner y."""
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("LOKY_MAX_CPU_COUNT", "1")
    c, s = np.cos(np.deg2rad(3.)), np.sin(np.deg2rad(3.))
    rotation = np.array([[1., 0., 0.], [0., c, -s], [0., s, c]])
    dirs = np.random.default_rng(123).normal(size=(36, 3))
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    bvecs = np.vstack([np.zeros((2, 3)), dirs, dirs])
    bvals = np.repeat([0., 1000., 2000.], [2, 36, 36])
    np.savetxt("dwi.bvec", bvecs.T)
    np.savetxt("dwi.bval", bvals[None])

    affine = np.diag([-2., 2., 2., 1.])
    nominal = bvecs * [-1., 1., 1.]  # FSL bvecs on this LAS grid to scanner RAS
    actual = nominal @ rotation.T
    eigenvalues = np.array([0.4, 1.5, 0.5])
    signal = 100 * np.exp(-bvals / 1000 * (actual**2 @ eigenvalues))
    shape = (3, 3, 3)
    nib.save(nib.Nifti1Image(np.broadcast_to(signal, shape + signal.shape), affine), "dwi.nii")
    # Reflection of x leaves a rotation about x unchanged. Native rows are coils.
    native = np.broadcast_to(rotation.T.ravel(), shape + (9,))
    nib.save(nib.Nifti1Image(native, affine), "preproc_graddev.nii")

    cases = [("nominal", [], 3.), ("native", ["-gnc", "preproc_graddev.nii", "-gnc_format", "tortoise"], 0.),
             ("legacy", ["-gnc", "preproc_graddev.nii", "-gnc_format", "scanner"], 6.)]
    for output, options, expected_angle in cases:
        subprocess.run(["tmi", "-DTI", "-DKI", "-fit_constraints", "0,0,0", "-n_cores", "1",
                        *options, "dwi.nii", output], check=True)
        for model in ("dti", "dki"):
            # FE stores abs(FA * principal eigenvector). Its y component gives
            # the unsigned angle to the true fiber after normalization.
            fe = nib.load(f"{output}/fe_{model}.nii").get_fdata()
            cos_angle = fe[..., 1] / np.linalg.norm(fe, axis=-1)
            np.testing.assert_allclose(np.degrees(np.arccos(np.clip(cos_angle, 0., 1.))),
                                       expected_angle, atol=1e-3)

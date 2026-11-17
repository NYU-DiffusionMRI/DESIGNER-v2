"""Component and spatial conventions for native and scanner-frame GNC fields."""

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

    mode = None if backend == "tortoise" else "gradunwarp"
    np.testing.assert_allclose(load_gnc_field("native field.nii", dwi, mode), expected_ras, atol=1e-7)
    _write_graddev("native field.nii", "tagged.nii.gz", backend)
    np.testing.assert_allclose(load_gnc_field("tagged.nii.gz", dwi), expected_ras, atol=1e-7)
    np.testing.assert_allclose(load_gnc_field("scanner.nii.gz", dwi, "scanner"), expected_ras)
    run.command(["mrconvert", "scanner.nii.gz", "scanner.mif"], show=False)
    np.testing.assert_allclose(load_gnc_field("scanner.mif", dwi), expected_ras)


def test_graddev_grid_must_match(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    dwi = read_dwi((3, 4, 5), np.eye(4))
    shifted = nib.affines.from_matvec(np.eye(3), [1., 0., 0.])
    field = np.broadcast_to(np.eye(3).ravel(), dwi.shape + (9,))
    nib.save(nib.Nifti1Image(field, shifted), "shifted.nii")
    with pytest.raises(MRtrixError, match="same spatial grid"):
        load_gnc_field("shifted.nii", dwi)

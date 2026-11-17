"""Check the Siemens displacement and graddev against one analytic coil."""

import shutil

import nibabel as nib
import numpy as np
import pytest

from lib.gnc.backends import compute_gnc_warp, derive_gnc_graddev
from lib.gnc.backends.tortoise import _center_siemens_displacement
from lib.gnc.graddev import load_gnc_field
from lib.io import load_mrtrix
from mrtrix3 import run


def test_voxel_center_offset_uses_physical_diagonal():
    # A permuted grid with both transverse diagonals zero has no offset.
    axes = np.array([[0., 0., -3.], [1., 0., 0.], [0., 2., 0.]])
    affine = nib.affines.from_matvec(axes, [30., 45., 50.])
    field = np.random.default_rng(0).normal(size=(4, 5, 6, 3))
    np.testing.assert_array_equal(_center_siemens_displacement(field, affine), field)


def test_siemens_displacement_interpolation_preserves_axes_and_edges():
    # A tilt about x puts part of the physical y offset along voxel k.
    axes = np.array([[-2., 0., 0.], [0., 1.6, -1.8], [0., 1.2, 2.4]])
    affine = nib.affines.from_matvec(axes, [30., 45., 50.])
    grid = np.moveaxis(np.indices((6, 7, 8)), 0, -1)
    centers = nib.affines.apply_affine(affine, grid)

    def displacement(ras):
        x, y, z = np.moveaxis(ras, -1, 0)
        return np.stack([x*y, y**2*z, z*x], axis=-1) / 1000

    sampled = displacement(centers + [-1., 0.8, 0.])
    np.testing.assert_allclose(_center_siemens_displacement(sampled, affine),
                               displacement(centers), atol=1e-10)


@pytest.mark.skipif(not shutil.which("CreateGradientNonlinearityBMatrix"), reason="requires TORTOISE")
@pytest.mark.parametrize("x_spacing", [-1.5, 1.5])
def test_siemens_warp_and_graddev_at_voxel_centers(tmp_path, monkeypatch, x_spacing):
    monkeypatch.chdir(tmp_path)
    # One Siemens A(3,0) harmonic on z has a closed form in scanner mm.
    coeff = tmp_path / "analytic.grad"
    coeff.write_text(" 0.2 m = R0\n  1 A( 3, 0)                 0.5                      z\n")
    axes = np.array([[0.8, -0.6, 0.], [0.6, 0.8, 0.], [0., 0., 1.]])
    affine = nib.affines.from_matvec(axes @ np.diag([x_spacing, 2., 2.5]), [40., 34., 52.])
    shape = (8, 9, 10)
    ref = nib.Nifti1Image(np.zeros(shape, np.float32), affine)
    nib.save(ref, "ref.nii")
    warp_file = compute_gnc_warp("ref.nii", coeff, "siemens", "test", "tortoise")
    grad_file = derive_gnc_graddev("ref.nii", coeff, "siemens", "test", "tortoise")

    ijk = np.moveaxis(np.indices(shape), 0, -1)
    x, y, z = np.moveaxis(nib.affines.apply_affine(affine, ijk), -1, 0)
    k = 0.5 / 200**2
    dz = k * (z**3 - 1.5*z*(x*x + y*y))
    # With this in-plane rotation, the FSL and RAS z axes coincide.
    warp = nib.load(warp_file)
    np.testing.assert_allclose(warp.affine, affine, atol=1e-5)
    np.testing.assert_allclose(warp.get_fdata()[..., 2], dz, atol=8e-5)
    np.testing.assert_allclose(warp.get_fdata()[..., :2], 0., atol=8e-5)

    # The gradients produced by the z coil are the derivatives of z + dz.
    L = np.broadcast_to(np.eye(3), shape + (3, 3)).copy()
    L[..., 0, 2] = -3*k*x*z
    L[..., 1, 2] = -3*k*y*z
    L[..., 2, 2] += k*(3*z*z - 1.5*(x*x + y*y))
    expected = nib.as_closest_canonical(nib.Nifti1Image(L, affine)).get_fdata()
    run.command(["mrconvert", "ref.nii", "ref.mif"], show=False)
    np.testing.assert_allclose(load_gnc_field(grad_file, load_mrtrix("ref.mif")), expected, atol=2e-6)

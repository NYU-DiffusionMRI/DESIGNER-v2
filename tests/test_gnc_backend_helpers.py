"""Unit tests for the pure-numpy GNC backend helpers.

No scanner data, no FSL/TORTOISE binaries, no scratch dir - just numpy/nibabel.
"""

import numpy as np
import nibabel as nib
import pytest

from lib.gnc.backends import _sanitize_graddev_matrix
from lib.gnc.backends.tortoise import _lps_disp_to_fsl_rel, resolve_is_ge


def test_lps_to_fsl_axis_aligned_matches_empirical_case():
    # MNI152-shaped affine diag(-1,1,1): empirically "flip Y only" against
    # gradunwarp's own output.
    affine = nib.affines.from_matvec(np.diag([-1, 1, 1]), [90, -126, -72])
    disp = np.random.default_rng(0).normal(size=(4, 4, 4, 3))
    np.testing.assert_allclose(_lps_disp_to_fsl_rel(disp, affine), disp * [1, -1, 1])


def test_lps_to_fsl_oblique_round_trips_to_world():
    # ~1 deg in-plane rotation, 2 mm voxels, det < 0 (like real scanner data).
    th = np.deg2rad(1.0)
    rot = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
    M = rot @ np.diag([-2.0, 2.0, 2.0])
    affine = nib.affines.from_matvec(M, [116, -120, -56])
    disp = np.random.default_rng(1).normal(size=(3, 3, 3, 3))
    fsl = _lps_disp_to_fsl_rel(disp, affine)
    # undo FSL's mm scaling (det < 0 -> no x flip), map voxel steps to world RAS
    world_ras = (fsl / np.linalg.norm(M, axis=0)) @ M.T
    np.testing.assert_allclose(world_ras, disp * [-1, -1, 1], atol=1e-12)


def test_lps_to_fsl_positive_determinant_negates_x():
    # RAS-identity storage (det > 0): +1 mm toward Left (LPS +x) is -1 voxel
    # along i, which FSL's radiological internal frame reports as +1 mm.
    fsl = _lps_disp_to_fsl_rel(np.array([[1.0, 0.0, 0.0]]), np.eye(4))
    np.testing.assert_allclose(fsl, [[1.0, 0.0, 0.0]])


def test_sanitize_graddev_matrix_replaces_full_block_3x3():
    mat = np.arange(2 * 2 * 2 * 9, dtype=float).reshape(2, 2, 2, 3, 3)
    mat[1, 0, 1] = np.nan
    out = _sanitize_graddev_matrix(mat)
    assert np.array_equal(out[1, 0, 1], np.eye(3))
    # every other voxel is untouched
    untouched = mat.copy()
    untouched[1, 0, 1] = np.eye(3)
    assert np.array_equal(out, untouched)


def test_sanitize_graddev_matrix_partial_nan_zeroes_whole_block():
    mat = np.ones((1, 1, 1, 3, 3))
    mat[0, 0, 0, 1, 2] = np.nan  # only one entry is NaN, rest are valid
    out = _sanitize_graddev_matrix(mat)
    assert np.array_equal(out[0, 0, 0], np.eye(3))


def test_sanitize_graddev_matrix_accepts_flat_and_square_shapes():
    rng = np.random.default_rng(0)
    square = rng.random((2, 2, 2, 3, 3))
    flat = square.reshape(2, 2, 2, 9)
    out_square = _sanitize_graddev_matrix(square)
    out_flat = _sanitize_graddev_matrix(flat).reshape(2, 2, 2, 3, 3)
    assert np.array_equal(out_square, out_flat)


@pytest.mark.parametrize('scanner,expected', [
    ('siemens', False),
    ('ge', True),
])
def test_resolve_is_ge(scanner, expected):
    assert resolve_is_ge(scanner) is expected

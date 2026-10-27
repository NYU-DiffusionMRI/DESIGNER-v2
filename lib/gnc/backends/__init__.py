"""GNC backend dispatch. Each backend module implements GncBackend."""

from pathlib import Path
from typing import Protocol

import nibabel as nib
import numpy as np

from lib.gnc.backends import gradunwarp, tortoise
from lib.utils import force_remove


class GncBackend(Protocol):
    """Interface every backend module satisfies (structurally - modules don't
    inherit from this). Both functions take the same arguments; a backend
    ignores any it doesn't need."""

    def compute_warp(
        self,
        ref_volume_nii: str | Path,
        coeff_file: str | Path,
        scanner: str,
        out_prefix: str,
    ) -> str:
        """
        Compute the GNC geometric warp on the reference volume's grid.

        Args:
            ref_volume_nii: Single 3D reference volume (NIfTI) defining the grid
            coeff_file: Scanner gradient coefficient file
            scanner: Scanner vendor ('siemens' or 'ge')
            out_prefix: Basename prefix for output/intermediate files (in CWD)

        Returns:
            Path to {out_prefix}_warp_rel.nii.gz (FSL relative-warp convention)
        """
        ...

    def derive_graddev(
        self,
        ref_volume_nii: str | Path,
        coeff_file: str | Path,
        scanner: str,
        out_prefix: str,
    ) -> str:
        """
        Produce the raw per-voxel gradient-deviation tensor L.

        Called right after compute_warp() with the same arguments; may
        consume files compute_warp() deliberately left on disk.

        Args:
            ref_volume_nii: Same reference volume passed to compute_warp()
            coeff_file: Same coefficient file passed to compute_warp()
            scanner: Same scanner vendor passed to compute_warp()
            out_prefix: Same basename prefix passed to compute_warp()

        Returns:
            Path to a raw (unsanitized, may contain NaN) L tensor NIfTI with
            last dimension 9 (row-major 3x3, ~Identity at undistorted voxels).
            The caller sanitizes it and then deletes it.
        """
        ...


_BACKENDS: dict[str, GncBackend] = {
    'gradunwarp': gradunwarp,
    'tortoise': tortoise,
}


def _get_backend(backend: str) -> GncBackend:
    """
    Look up a backend module by its -gnc_backend name.

    Args:
        backend: Backend name, a key of _BACKENDS ('gradunwarp' or 'tortoise')

    Returns:
        The backend module

    Raises:
        MRtrixError: If backend is not a known name
    """
    from mrtrix3 import MRtrixError

    try:
        return _BACKENDS[backend]
    except KeyError:
        raise MRtrixError(f"Unknown GNC backend '{backend}'") from None


def compute_gnc_warp(
    ref_volume_nii: str | Path,
    coeff_file: str | Path,
    scanner: str,
    out_prefix: str,
    backend: str = 'gradunwarp',
) -> str:
    """
    Compute the relative FSL warp for gradient nonlinearity correction.

    Dispatches to the requested backend; both produce a warp in FSL's
    relative-warp convention (for use with applywarp), so nothing downstream
    (apply_warp_to_series, finalize_single_interpolation) needs to know which
    backend ran. Must be followed by a derive_gnc_graddev() call with the
    same ref_volume_nii/out_prefix/backend before any other GNC intermediates
    are produced - the tortoise backend deliberately leaves a byproduct file
    behind for that call to consume (the same kind of fixed-scratch-name
    coupling run_gnc_precorrection's 'gnc' basename already relies on).

    Args:
        ref_volume_nii: Path to the single reference volume to unwarp (the
            warp is static in scanner space, so one volume is enough)
        coeff_file: Path to the scanner gradient coefficient file
        scanner: Scanner vendor ('siemens' or 'ge')
        out_prefix: Basename prefix for the intermediate/output warp files
        backend: 'gradunwarp' (default) or 'tortoise'

    Returns:
        Path to the relative warp file (FSL convention, for use with applywarp)
    """
    return _get_backend(backend).compute_warp(ref_volume_nii, coeff_file, scanner, out_prefix)


def derive_gnc_graddev(
    ref_volume_nii: str | Path,
    coeff_file: str | Path,
    scanner: str,
    out_prefix: str,
    backend: str,
) -> str:
    """
    Produce {out_prefix}_graddev.nii.gz: the per-voxel gradient-deviation
    tensor L (L ~= Identity at undistorted voxels) consumed by
    designer2/tmi.py's -gnc option. Must be called immediately after
    compute_gnc_warp() with the same arguments.

    Args:
        ref_volume_nii: Same reference volume passed to compute_gnc_warp()
        coeff_file: Same coefficient file passed to compute_gnc_warp()
        scanner: Same scanner vendor passed to compute_gnc_warp()
        out_prefix: Same basename prefix passed to compute_gnc_warp()
        backend: 'gradunwarp' or 'tortoise' - same as compute_gnc_warp()

    Returns:
        Path to the sanitized graddev tensor field
    """
    raw = _get_backend(backend).derive_graddev(ref_volume_nii, coeff_file, scanner, out_prefix)

    out_path = f'{out_prefix}_graddev.nii.gz'
    _sanitize_graddev_nan(raw, out_path)
    force_remove(raw)
    return out_path


def _sanitize_graddev_matrix(mat):
    """Replace any 3x3 gradient-deviation block containing a NaN with the
    identity matrix as a whole block (not entrywise). Works around a
    confirmed TORTOISE bug: CreateGradientNonlinearityBMatrix emits NaN for
    the entire 3x3 block exactly at the isocenter voxel, where the true
    value is Identity. Accepts (...,9) or (...,3,3); returns the same shape.
    """
    flat = mat.reshape(*mat.shape[:-1], 3, 3) if mat.shape[-1] == 9 else mat
    bad = np.isnan(flat).any(axis=(-2, -1))
    flat = np.where(bad[..., None, None], np.eye(3, dtype=flat.dtype), flat)
    return flat.reshape(mat.shape)


def _sanitize_graddev_nan(raw_path: str | Path, out_path: str | Path) -> None:
    """
    Apply _sanitize_graddev_matrix() to a graddev NIfTI on disk.

    Args:
        raw_path: Raw L tensor NIfTI from a backend's derive_graddev()
        out_path: Where to write the sanitized tensor (float32, same affine)
    """
    img = nib.load(str(raw_path))
    data = _sanitize_graddev_matrix(np.asarray(img.dataobj)).astype(np.float32)
    nib.save(nib.Nifti1Image(data, img.affine), str(out_path))

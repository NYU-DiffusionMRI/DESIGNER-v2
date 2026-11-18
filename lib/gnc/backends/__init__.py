"""GNC backend dispatch. Each backend module implements GncBackend."""

from pathlib import Path
from typing import Protocol

import nibabel as nib
import numpy as np

from lib.gnc.backends import gradunwarp, tortoise
from lib.gnc.graddev import GRADDEV_DESCRIPTIONS, sanitize_graddev_matrix
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
        Produce the raw per-voxel gradient Jacobian J in backend axes.

        Called right after compute_warp() with the same arguments; may
        consume files compute_warp() deliberately left on disk.

        Args:
            ref_volume_nii: Same reference volume passed to compute_warp()
            coeff_file: Same coefficient file passed to compute_warp()
            scanner: Same scanner vendor passed to compute_warp()
            out_prefix: Same basename prefix passed to compute_warp()

        Returns:
            Path to a raw (unsanitized, may contain NaN) Jacobian NIfTI with
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
    backend: str,
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
        backend: 'tortoise' or 'gradunwarp'

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
    Export the native gradient Jacobian to {out_prefix}_graddev.nii.gz
    with the backend tag consumed by
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
    _write_graddev(raw, out_path, backend)
    force_remove(raw)
    return out_path


def _write_graddev(raw_path: str | Path, out_path: str | Path, backend: str) -> None:
    """Export native components and identify their convention for TMI."""
    img = nib.load(str(raw_path))
    data = sanitize_graddev_matrix(np.asarray(img.dataobj)).astype(np.float32)
    out = nib.Nifti1Image(data, img.affine)
    out.header['descrip'] = GRADDEV_DESCRIPTIONS[backend]
    nib.save(out, str(out_path))

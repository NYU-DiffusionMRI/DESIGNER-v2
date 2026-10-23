"""Gradient nonlinearity correction (GNC): spatial unwarping, with a choice of
backend ('gradunwarp' or 'tortoise').

Computes a single geometric warp (per scan) from the scanner's gradient
coefficient file and applies it either as a standalone pre-eddy correction
(for topup/eddy to operate on geometrically-correct data) or composed with
eddy's own per-volume displacement fields so the combined GNC+eddy correction
is resampled in a single interpolation pass.

Both backends also produce a per-voxel gradient-deviation ("L") tensor field
(see derive_gnc_graddev), consumed by designer2/tmi.py's -gnc option to
correct the b-matrix itself for gradient nonlinearity, not just image
geometry. Each backend (lib/gnc/backends/) only has to satisfy two contracts
so the rest of this package (and eddy composition in
finalize_single_interpolation) never needs to know which one ran:
compute_warp() returns an FSL relative-warp file, and derive_graddev()
returns a NIfTI with last dimension 9 (a flattened 3x3 matrix per voxel,
L ~= Identity at undistorted voxels).
"""

from lib.gnc.backends import compute_gnc_warp, derive_gnc_graddev
from lib.gnc.pipeline import apply_warp_to_series, finalize_single_interpolation, run_gnc_precorrection
from lib.gnc.request import resolve_gnc_request

__all__ = [
    'apply_warp_to_series',
    'compute_gnc_warp',
    'derive_gnc_graddev',
    'finalize_single_interpolation',
    'resolve_gnc_request',
    'run_gnc_precorrection',
]

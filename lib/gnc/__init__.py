"""Gradient nonlinearity correction (GNC): spatial unwarping, with a choice of
backend ('gradunwarp' or 'tortoise').

Computes a single geometric warp (per scan) from the scanner's gradient
coefficient file and applies it either as a standalone pre-eddy correction
(for topup/eddy to operate on geometrically-correct data) or composed with
eddy's own per-volume displacement fields so the combined GNC+eddy correction
is resampled in a single interpolation pass.

Each backend returns an FSL relative warp from compute_warp() and a native
gradient Jacobian from derive_graddev(). derive_gnc_graddev() exports the
nine row-ordered components and tags the NIfTI with the backend convention.
TMI converts this field to the scanner-frame coil tensor before fitting.
"""

from lib.gnc.backends import compute_gnc_warp, derive_gnc_graddev
from lib.gnc.coeff import strip_coeff_header
from lib.gnc.pipeline import apply_warp_to_series, finalize_single_interpolation, run_gnc_precorrection
from lib.gnc.request import resolve_gnc_request

__all__ = [
    'apply_warp_to_series',
    'compute_gnc_warp',
    'derive_gnc_graddev',
    'finalize_single_interpolation',
    'resolve_gnc_request',
    'run_gnc_precorrection',
    'strip_coeff_header',
]

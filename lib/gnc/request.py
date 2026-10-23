"""Decide whether a requested GNC should run, from the inputs' BIDS sidecars."""

import json
import warnings
from collections.abc import Sequence
from pathlib import Path

from colorama import Fore


def _read_nonlinear_gradient_correction(sidecar: str | Path) -> bool | None:
    """Return the BIDS NonlinearGradientCorrection value, or None if unavailable."""
    if not Path(sidecar).is_file():
        return None
    with open(sidecar) as f:
        return json.load(f).get('NonlinearGradientCorrection')


def resolve_gnc_request(sidecars: Sequence[str | Path]) -> bool:
    """
    Decide whether a requested GNC should run, based on the BIDS
    NonlinearGradientCorrection field of every input that GNC would warp.

    All sidecars must agree: a sidecar that is missing (or lacks the field)
    is assumed to match the others, and if none report it the user's request
    is followed. GNC is skipped if the scanner already applied it.

    Args:
        sidecars: BIDS .json paths of the main DWI series (and the RPE image,
            when it is used by eddy)

    Returns:
        True if GNC should run, False if it should be skipped
    """
    from mrtrix3 import MRtrixError

    states = {str(s): _read_nonlinear_gradient_correction(s) for s in sidecars}
    known = {state for state in states.values() if state is not None}

    if len(known) > 1:
        details = ', '.join(f'{s}: {state}' for s, state in states.items() if state is not None)
        raise MRtrixError(
            'Inputs have inconsistent gradient nonlinearity correction states '
            f'(NonlinearGradientCorrection - {details}); cannot apply -gnc_grad_coeff.')

    if not known:
        warnings.warn(Fore.RED + '[WARNING] NonlinearGradientCorrection not found in BIDS sidecar(s); '
                      'cannot verify whether the scanner already applied gradient nonlinearity '
                      'correction. Proceeding with -gnc_grad_coeff.')
        return True

    already_corrected = known.pop()
    unknown = [s for s, state in states.items() if state is None]
    if unknown:
        warnings.warn(Fore.RED + f'[WARNING] NonlinearGradientCorrection not found for {", ".join(unknown)}; '
                      f'assuming {already_corrected} to match the other input(s).')

    if already_corrected:
        warnings.warn(Fore.RED + '[WARNING] BIDS sidecar(s) report NonlinearGradientCorrection=true '
                      '(already corrected on the scanner); ignoring -gnc_grad_coeff.')
    return not already_corrected

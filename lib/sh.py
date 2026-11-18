"""Standalone rotational invariants using TMI's existing SH estimator."""
import numpy as np

from lib.smi import SMI


def fit_rotational_invariants(dwi, bval, bvec, *, mask=None, sigma=None,
                             lmax=None, beta=None, echo_time=None,
                             rician=False):
    """Return one 4D image per even SH order, plus ordered shell metadata.

    ``lmax`` sets the SH fitting order on every directional diffusion shell;
    None retains TMI's automatic per-shell orders. b-values are in ms/um^2
    and echo times in ms. b=0/isotropic shells retain their arithmetic means
    and have zero higher-order invariants. No SMI priors, noise estimation,
    training, or microstructure regression are run. The shared SH backend
    retains SMI's normalization and rank-one processing across shells.
    """
    if lmax is not None and (not isinstance(lmax, (int, np.integer))
                             or lmax < 0 or lmax % 2):
        raise ValueError("SH lmax must be a nonnegative even integer")
    dwi = np.asarray(dwi)
    if dwi.ndim != 4:
        raise ValueError("SH input must be a 4D DWI")
    if mask is None:
        mask = np.ones(dwi.shape[:3], dtype=bool)
    mask = np.asarray(mask, dtype=bool)
    if mask.shape != dwi.shape[:3] or not np.any(mask):
        raise ValueError("SH mask must match the DWI spatial shape and contain voxels")

    estimator = SMI(np.asarray(bval), np.asarray(bvec), beta=beta,
                    echo_time=echo_time, rician_sh=rician, seed=None)
    automatic = estimator.get_default_lmax()  # Also populates the sorted shell table.
    table = estimator.table_4d
    orders = (np.maximum(automatic, 0).astype(int) if lmax is None
              else np.full(table.shape[1], lmax, dtype=int))
    isotropic = (table[0] <= .05) | (np.abs(table[1]) <= .045)
    orders[isotropic] = 0
    output_lmax = int(np.max(orders)) if lmax is None else lmax

    # Explicit orders must be identifiable; do not silently lower the order.
    for i, order in enumerate(orders):
        ids = ((abs(estimator.b - table[0, i]) < estimator.merge_distance)
               & (abs(estimator.beta - table[1, i]) < estimator.merge_distance)
               & (abs(estimator.echo_time - table[3, i]) < .1))
        n_coefficients = (int(order) + 1) * (int(order) + 2) // 2
        if n_coefficients > np.count_nonzero(ids):
            raise ValueError(f"SH order {order} needs at least {n_coefficients} "
                             f"directions in shell {i}; found {np.count_nonzero(ids)}")
        if not isotropic[i]:
            directions = estimator.dirs
            if directions.shape != (len(estimator.b), 3):
                directions = directions.T
            basis = estimator.get_even_sh(directions[ids], order)
            if np.linalg.matrix_rank(basis) < n_coefficients:
                raise ValueError(f"SH order {order} has a rank-deficient direction "
                                 f"set in shell {i}; reduce -lmax")

    estimator.l_max = orders
    estimator.rotinv_lmax = output_lmax
    invariants = estimator.fit2D4D_LLS_RealSphHarm_wSorting_norm_var(
        dwi.astype(np.float32), mask, sigma=sigma)
    n_shells = table.shape[1]
    images = {}
    for order in range(0, output_lmax + 1, 2):
        start = (order // 2) * n_shells
        images[f"S{order}"] = (invariants[..., start:start + n_shells]
                               if start < invariants.shape[-1]
                               else np.zeros(dwi.shape[:3] + (n_shells,)))
    shells = [dict(volume=i, b_value_s_per_mm2=float(table[0, i] * 1000),
                   b_shape=float(table[1, i]), echo_time_ms=float(table[3, i]),
                   n_measurements=int(table[2, i]), lmax=int(orders[i]))
              for i in range(n_shells)]
    return images, shells

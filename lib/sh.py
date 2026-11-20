"""Standalone rotational invariants using TMI's existing SH estimator."""
import numpy as np

from lib.smi import SMI


def fit_rotational_invariants(dwi, bval, bvec, *, mask=None, sigma=None,
                             lmax=None, beta=None, echo_time=None,
                             rician=False, rician_maxiter=None,
                             diagnostics=None, rician_clip=False, n_jobs=1,
                             rician_clip_mode='maxzero', nonnegative=False):
    """Return one 4D image per even SH order, plus ordered shell metadata.

    ``lmax`` is a common SH fitting order or a sequence with one order per
    sorted shell (ascending b-value, b-shape, then echo time, including b=0).
    None retains TMI's automatic per-shell orders. b-values are in ms/um^2
    and echo times in ms. b=0/isotropic shells use scalar estimation with the
    selected correction and have zero higher-order invariants. No SMI priors, noise estimation,
    training, or microstructure regression are run. The shared SH backend
    retains SMI's normalization and rank-one processing across shells.
    rician_maxiter optionally limits constrained optimization. When supplied,
    diagnostics can receive image-grid iteration counts and early-stop masks.
    rician_clip selects clipped-corrected ordinary LLS instead of Rician
    expectation fitting. rician_clip_mode chooses 'maxzero' or 'abs' for that
    correction; it never changes Rician-expectation initialization.
    n_jobs controls process workers for either correction.
    nonnegative optionally constrains the SH signal; only valid with rician=True.
    """
    requested = None if lmax is None else np.asarray(lmax)
    if requested is not None and (requested.ndim > 1 or requested.size == 0
                                  or requested.dtype.kind not in 'iu'
                                  or np.any(requested < 0)
                                  or np.any(requested % 2)):
        raise ValueError("SH lmax must be a nonnegative even integer or a "
                         "one-dimensional list of nonnegative even integers")
    dwi = np.asarray(dwi)
    if dwi.ndim != 4:
        raise ValueError("SH input must be a 4D DWI")
    if mask is None:
        mask = np.ones(dwi.shape[:3], dtype=bool)
    mask = np.asarray(mask, dtype=bool)
    if mask.shape != dwi.shape[:3] or not np.any(mask):
        raise ValueError("SH mask must match the DWI spatial shape and contain voxels")

    estimator = SMI(np.asarray(bval), np.asarray(bvec), beta=beta,
                    echo_time=echo_time, rician_sh=rician, seed=None,
                    rician_sh_maxiter=rician_maxiter, rician_sh_clip=rician_clip,
                    sh_n_jobs=n_jobs, rician_clip_mode=rician_clip_mode,
                    rician_sh_nonnegative=nonnegative)
    automatic = estimator.get_default_lmax()  # Also populates the sorted shell table.
    table = estimator.table_4d
    if requested is None:
        orders = np.maximum(automatic, 0).astype(int)
    elif requested.ndim == 0:
        orders = np.full(table.shape[1], int(requested), dtype=int)
    else:
        if requested.size != table.shape[1]:
            raise ValueError(f"SH lmax list needs one order per sorted shell "
                             f"(including b=0/isotropic shells): expected "
                             f"{table.shape[1]}, got {requested.size}")
        orders = requested.astype(int, copy=True)
    isotropic = (table[0] <= .05) | (np.abs(table[1]) <= .045)
    orders[isotropic] = 0
    output_lmax = int(np.max(orders if requested is None else requested))

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
                                 f"set in shell {i}; reduce -lmaxSH")

    estimator.l_max = orders
    estimator.rotinv_lmax = output_lmax
    invariants = estimator.fit2D4D_LLS_RealSphHarm_wSorting_norm_var(
        dwi.astype(np.float32), mask, sigma=sigma)
    if diagnostics is not None and estimator.rician_sh_diagnostics is not None:
        from lib.rician_sh import diagnostics_to_images
        diagnostics.update(diagnostics_to_images(estimator.rician_sh_diagnostics, mask))
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

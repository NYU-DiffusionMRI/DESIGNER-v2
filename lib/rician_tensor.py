"""Unconstrained DTI/DKI least squares with a Rician-mean observation model.

This models the retained magnitude noise floor, not the Rician likelihood.
Sigma is the original complex-component noise SD, even after denoising.
"""
import numpy as np
from joblib import Parallel, delayed, effective_n_jobs, parallel_config
from scipy.linalg import pinv, solve_triangular
from scipy.optimize import least_squares

from lib.rician_sh import clipped_magnitude, rician_mean_and_derivative


def tensor_design(grad, kurtosis=False):
    """Design for [log S0, D6, C15], where C = MD**2 W and b is ms/µm²."""
    from lib.tensor import TensorFitting
    g, b = grad[:, :3], grad[:, 3]
    dc, di = TensorFitting.create_tensor_order(None, 2)
    columns = [np.ones((len(b), 1)), -b[:, None] * np.prod(g[:, di], axis=-1) * dc]
    if kurtosis:
        wc, wi = TensorFitting.create_tensor_order(None, 4)
        columns.append(b[:, None]**2 / 6 * np.prod(g[:, wi], axis=-1) * wc)
    return np.column_stack(columns)


def _factor(design):
    if np.linalg.matrix_rank(design) < design.shape[1]:
        raise ValueError("RicianExp tensor fitting requires a full-rank DTI/DKI design")
    q, r = np.linalg.qr(design, mode='reduced')
    return q, r, pinv(design)


def _fit_voxel(y, sigma, design, factors):
    """Start both the first fit and any retry from absolute-value-corrected LLS."""
    q, r, inverse = factors
    initial = inverse @ np.log(np.maximum(clipped_magnitude(y, sigma, mode='abs'), np.finfo(float).eps))
    baseline = np.all(design[:, 1:] == 0, axis=1)
    scale = max(float(np.mean(y[baseline])) if np.any(baseline) else float(np.max(y)),
                float(sigma), np.finfo(float).tiny)
    observed, noise = y / scale, sigma / scale
    # A voxelwise constant changes objective units, not the minimizer.
    residual_scale = noise if noise > 0 else 1.0
    start = initial.copy()
    start[0] -= np.log(scale)
    start = r @ start

    def evaluate(z, jac=False):
        logs = q @ z
        signal = np.exp(np.clip(logs, -700, 10))
        ds = signal.copy()
        # Numerical continuation for extreme clipped starts, not a tensor bound.
        # Accepted fits are checked to be in the exact exponential region below.
        high = logs > 10
        signal[high] *= 1 + logs[high] - 10
        ds[logs < -700] = 0
        mean, derivative = rician_mean_and_derivative(signal, noise)
        if jac:
            return (derivative * ds / residual_scale)[:, None] * q
        return (mean - observed) / residual_scale

    initial_cost = float(np.sum(evaluate(start)**2))
    results = []
    evaluations = 0
    for budget, tol in [(200, 1e-8), (1000, 1e-10)]:
        result = least_squares(evaluate, start, jac=lambda z: evaluate(z, True),
                               method='trf', max_nfev=budget, ftol=tol, xtol=tol, gtol=tol)
        evaluations += result.nfev
        results.append(result)
        if result.success:
            break
    # A retry must never replace a better first solution.
    result = min(results, key=lambda fit: fit.cost)
    logs = q @ result.x
    if np.any(logs > 10) or not np.all(np.isfinite(result.x)):
        raise RuntimeError("RicianExp tensor fit did not reach a finite exponential-signal solution")
    theta = solve_triangular(r, result.x)
    theta[0] += np.log(scale)
    diagnostics = dict(status=int(result.status), nfev=int(evaluations),
                       cost=float(2 * result.cost), initial_cost=initial_cost,
                       retried=len(results) > 1, nonconverged=not result.success,
                       signal_underflow=bool(np.any(logs < -700)))
    return theta, diagnostics


def _fit_block(data, sigma, grad, kurtosis, design, factors, local_fields):
    from lib.tensor import TensorFitting
    fits, diagnostics = [], []
    for i, y in enumerate(data):
        x, qr = design, factors
        if local_fields is not None:
            local_grad = TensorFitting._voxel_grad_from_L(None, grad, local_fields[i])
            x = tensor_design(local_grad, kurtosis)
            qr = _factor(x)
        theta, info = _fit_voxel(y, sigma[i], x, qr)
        fits.append(theta)
        diagnostics.append(info)
    return np.asarray(fits), diagnostics


def fit_rician_tensor(dwi, mask, grad, sigma, *, kurtosis=False, L=None, n_jobs=1):
    """Return existing tmi tensor/S0/design conventions plus voxel diagnostics.

    No constraints on diffusion eigenvalues, kurtosis, or signal monotonicity.
    Signal amplitudes are exponential by definition of the DTI/DKI model.
    The input observations, gradients, sigma and mask are never modified.
    """
    dwi = np.asarray(dwi, dtype=float)
    mask = np.asarray(mask, dtype=bool)
    sigma = np.asarray(sigma, dtype=float)
    if dwi.ndim != 4 or mask.shape != dwi.shape[:3] or sigma.shape != mask.shape:
        raise ValueError("RicianExp requires DWI, mask and sigma on the same 3D spatial grid")
    if not np.any(mask):
        raise ValueError("RicianExp tensor mask is empty")
    data, noise = dwi[mask], sigma[mask]
    if not np.all(np.isfinite(data)) or not np.all(np.isfinite(noise)) or np.any(noise < 0):
        raise ValueError("RicianExp requires finite observations and finite nonnegative sigma inside the mask")
    grad = np.array(grad, dtype=float, copy=True)
    if grad.shape != (dwi.shape[-1], 4) or not np.all(np.isfinite(grad)) or np.any(grad[:, 3] < 0):
        raise ValueError("RicianExp requires a finite N-by-4 gradient table with nonnegative b values")
    if np.max(grad[:, 3]) >= 99:
        grad[:, 3] /= 1000
    norms = np.linalg.norm(grad[:, :3], axis=1)
    grad[:, :3] /= np.where(norms > 0, norms, 1)[:, None]
    design = tensor_design(grad, kurtosis)
    fields = None
    if L is not None:
        L = np.asarray(L, dtype=float)
        if L.shape != mask.shape + (3, 3) or not np.all(np.isfinite(L[mask])):
            raise ValueError("RicianExp gradient field must be finite with shape (X,Y,Z,3,3)")
        fields = L[mask]
    factors = _factor(design) if fields is None else None
    workers = min(effective_n_jobs(n_jobs), max(1, (len(data) + 31) // 32))
    if workers == 1:
        blocks = [_fit_block(data, noise, grad, kurtosis, design, factors, fields)]
    else:
        block_size = max(32, min(256, (len(data) + 4 * workers - 1) // (4 * workers)))
        with parallel_config(backend='loky', inner_max_num_threads=1):
            blocks = Parallel(n_jobs=workers)(
                delayed(_fit_block)(data[start:start+block_size], noise[start:start+block_size],
                    grad, kurtosis, design, factors,
                    None if fields is None else fields[start:start+block_size])
                for start in range(0, len(data), block_size))
    theta = np.concatenate([block[0] for block in blocks])
    infos = [info for block in blocks for info in block[1]]
    diagnostics = {key: np.array([info[key] for info in infos]) for key in infos[0]}
    matrices = theta[:, 1:7][:, [0, 1, 2, 1, 3, 4, 2, 4, 5]].reshape(-1, 3, 3)
    diagnostics['nonpositive_D'] = np.linalg.eigvalsh(matrices)[:, 0] <= 0
    dt = theta[:, 1:].T.copy()
    if kurtosis:
        md = dt[[0, 3, 5]].mean(axis=0)
        with np.errstate(divide='ignore', invalid='ignore'):
            dt[6:] /= md[None, :]**2
    return dt, np.exp(theta[:, 0]), design if fields is None else None, diagnostics

"""Rician-expectation SH fitting, adapted from S. Coelho's SHestimation.

The underlying signal is Y @ coefficients; the observations are magnitude data.
Fit the Rician *mean* by nonlinear least squares, not the Rician likelihood.
"""
import warnings

import numpy as np
from scipy.special import i0e, i1e
from scipy.optimize import least_squares


def rician_mean_and_derivative(signal, sigma):
    """Return E[magnitude | signal, sigma] and its analytic signal derivative.

    Sigma is the standard deviation of each complex Gaussian noise component.
    Scaled Bessel functions avoid exponential overflow. At SNR > 1e4 use
    the asymptote (relative error O(SNR**-4)) to avoid squaring huge SNRs.
    """
    signal, sigma = np.broadcast_arrays(np.asarray(signal, float), np.asarray(sigma, float))
    amplitude = np.abs(signal)
    mean = np.array(amplitude, copy=True)
    deriv = np.array(np.sign(signal), copy=True)
    regular = (sigma > 0) & (amplitude <= 1e4 * sigma)
    ratio = signal[regular] / sigma[regular]
    z = (ratio / 2)**2
    bessel_sum = i0e(z) + i1e(z)
    c = np.sqrt(np.pi / 2)
    mean[regular] = sigma[regular] * c * (i0e(z) + 2 * z * bessel_sum)
    deriv[regular] = c * ratio * bessel_sum / 2
    high = (sigma > 0) & ~regular
    inv_snr = sigma[high] / amplitude[high]
    mean[high] = amplitude[high] + sigma[high] * inv_snr / 2
    deriv[high] *= 1 - inv_snr**2 / 2
    return mean, deriv


def fit_rician_sh(basis, measurements, sigma, initial, chunk_size=256):
    """Refine an LLS fit with bounded-memory batched Levenberg--Marquardt.

    Shapes: basis (directions, coefficients), measurements (directions, voxels),
    sigma (voxels,), initial (coefficients, voxels). Zero-noise voxels retain
    their LLS coefficients exactly. Fits are scaled per voxel so tolerances do
    not depend on signal units. No spatial smoothing or intensity correction
    is applied to the input data.
    """
    result = np.array(initial, dtype=float, copy=True)
    indices = np.flatnonzero(sigma > 0)
    unfinished = 0
    for start in range(0, indices.size, chunk_size):
        ids = indices[start:start + chunk_size]
        fitted, remaining = _fit_batch(
            basis, measurements[:, ids].T, sigma[ids], result[:, ids].T)
        result[:, ids] = fitted.T
        unfinished += remaining
    if unfinished:
        warnings.warn(
            f"Rician SH fit reached the iteration limit in {unfinished} voxel(s) "
            "in this shell; returning the best accepted coefficients.",
            RuntimeWarning, stacklevel=2)
    return result


def _fit_batch(basis, measurements, sigma, initial, max_iter=200):
    # This follows SHestimation's batched LM: accept only cost-decreasing
    # updates, adjust damping per voxel, and drop converged voxels from work.
    scale = np.maximum(np.max(np.abs(measurements), axis=1), sigma)
    measurements = measurements / scale[:, None]
    sigma = sigma[:, None] / scale[:, None]
    coefficients = initial / scale[:, None]
    nvox, ncoeff = coefficients.shape
    active = np.ones(nvox, dtype=bool)
    damping = np.full(nvox, 1e-3)
    stall = np.zeros(nvox, dtype=int)
    eye = np.eye(ncoeff)

    def pieces(x, ids):
        mean, derivative = rician_mean_and_derivative(x @ basis.T, sigma[ids])
        residual = mean - measurements[ids]
        cost = np.einsum('vi,vi->v', residual, residual) / 2
        return residual, derivative, cost

    for _ in range(max_iter):
        ids = np.flatnonzero(active)
        if not ids.size:
            break
        x = coefficients[ids]
        residual, derivative, cost = pieces(x, ids)
        jacobian = derivative[:, :, None] * basis[None, :, :]
        jtj = np.einsum('vki,vkj->vij', jacobian, jacobian)
        jtr = np.einsum('vki,vk->vi', jacobian, residual)
        diagonal = np.maximum(np.einsum('vii->vi', jtj), np.finfo(float).eps)
        system = jtj + damping[ids, None, None] * diagonal[:, :, None] * eye
        try:
            step = np.linalg.solve(system, -jtr[:, :, None])[:, :, 0]
        except np.linalg.LinAlgError:
            step = np.stack([np.linalg.lstsq(a, -b, rcond=None)[0]
                             for a, b in zip(system, jtr)])
        candidate = x + step
        _, _, new_cost = pieces(candidate, ids)
        new_cost = np.where(np.isfinite(new_cost), new_cost, np.inf)
        improved = new_cost < cost
        coefficients[ids[improved]] = candidate[improved]
        damping[ids[improved]] = np.maximum(damping[ids[improved]] * 0.3, 1e-9)
        damping[ids[~improved]] = np.minimum(damping[ids[~improved]] * 3, 1e9)
        tiny_step = np.linalg.norm(step, axis=1) <= 1e-8 * (1 + np.linalg.norm(candidate, axis=1))
        tiny_cost = np.abs(cost - new_cost) <= 1e-12 * np.maximum(cost, 1e-30)
        # Rejected steps must not stop a fit while it is searching for damping.
        converged = (improved & (tiny_step | tiny_cost)) | (
            ~improved & tiny_step & (damping[ids] >= 1e9))
        stall[ids] = np.where(converged, stall[ids] + 1, 0)
        active[ids[stall[ids] >= 2]] = False
    # Difficult low-SNR voxels can need more steps than a shared batch. Finish
    # just those voxels with SciPy's independently controlled trust region.
    for voxel in np.flatnonzero(active):
        def residual(x):
            return rician_mean_and_derivative(basis @ x, sigma[voxel, 0])[0] - measurements[voxel]

        def jacobian(x):
            derivative = rician_mean_and_derivative(basis @ x, sigma[voxel, 0])[1]
            return derivative[:, None] * basis

        fit = least_squares(residual, coefficients[voxel], jac=jacobian,
                            method='lm' if basis.shape[0] >= ncoeff else 'trf',
                            max_nfev=2000, ftol=1e-10, xtol=1e-10, gtol=1e-10)
        previous = residual(coefficients[voxel])
        if np.isfinite(fit.cost) and fit.cost <= previous @ previous / 2:
            coefficients[voxel] = fit.x
        active[voxel] = not fit.success
    return coefficients * scale[:, None], int(active.sum())

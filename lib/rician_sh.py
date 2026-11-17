"""Rician-expectation SH fitting, adapted from S. Coelho's SHestimation.

The underlying signal is Y @ coefficients; the observations are magnitude data.
Fit the Rician *mean* by nonlinear least squares, not the Rician likelihood.
"""
import warnings

import numpy as np
from scipy.special import i0e, i1e
from scipy.optimize import least_squares, minimize, LinearConstraint, brentq


def nonnegative_directions(count=65536):
    """Deterministic sphere grid used to constrain the underlying SH signal.

    A finite grid is not a continuous-sphere nonnegativity certificate.
    Acquired directions are included separately by the constrained fitter.
    """
    z = 1 - 2 * (np.arange(count) + .5) / count
    phi = np.arange(count) * np.pi * (3 - np.sqrt(5))
    radius = np.sqrt(1 - z*z)
    return np.column_stack((radius*np.cos(phi), radius*np.sin(phi), z))


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


def fit_rician_sh(basis, measurements, sigma, initial, chunk_size=256,
                  constraint_basis=None):
    """Refine an LLS fit with bounded-memory batched Levenberg--Marquardt.

    Shapes: basis (directions, coefficients), measurements (directions, voxels),
    sigma (voxels,), initial (coefficients, voxels). With constraint_basis,
    require nonnegative signal at its rows and all acquired directions using
    linear constraints in coefficient space. The unconstrained fit is used
    when already feasible; otherwise solve constrained mean least squares.
    Zero-noise fits reduce to constrained LLS. Without constraints, zero-noise
    voxels retain their LLS coefficients exactly. Fits are scaled so tolerances do
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
    if constraint_basis is not None:
        result = _enforce_nonnegative(basis, measurements, sigma, initial,
                                      result, constraint_basis, chunk_size)
    return result


def _enforce_nonnegative(basis, measurements, sigma, initial, fitted,
                         constraint_basis, chunk_size):
    """Constraint generation avoids passing 65k inequalities to every QP.

    Each subproblem contains all acquired directions, an initial sphere
    subset, and its worst violated grid constraints. Stop only after checking
    the entire supplied grid. Unfinished fits return the best feasible point
    with a warning, never an unconstrained fallback.
    """
    constraint_basis = np.asarray(constraint_basis, dtype=float)
    if (constraint_basis.ndim != 2 or constraint_basis.shape[1] != basis.shape[1]
            or not np.all(np.isfinite(constraint_basis))):
        raise ValueError('constraint_basis must be finite with one column per SH coefficient')
    grid = np.vstack((basis, constraint_basis))
    # Feasibility restoration adjusts only the constant SH coefficient.
    if not np.allclose(grid[:, 0], grid[0, 0], rtol=1e-12, atol=0) or grid[0, 0] <= 0:
        raise ValueError('Nonnegative SH fitting requires a positive constant first basis column')
    constant = grid[0, 0]
    seed = np.unique(np.r_[np.arange(len(basis)),
                           np.linspace(len(basis), len(grid)-1,
                                       min(256, len(constraint_basis)), dtype=int)])
    unfinished = 0
    for start in range(0, fitted.shape[1], chunk_size):
        stop = min(start + chunk_size, fitted.shape[1])
        minima = np.min(grid @ fitted[:, start:stop], axis=0)
        for voxel in start + np.flatnonzero(minima < 0):
            scale = max(np.max(np.abs(measurements[:, voxel])), sigma[voxel], np.finfo(float).tiny)
            observed = measurements[:, voxel] / scale
            noise = sigma[voxel] / scale

            def objective(coef):
                signal = basis @ coef
                if noise == 0:
                    # On the feasible domain, the zero-noise model is linear.
                    mean, derivative = signal, np.ones_like(signal)
                else:
                    mean, derivative = rician_mean_and_derivative(signal, noise)
                residual = mean - observed
                return .5*(residual @ residual), basis.T @ (residual*derivative)

            def feasible(coef, minimum=None):
                coef = np.array(coef, copy=True)
                if minimum is None:
                    minimum = np.min(grid @ coef)
                if minimum < 0:
                    coef[0] += (-minimum + 32*np.finfo(float).eps) / constant
                return coef

            # Start from feasible LLS, not a folded negative-signal solution.
            best = feasible(initial[:, voxel] / scale)
            best_cost = objective(best)[0]
            active = seed.copy()
            current = best.copy()
            converged = False
            restarted_zero = False
            for _ in range(24):
                fit = minimize(objective, current, jac=True, method='SLSQP',
                               constraints=[LinearConstraint(grid[active], 0, np.inf)],
                               options={'maxiter':500, 'ftol':1e-11})
                if not np.all(np.isfinite(fit.x)):
                    break
                signal = grid @ fit.x
                minimum = signal.min()
                candidate = feasible(fit.x, minimum)
                cost = objective(candidate)[0]
                if np.isfinite(cost) and cost < best_cost:
                    best, best_cost = candidate, cost
                if fit.success and minimum >= -1e-10:
                    if noise > 0 and best[0]*constant < 1e-8 and not restarted_zero:
                        # mu'(0)=0 makes the all-zero signal a stationary point
                        # even when it is not a minimum. Restart from a positive
                        # isotropic signal before accepting that boundary point.
                        mean_observed = observed.mean()
                        floor = noise*np.sqrt(np.pi/2)
                        amplitude = (.1*noise if mean_observed <= floor else brentq(
                            lambda value: rician_mean_and_derivative(value, noise)[0]-mean_observed,
                            0., mean_observed, xtol=1e-12))
                        current = np.zeros_like(best)
                        current[0] = amplitude/constant
                        cost = objective(current)[0]
                        if cost < best_cost:
                            best, best_cost = current.copy(), cost
                        restarted_zero = True
                        continue
                    converged = True
                    break
                # Keep the most violated rows. A failed SLSQP solve without
                # violations gets another attempt from its feasible result.
                violated = np.flatnonzero(signal < -1e-10)
                if violated.size:
                    worst = violated[np.argsort(signal[violated])[:32]]
                    active = np.union1d(active, worst)
                current = candidate
            fitted[:, voxel] = best * scale
            unfinished += not converged
    if unfinished:
        warnings.warn(f'Nonnegative Rician SH fit did not converge in {unfinished} voxel(s) '
                      'in this shell; returning the best grid-feasible coefficients.',
                      RuntimeWarning, stacklevel=2)
    return fitted


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

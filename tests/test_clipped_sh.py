"""Clipped-magnitude constrained LLS has a convex, independent reference."""
import numpy as np
import pytest
from scipy.linalg import lstsq
from scipy.optimize import minimize, LinearConstraint

from lib.rician_sh import fit_rician_sh, nonnegative_directions
from lib.smi import SMI


@pytest.mark.parametrize('scale', [1e-6, 1., 1e6])
def test_clipped_constrained_known_solution_without_nonlinear_search(monkeypatch, scale):
    import lib.rician_sh as module
    def unexpected(*args, **kwargs):
        raise AssertionError('Clipped LLS must not evaluate the Rician expectation')
    monkeypatch.setattr(module, '_fit_batch', unexpected)
    monkeypatch.setattr(module, 'rician_mean_and_derivative', unexpected)
    basis = np.array([[1., -1.], [1., 0.], [1., 1.]])
    grid = np.array([[1., -2.], [1., 2.]])
    raw = np.sqrt(np.array([[0.], [1.], [4.]]) + .25)*scale
    fitted = fit_rician_sh(basis, raw, np.array([.5])*scale,
                           constraint_basis=grid, clipped=True)/scale
    np.testing.assert_allclose(fitted[:, 0], [8/7, 4/7], atol=1e-7)
    assert np.min(grid @ fitted) >= -1e-12


@pytest.mark.parametrize('mode', ['maxzero', 'abs'])
def test_clipped_fit_matches_independent_quadratic_program(mode):
    rng = np.random.default_rng(318)
    dirs = nonnegative_directions(40)
    smi = SMI(np.ones(40), dirs)
    basis = smi.get_even_sh(dirs, 4)
    grid = np.vstack([basis, smi.get_even_sh(nonnegative_directions(512), 4)])
    raw = np.hypot(.2 + .4*rng.normal(size=(40, 3)), .4*rng.normal(size=(40, 3)))
    if mode == 'abs':
        raw[0, :] = 4.  # Force an infeasible SH projection even with the raised low-signal floor.
    difference = raw**2 - .4**2
    corrected = np.sqrt(abs(difference) if mode == 'abs' else np.maximum(0, difference))
    initial = lstsq(basis, corrected)[0]
    assert np.min(grid @ initial) < -.01
    fitted = fit_rician_sh(basis, raw, np.full(3, .4), constraint_basis=grid, clipped=True, clip_mode=mode)
    for v in range(3):
        def objective(c):
            residual = basis @ c - corrected[:, v]
            return .5*(residual @ residual), basis.T @ residual
        reference = minimize(objective, np.r_[1., np.zeros(14)], jac=True,
                             method='SLSQP', constraints=[LinearConstraint(grid, 0, np.inf)],
                             options={'ftol':1e-12, 'maxiter':1000})
        assert reference.success
        np.testing.assert_allclose(fitted[:, v], reference.x, atol=3e-6)
    assert np.min(grid @ fitted) >= -1e-12


def test_clipped_feasible_fit_is_exact_lls_and_zero_data_stays_zero():
    basis = np.array([[1., -1.], [1., 0.], [1., 1.]])
    raw = np.array([[2., .2], [2.2, .3], [2.4, .4]])
    noise = np.array([.4, .5])
    corrected = np.sqrt(np.maximum(0, raw**2-noise**2))
    expected = lstsq(basis, corrected)[0]
    fitted = fit_rician_sh(basis, raw, noise, constraint_basis=basis, clipped=True)
    np.testing.assert_array_equal(fitted, expected)
    np.testing.assert_array_equal(fitted[:, 1], 0)


def test_smi_rician_precedence_and_rejection_of_noise_bias():
    estimator = SMI(np.ones(6), nonnegative_directions(6), rician_sh=True, rician_sh_clip=True)
    assert estimator.rician_sh and not estimator.rician_sh_clip
    with pytest.raises(ValueError, match='noise_bias'):
        SMI(np.ones(6), nonnegative_directions(6), rician_sh_clip=True, noise_bias=True)


def test_correction_formulas_broadcast_and_preserve_input():
    from lib.rician_sh import clipped_magnitude
    raw = np.array([[0., 1., 2.], [0., 2., 4.]])
    sigma = np.array([[1.], [2.]])
    original = raw.copy()
    np.testing.assert_allclose(clipped_magnitude(raw, sigma), [[0, 0, np.sqrt(3)], [0, 0, np.sqrt(12)]])
    np.testing.assert_allclose(clipped_magnitude(raw, sigma, mode='abs'), [[1, 0, np.sqrt(3)], [2, 0, np.sqrt(12)]])
    np.testing.assert_array_equal(raw, original)
    with pytest.raises(ValueError, match='mode'):
        clipped_magnitude(raw, sigma, mode='invalid')


def test_absolute_correction_is_forwarded_to_parallel_workers():
    rng = np.random.default_rng(813)
    basis = np.array([[1., -1.], [1., 0.], [1., 1.]])
    grid = np.array([[1., -2.], [1., 2.]])
    raw = abs(rng.normal(size=(3, 65)))
    noise = np.full(65, .5)
    serial = fit_rician_sh(basis, raw, noise, constraint_basis=grid, clipped=True, clip_mode='abs')
    parallel = fit_rician_sh(basis, raw, noise, constraint_basis=grid, clipped=True, clip_mode='abs', n_jobs=2)
    np.testing.assert_allclose(parallel, serial, atol=1e-10)
    assert np.min(grid @ parallel) >= -1e-12


def test_expectation_initialization_ignores_absolute_correction_setting(monkeypatch):
    import lib.rician_sh as module
    basis = np.array([[1., -1.], [1., 0.], [1., 1.]])
    raw = np.array([[.1], [.2], [1.]])
    noise = np.array([.5])
    expected = lstsq(basis, np.sqrt(np.maximum(0, raw**2 - noise**2)))[0]
    def capture(basis, measured, sigma, initial):
        np.testing.assert_array_equal(initial.T, expected)
        return initial.copy(), 0
    monkeypatch.setattr(module, '_fit_batch', capture)
    fit_rician_sh(basis, raw, noise, clipped=False, clip_mode='abs')

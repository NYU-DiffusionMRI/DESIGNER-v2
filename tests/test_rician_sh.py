"""Numerical and SMI integration checks; no external imaging data required."""
import numpy as np
import pytest
from scipy.optimize import least_squares
from scipy.stats import rice

from lib.rician_sh import fit_rician_sh, rician_mean_and_derivative, nonnegative_directions
from lib.smi import SMI


def protocol():
    rng = np.random.default_rng(7)
    dirs = rng.normal(size=(90, 3))
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    return np.repeat([0., 1., 2.], [6, 90, 90]), np.vstack((np.zeros((6, 3)), dirs, dirs))


def model(rician=False, **kwargs):
    b, g = protocol()
    obj = SMI(b, g, rician_sh=rician, compartments=['IAS', 'EAS'], **kwargs)
    obj.set_lmax_rotinv_smi()
    return obj


def test_mean_and_analytic_jacobian():
    value, derivative = rician_mean_and_derivative(1., 1.)
    np.testing.assert_allclose(value, rice.mean(1.), rtol=1e-14)
    assert 0 < derivative < 1
    s = np.r_[-30., -3., -.1, 0., np.geomspace(.001, 30, 50)]
    value, derivative = rician_mean_and_derivative(s, 1.)
    np.testing.assert_allclose(value, rice.mean(np.abs(s)), rtol=1e-13)
    h = 1e-5
    numerical = (rician_mean_and_derivative(s + h, 1.)[0] -
                 rician_mean_and_derivative(s - h, 1.)[0]) / (2*h)
    np.testing.assert_allclose(derivative, numerical, atol=1e-9)
    value, derivative = rician_mean_and_derivative(s, 1e-100)
    np.testing.assert_allclose(value, np.abs(s), atol=1e-99)
    np.testing.assert_allclose(derivative, np.sign(s), atol=1e-10)
    np.testing.assert_array_equal(rician_mean_and_derivative(s, 0)[0], np.abs(s))


@pytest.mark.parametrize('scale', [1e-6, 1., 1e6])
def test_fit_matches_independent_scipy_solver_and_rescaling(scale):
    rng = np.random.default_rng(17)
    obj = model()
    y = obj.get_even_sh(obj.dirs[6:96], 4)
    true = rng.normal(0, .03, (15, 7))
    true[0] = 2
    sigma = np.linspace(.05, .3, 7)
    m = np.hypot(y @ true + rng.normal(size=(90,7))*sigma,
                 rng.normal(size=(90,7))*sigma)
    initial = np.linalg.lstsq(y, m, rcond=None)[0]
    fitted = fit_rician_sh(y, m*scale, sigma*scale, initial*scale, chunk_size=3)/scale
    reference = np.column_stack([
        least_squares(lambda c: rician_mean_and_derivative(y@c, sigma[i])[0]-m[:,i],
                      initial[:,i], method='lm', ftol=1e-12, xtol=1e-12, gtol=1e-12).x
        for i in range(7)])
    np.testing.assert_allclose(fitted, reference, atol=2e-5, rtol=2e-4)


def test_known_expectation_recovery_and_zero_noise():
    obj = model()
    y = obj.get_even_sh(obj.dirs[6:96], 4)
    coefficients = np.zeros((15, 4)); coefficients[0] = 2; coefficients[3] = .2
    m = rician_mean_and_derivative(y @ coefficients, .3)[0]
    initial = np.linalg.lstsq(y, m, rcond=None)[0]
    fitted = fit_rician_sh(y, m, np.full(4, .3), initial)
    np.testing.assert_allclose(fitted, coefficients, atol=1e-7)
    assert np.linalg.norm(initial-coefficients) > .1
    np.testing.assert_array_equal(fit_rician_sh(y, m, np.zeros(4), initial), initial)
    np.testing.assert_allclose(fit_rician_sh(y, m, np.full(4, 1e-12), initial), initial, atol=1e-10)


@pytest.mark.parametrize('shape', [(2, 2, 2), (2, 2, 1)])
def test_smi_shell_invariants_mask_and_unchanged_zero_noise(shape):
    b, g = protocol()
    # l=0 + l=2 signal lies exactly within both shells' SH spaces.
    signal = np.r_[np.ones(6), .7+.1*g[6:96,2]**2, .4+.1*g[96:,2]**2]
    dwi = np.broadcast_to(signal, shape+(b.size,)).copy()
    mask = np.ones(shape, bool); mask[0,0,0] = False
    sigma = np.full(shape, .2); sigma[~mask] = np.nan
    measured = rician_mean_and_derivative(dwi, .2)[0]
    measured[~mask] = np.nan
    outputs = []
    for use_rician in (False, True):
        obj = model(use_rician); obj.mask = mask
        outputs.append(obj.fit2D4D_LLS_RealSphHarm_wSorting_norm_var(
            measured, mask, rank1sl=False, sigma=sigma))
    # Reference implementation deliberately preserves b0 shell means.
    np.testing.assert_array_equal(outputs[0][...,0], outputs[1][...,0])
    clean = model(); clean.mask=mask
    truth = clean.fit2D4D_LLS_RealSphHarm_wSorting_norm_var(dwi, mask, rank1sl=False)
    np.testing.assert_allclose(outputs[1][...,1:3][mask], truth[...,1:3][mask], atol=1e-7)
    assert np.all(outputs[1][~mask] == 0)
    obj = model(True); obj.mask=mask
    zero = obj.fit2D4D_LLS_RealSphHarm_wSorting_norm_var(
        measured, mask, rank1sl=False, sigma=np.zeros(shape))
    np.testing.assert_array_equal(zero, outputs[0])


@pytest.mark.parametrize('sigma, message', [
    (None, 'explicit sigma'), (np.ones((3,3,3)), 'shape'),
    (np.full((2,2,2), -1.), 'nonnegative'),
    (np.full((2,2,2), np.nan), 'finite'),
    (np.full((2,2,2), np.inf), 'finite')])
def test_invalid_sigma(sigma, message):
    obj = model(True); obj.mask = np.ones((2,2,2), bool)
    dwi = np.ones((2,2,2,186))
    with pytest.raises(ValueError, match=message):
        obj.fit2D4D_LLS_RealSphHarm_wSorting_norm_var(dwi, obj.mask, sigma=sigma)


def test_no_implicit_sigma_or_double_correction():
    obj = model(True)
    with pytest.raises(ValueError, match='explicit sigma'):
        obj.fit(np.ones((2,2,2,186)))
    with pytest.raises(ValueError, match='combined'):
        model(True, noise_bias='rician')


@pytest.mark.parametrize('value', [-1., np.nan, np.inf])
def test_invalid_magnitude_data(value):
    obj=model(True);obj.mask=np.ones((2,2,2),bool)
    dwi=np.ones((2,2,2,186));dwi[0,0,0,8]=value
    with pytest.raises(ValueError, match='magnitude DWI'):
        obj.fit2D4D_LLS_RealSphHarm_wSorting_norm_var(dwi,obj.mask,sigma=np.ones((2,2,2)))


def test_isotropic_shells_and_multiple_echo_times():
    b,g=protocol()
    b=np.tile(b,2);g=np.tile(g,(2,1))
    te=np.repeat([60.,90.],186)
    beta=np.ones_like(b);beta[96:186]=0;beta[282:]=0
    clean=np.broadcast_to((np.exp(-b*.7-te/100))[None,None,None,:],(2,2,2,b.size))
    data=rician_mean_and_derivative(clean,.2)[0]
    outputs=[]
    for use_rician in (False,True):
        obj=SMI(b,g,beta=beta.copy(),echo_time=te,rician_sh=use_rician)
        obj.set_lmax_rotinv_smi();obj.mask=np.ones((2,2,2),bool)
        outputs.append(obj.fit2D4D_LLS_RealSphHarm_wSorting_norm_var(
            data,obj.mask,rank1sl=False,sigma=np.full((2,2,2),.2)))
    for shell in range(obj.table_4d.shape[1]):
        if obj.table_4d[0,shell] == 0 or obj.table_4d[1,shell] == 0:
            np.testing.assert_array_equal(outputs[0][...,shell],outputs[1][...,shell])
        else:
            expected=np.exp(-obj.table_4d[0,shell]*.7-obj.table_4d[3,shell]/100)
            np.testing.assert_allclose(outputs[1][...,shell],expected,atol=1e-7)



def test_default_float32_data_keep_double_precision_shell_means():
    # Preserve the original vectorize() float64 accumulation, even when the
    # DWI loader supplies float32. Float32 means perturb the downstream fit.
    obj=model();obj.mask=np.ones((2,2,2),bool)
    rng=np.random.default_rng(29)
    dwi=rng.uniform(100,1000,(2,2,2,186)).astype(np.float32)
    actual=obj.fit2D4D_LLS_RealSphHarm_wSorting_norm_var(dwi,obj.mask)
    expected=dwi[...,:6].astype(float).mean(axis=-1)
    np.testing.assert_allclose(actual[...,0],expected,rtol=1e-15)


@pytest.mark.parametrize('scale', [1e-6, 1., 1e6])
def test_nonnegative_zero_noise_has_known_constrained_linear_solution(scale):
    # Minimize 3*(c0-1)^2 + 2*(c1-1)^2 subject to c0 >= 2*abs(c1).
    # The solution is c0=8/7, c1=4/7; ordinary LLS is [1,1].
    basis = np.array([[1., -1.], [1., 0.], [1., 1.]])
    constraints = np.array([[1., -2.], [1., 2.]])
    measured = np.array([[0.], [1.], [2.]])*scale
    initial = np.array([[1.], [1.]])*scale
    result = fit_rician_sh(basis, measured, np.zeros(1), initial,
                           constraint_basis=constraints)/scale
    np.testing.assert_allclose(result[:,0], [8/7, 4/7], atol=1e-7)
    assert np.min(constraints @ result) >= -1e-12
    # Default unconstrained zero-noise behavior is preserved for diagnostics.
    np.testing.assert_array_equal(fit_rician_sh(basis, measured, np.zeros(1), initial), initial)


def test_nonnegative_low_snr_signal_on_acquired_and_dense_grid():
    rng = np.random.default_rng(20261120)
    obj = model()
    basis = obj.get_even_sh(obj.dirs[6:66], 4)
    grid = obj.get_even_sh(nonnegative_directions(4096), 4)
    sigma = np.full(12, .4)
    measured = np.hypot(.3 + .4*rng.normal(size=(60,12)),
                        .4*rng.normal(size=(60,12)))
    initial = np.linalg.lstsq(basis, measured, rcond=None)[0]
    old = fit_rician_sh(basis, measured, sigma, initial)
    assert (grid @ old).min() < -.1
    result = fit_rician_sh(basis, measured, sigma, initial,
                           constraint_basis=grid, chunk_size=3)
    assert (grid @ result).min() >= -1e-10
    assert (basis @ result).min() >= -1e-10
    assert np.isfinite(result).all()


def test_nonnegative_constraint_does_not_change_feasible_fit():
    obj = model()
    basis = obj.get_even_sh(obj.dirs[6:96], 4)
    grid = obj.get_even_sh(nonnegative_directions(1024), 4)
    true = np.zeros((15,2)); true[0] = 3; true[2] = .1
    measured = rician_mean_and_derivative(basis @ true, .1)[0]
    initial = np.linalg.lstsq(basis, measured, rcond=None)[0]
    old = fit_rician_sh(basis, measured, np.full(2,.1), initial)
    result = fit_rician_sh(basis, measured, np.full(2,.1), initial, constraint_basis=grid)
    np.testing.assert_array_equal(result, old)
    np.testing.assert_allclose(result, true, atol=1e-7)


def test_nonnegative_failed_solver_returns_feasible_result_with_warning(monkeypatch):
    from types import SimpleNamespace
    import lib.rician_sh as module
    basis = np.array([[1.,-1.],[1.,0.],[1.,1.]])
    grid = np.array([[1.,-2.],[1.,2.]])
    monkeypatch.setattr(module, 'minimize', lambda fun, x, **kwargs:
                        SimpleNamespace(x=np.array([1.,1.]), success=False))
    with pytest.warns(RuntimeWarning, match='grid-feasible'):
        result = fit_rician_sh(basis, np.array([[0.],[1.],[2.]]), np.zeros(1),
                               np.array([[1.],[1.]]), constraint_basis=grid)
    assert (grid @ result).min() >= -1e-12


def test_smi_uses_nonnegative_signal_constraints_in_forward_projection():
    obj = model(True)
    obj.mask = np.ones((2,2,1), bool)
    obj.project_forward_SH = True
    rng = np.random.default_rng(31)
    data = np.hypot(.3 + .4*rng.normal(size=(2,2,1,186)),
                    .4*rng.normal(size=(2,2,1,186)))
    basis = obj.get_even_sh(obj.dirs[96:],4)
    raw = data[...,96:].reshape(-1,90).T
    initial = np.linalg.lstsq(basis,raw,rcond=None)[0]
    unconstrained = fit_rician_sh(basis,raw,np.full(4,.4),initial)
    assert (basis @ unconstrained).min() < -.01
    obj.fit2D4D_LLS_RealSphHarm_wSorting_norm_var(
        data,obj.mask,rank1sl=False,sigma=np.full(obj.mask.shape,.4))
    assert np.isfinite(obj.signal_4d).all()
    assert obj.signal_4d[...,6:].min() >= -1e-10

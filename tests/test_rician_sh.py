"""Numerical and SMI integration checks; no external imaging data required."""
import numpy as np
import pytest
from scipy.optimize import least_squares
from scipy.stats import rice

from lib.rician_sh import fit_rician_sh, rician_mean_and_derivative
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

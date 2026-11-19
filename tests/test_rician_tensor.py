"""Independent Rician means, known tensors, and nonlinear-fit regression checks."""
import numpy as np
import pytest
from scipy.optimize import least_squares
from scipy.special import hyp1f1

from lib.tensor import TensorFitting
from lib.rician_tensor import fit_rician_tensor, tensor_design, _factor, _fit_voxel
from tests.test_rician_sh import protocol


def independent_mean(signal, sigma):
    return signal if sigma == 0 else sigma*np.sqrt(np.pi/2)*hyp1f1(-.5, 1, -.5*(signal/sigma)**2)


@pytest.mark.parametrize('kurtosis', [False, True])
@pytest.mark.parametrize('sigma', [0., 10., 100.])
def test_recovers_exact_mean_without_tensor_constraints(kurtosis, sigma):
    b, g = protocol()
    d = np.array([-.1, .7, 1.5])  # Deliberately not positive definite.
    cumulant = -.25 if kurtosis else 0.  # Deliberately negative kurtosis.
    signal = 1000*np.exp(-b*(g*g@d) + b*b*cumulant/6)
    y = independent_mean(signal, sigma).reshape(1, 1, 1, -1)
    original = y.copy()
    mask = np.ones((1, 1, 1), bool);noise = np.full(mask.shape, sigma)
    fitter = TensorFitting(np.column_stack([g, b*1000]), 1)
    if kurtosis:
        dt, s0, x = fitter.dki_fit(y, mask, constraints=[1, 1, 1], sigma=noise)
    else:
        dt, s0, x = fitter.dti_fit(y, mask, sigma=noise)
    np.testing.assert_allclose(dt[[0, 3, 5], 0], d, atol=2e-6)
    np.testing.assert_allclose(dt[[1, 2, 4], 0], 0, atol=2e-6)
    np.testing.assert_allclose(s0, 1000, atol=2e-4)
    if kurtosis:
        c = dt[6:, 0]*dt[[0, 3, 5], 0].mean()**2
        np.testing.assert_allclose(x[:, 7:]@c, b*b*cumulant/6, atol=2e-6)
    assert fitter.rician_diagnostics['nonpositive_D'][0]
    assert not fitter.rician_diagnostics['nonconverged'][0]
    np.testing.assert_array_equal(y, original)


def test_voxel_specific_gradients_and_singleton_axes():
    b, g = protocol();grad = np.column_stack([g, b])
    mask = np.ones((2, 1, 1), bool)
    fields = np.array([np.eye(3), [[.9, .1, 0], [0, 1.1, .03], [0, 0, 1.03]]]).reshape(2, 1, 1, 3, 3)
    signals = []
    for field in fields[:, 0, 0]:
        q = (np.sqrt(b)[:, None]*g)@field.T
        attenuation = q*q@np.array([.4, .8, 1.7])
        clean = 900*np.exp(-attenuation + .6*attenuation**2/6)
        signals.append(independent_mean(clean, 35))
    dt, s0, design, info = fit_rician_tensor(np.array(signals).reshape(2,1,1,-1), mask, grad,
        np.full(mask.shape, 35.), kurtosis=True, L=fields)
    np.testing.assert_allclose(dt[[0,3,5]], np.array([[.4,.8,1.7]]).T*np.ones((1,2)), atol=1e-6)
    np.testing.assert_allclose(s0, 900, atol=1e-4)
    assert design is None and not np.any(info['nonconverged'])


def test_noisy_objective_matches_independent_numerical_jacobian_fit():
    b, g = protocol();x = tensor_design(np.column_stack([g,b]), True)
    clean = 1000*np.exp(-b*(g*g@np.array([.4,.7,1.6])) + .4*b*b/6)
    sigma = 50.
    y = independent_mean(clean, sigma) + sigma/2*np.random.default_rng(728).normal(size=len(b))
    theta, info = _fit_voxel(y, sigma, x, _factor(x))
    initial = np.linalg.pinv(x)@np.log(np.maximum(np.sqrt(np.maximum(0,y*y-sigma*sigma)),np.finfo(float).eps))
    def residual(t):
        return (independent_mean(np.exp(x@t), sigma)-y)/sigma
    reference = least_squares(residual, initial, jac='3-point', max_nfev=2000,
                              ftol=1e-11, xtol=1e-11, gtol=1e-11)
    np.testing.assert_allclose(np.sum(residual(theta)**2), 2*reference.cost, atol=1e-7)
    np.testing.assert_allclose(info['cost'], np.sum(residual(theta)**2), atol=1e-10)
    np.testing.assert_allclose(info['initial_cost'], np.sum(residual(initial)**2), atol=1e-10)
    assert info['cost'] <= info['initial_cost']


def test_parallel_is_deterministic_and_preserves_voxel_order():
    b,g=protocol();grad=np.column_stack([g,b])
    clean=np.exp(-b*(g*g@np.array([.4,.7,1.6]))+.3*b*b/6)
    amplitudes=np.linspace(700,1400,65)
    y=np.array([independent_mean(s*clean,30) for s in amplitudes]).reshape(65,1,1,-1)
    mask=np.ones(y.shape[:3],bool);sigma=np.full(mask.shape,30.)
    serial=fit_rician_tensor(y,mask,grad,sigma,kurtosis=True,n_jobs=1)
    parallel=fit_rician_tensor(y,mask,grad,sigma,kurtosis=True,n_jobs=2)
    for i in range(3):np.testing.assert_allclose(parallel[i],serial[i],rtol=1e-10,atol=1e-10)
    np.testing.assert_allclose(parallel[1],amplitudes,atol=1e-3)


@pytest.mark.parametrize('problem', ['negative_sigma','nan_sigma','nan_dwi','rank'])
def test_invalid_inputs_fail_clearly(problem):
    b,g=protocol();grad=np.column_stack([g,b]);y=np.ones((1,1,1,len(b)));sigma=np.ones((1,1,1))
    if problem=='negative_sigma':sigma[:]=-1
    if problem=='nan_sigma':sigma[:]=np.nan
    if problem=='nan_dwi':y[...,0]=np.nan
    if problem=='rank':grad[:,:3]=0
    with pytest.raises(ValueError):fit_rician_tensor(y,np.ones((1,1,1),bool),grad,sigma,kurtosis=True)

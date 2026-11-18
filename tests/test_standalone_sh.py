"""Standalone SH fits: analytic invariants, shell bookkeeping, no SM training."""
import numpy as np
import pytest
from scipy.special import eval_legendre
from scipy.stats import rice

from lib.sh import fit_rotational_invariants
from lib.smi import SMI
from tests.test_rician_sh import protocol


@pytest.mark.parametrize('rician', [False, True])
def test_explicit_order_recovers_analytic_invariants_without_training(monkeypatch, rician):
    def forbidden(*args, **kwargs):
        raise AssertionError('Standalone SH must not estimate noise or fit the SM')
    for method in ['fit', 'set_priors', 'set_sigma', 'standard_model_mlfit_rot_invs']:
        monkeypatch.setattr(SMI, method, forbidden)
    _, directions = protocol()
    directions = directions[6:96]
    amplitudes = [.6, .12, .08, .04]
    clean = sum(a * eval_legendre(l, directions[:, 2])
                for l, a in zip(range(0, 7, 2), amplitudes))
    measured = rice.mean(clean/.2)*.2 if rician else clean
    dwi = np.stack([measured, np.full_like(measured, np.nan)]).reshape(2, 1, 1, -1)
    mask = np.array([True, False]).reshape(2, 1, 1)
    images, shells = fit_rotational_invariants(
        dwi, np.ones(len(clean)), directions, lmax=6, mask=mask,
        sigma=np.full(mask.shape, .2) if rician else None, rician=rician)
    assert list(images) == ['S0', 'S2', 'S4', 'S6']
    assert shells == [dict(volume=0, b_value_s_per_mm2=1000., b_shape=1.,
                           echo_time_ms=0., n_measurements=90, lmax=6)]
    for l, a in zip(range(0, 7, 2), amplitudes):
        assert images[f'S{l}'].shape == (2, 1, 1, 1)
        np.testing.assert_allclose(images[f'S{l}'][0, 0, 0, 0], a/(2*l+1), atol=2e-7)
        assert images[f'S{l}'][1, 0, 0, 0] == 0


def test_zero_order_is_explicit_and_default_uses_all_fitted_orders():
    b, g = protocol()
    data = (1 + .2*g[:, 2]**2).reshape(1, 1, 1, -1)
    images, shells = fit_rotational_invariants(data, b, g, lmax=0)
    assert list(images) == ['S0']
    np.testing.assert_allclose(images['S0'].ravel(), [data[..., b == x].mean() for x in [0, 1, 2]])
    assert all(s['lmax'] == 0 for s in shells)
    images, shells = fit_rotational_invariants(data, b, g)
    assert list(images) == ['S0', 'S2', 'S4']
    assert [s['lmax'] for s in shells] == [0, 2, 4]
    assert images['S4'][0, 0, 0, 1] == 0


def test_multiple_echo_times_and_shapes_keep_shell_mapping():
    _, g = protocol()
    g = g[6:96]
    b = np.ones(270)
    beta = np.repeat([1., 0., 1.], 90)
    te = np.repeat([90., 60., 60.], 90)
    data = np.repeat([.4, .7, .8], 90).reshape(1, 1, 1, -1)
    images, shells = fit_rotational_invariants(data, b, np.tile(g, (3, 1)),
                                              beta=beta, echo_time=te, lmax=4)
    assert [(s['b_shape'], s['echo_time_ms'], s['lmax']) for s in shells] == [
        (0., 60., 0), (1., 60., 4), (1., 90., 4)]
    np.testing.assert_allclose(images['S0'].ravel(), [.7, .8, .4], atol=1e-7)
    assert images['S2'][0, 0, 0, 0] == images['S4'][0, 0, 0, 0] == 0


def test_unidentifiable_order_fails_before_nonlinear_fit():
    _, g = protocol()
    g = g[6:16]
    with pytest.raises(ValueError, match='needs at least 28 directions'):
        fit_rotational_invariants(np.ones((1, 1, 1, 10)), np.ones(10), g, lmax=6)
    with pytest.raises(ValueError, match='rank-deficient'):
        fit_rotational_invariants(np.ones((1, 1, 1, 90)), np.ones(90),
                                  np.tile([0., 0., 1.], (90, 1)), lmax=2)

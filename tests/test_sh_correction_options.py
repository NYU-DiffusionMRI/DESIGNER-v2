"""Optional SH positivity, scalar-shell corrections, and CLI integration."""
import json

import nibabel as nib
import numpy as np
import pytest
from scipy.special import hyp1f1

from lib.rician_sh import invert_rician_mean, nonnegative_directions
from lib.sh import fit_rotational_invariants
from lib.smi import SMI
from tests.test_tmi_rician_sh import inputs, run_cli


@pytest.mark.parametrize('scale', [1e-6, 1., 1e6])
def test_scalar_mean_inverse_matches_independent_formula(scale):
    signal = scale*np.array([0., .01, .2, 1., 2., 10., 100.])
    sigma = scale*np.array([1., .2, .2, .3, 1., 2., 3.])
    mean = sigma*np.sqrt(np.pi/2)*hyp1f1(-.5, 1, -.5*(signal/sigma)**2)
    np.testing.assert_allclose(invert_rician_mean(mean, sigma), signal, atol=1e-10*scale)
    np.testing.assert_array_equal(invert_rician_mean([-1., 0., 1.], 1), 0)
    np.testing.assert_array_equal(invert_rician_mean([-1., 0., 2.], 0), [0., 0., 2.])


@pytest.mark.parametrize('mode', ['maxzero', 'abs'])
def test_b0_and_ste_clip_individual_measurements_before_averaging(mode):
    # Unequal measurements distinguish mean(clip(M)) from clip(mean(M)).
    data = np.array([.1, 1., 2., .2, .5, 1.5]).reshape(1, 1, 1, 6)
    b = np.repeat([0., 2.], 3)
    images, shells = fit_rotational_invariants(data, b, np.zeros((6, 3)),
        beta=np.zeros(6), sigma=np.ones((1, 1, 1)), lmax=4,
        rician_clip=True, rician_clip_mode=mode)
    diff = data.astype(np.float32).astype(float)**2-1
    corrected = np.sqrt(abs(diff) if mode == 'abs' else np.maximum(0, diff))
    expected = corrected.reshape(1, 1, 1, 2, 3).mean(axis=-1)
    np.testing.assert_allclose(images['S0'], expected, atol=1e-12)
    assert all(s['lmax'] == 0 for s in shells)
    assert np.all(images['S2'] == 0) and np.all(images['S4'] == 0)


@pytest.mark.parametrize('constrained', [False, True])
def test_scalar_ste_and_b0_match_order_zero_lte_expectation(constrained):
    # The mean model uses per-measurement sigma, not sigma/sqrt(n).
    noise = .4
    truth = np.repeat([1., .3, .3], 8)
    mean = noise*np.sqrt(np.pi/2)*hyp1f1(-.5, 1, -.5*(truth/noise)**2)
    offsets = np.tile(np.linspace(-.08, .08, 8), 3)
    data = (mean+offsets).reshape(1, 1, 1, -1)
    b = np.repeat([0., 2., 2.], 8)
    beta = np.repeat([1., 0., 1.], 8)
    images, _ = fit_rotational_invariants(data, b, nonnegative_directions(24),
        beta=beta, sigma=np.full((1, 1, 1), noise), lmax=0,
        rician=True, nonnegative=constrained)
    np.testing.assert_allclose(images['S0'].ravel(), [1., .3, .3], atol=1e-6)


def test_positivity_is_opt_in_for_nontruncated_signal():
    # A positive full stick signal can have a negative finite SH projection.
    directions = nonnegative_directions(90)
    signal = np.exp(-12*directions[:, 2]**2)
    data = signal.reshape(1, 1, 1, -1)
    b = np.full(90, 6.)
    sigma = np.zeros((1, 1, 1))
    lls, _ = fit_rotational_invariants(data, b, directions, lmax=4)
    for options in [dict(rician=True), dict(rician_clip=True),
                    dict(rician_clip=True, rician_clip_mode='abs')]:
        actual, _ = fit_rotational_invariants(data, b, directions, lmax=4,
                                             sigma=sigma, **options)
        for key in lls:
            np.testing.assert_allclose(actual[key], lls[key], atol=1e-12)
    obj = SMI(b, directions)
    basis = obj.get_even_sh(directions, 4)
    coefficients = np.linalg.lstsq(basis, signal, rcond=None)[0]
    assert np.min(basis@coefficients) < -.01
    positive, _ = fit_rotational_invariants(data, b, directions, lmax=4,
        sigma=sigma, rician=True, nonnegative=True)
    assert abs(positive['S0'].item()-lls['S0'].item()) > .001


@pytest.mark.parametrize('clip', [False, True])
def test_nonnegative_requires_expectation_api(clip):
    with pytest.raises(ValueError, match='requires Rician-expectation'):
        SMI(np.ones(6), nonnegative_directions(6), rician_sh_clip=clip,
            rician_sh_nonnegative=True)


def test_cli_nonneg_changes_high_b_projection(inputs):
    directions = nonnegative_directions(90)
    signal = np.exp(-12*directions[:, 2]**2)
    data = np.broadcast_to(signal, (2, 2, 1, 90)).copy()
    nib.save(nib.Nifti1Image(data, np.eye(4)), inputs/'dwi.nii')
    nib.save(nib.Nifti1Image(np.zeros((2, 2, 1)), np.eye(4)), inputs/'sigma.nii')
    np.savetxt(inputs/'dwi.bval', [np.full(90, 6000.)])
    np.savetxt(inputs/'dwi.bvec', directions.T)
    common = ['-SH', '-lmaxSH', '4', '-RicianExp', '-sigma', 'sigma.nii']
    for name, extra in [('default', []), ('positive', ['-nonNeg'])]:
        result = run_cli(inputs, name, common+extra)
        assert result.returncode == 0, result.stdout
        meta = json.loads((inputs/name/'SH_shells.json').read_text())
        assert meta['nonnegative'] == bool(extra)
        assert ('with nonnegative signal constraints' in result.stdout) == bool(extra)
    a = nib.load(inputs/'default/S0.nii').get_fdata()
    c = nib.load(inputs/'positive/S0.nii').get_fdata()
    assert np.all(np.abs(c-a) > .001)


@pytest.mark.parametrize('estimator', ['RicianClipMax', 'RicianClipAbs', 'RicianExp'])
def test_cli_corrects_ste_shells(inputs, estimator):
    data = np.array([1., 1.1, 1.2, .15, .3, .5])
    nib.save(nib.Nifti1Image(np.broadcast_to(data, (2, 2, 1, 6)).copy(), np.eye(4)), inputs/'dwi.nii')
    np.savetxt(inputs/'dwi.bval', [[0, 0, 0, 2000, 2000, 2000]])
    np.savetxt(inputs/'dwi.bvec', np.zeros((3, 6)))
    np.savetxt(inputs/'dwi.bshape', [1, 1, 1, 0, 0, 0])
    result = run_cli(inputs, 'ste', ['-SH', '-'+estimator, '-sigma', 'sigma.nii', '-lmaxSH', '4'])
    assert result.returncode == 0, result.stdout
    observed = data.astype(np.float32).astype(float)
    if estimator == 'RicianExp':
        expected = invert_rician_mean(observed.reshape(2, 3).mean(axis=-1), .12)
    else:
        diff = observed**2-.12**2
        expected = np.sqrt(abs(diff) if estimator.endswith('Abs') else np.maximum(0, diff)).reshape(2, 3).mean(axis=-1)
    actual = nib.load(inputs/'ste/S0.nii').get_fdata()
    np.testing.assert_allclose(actual, np.broadcast_to(expected, actual.shape), atol=1e-7)
    assert np.all(nib.load(inputs/'ste/S2.nii').get_fdata() == 0)
    assert np.all(nib.load(inputs/'ste/S4.nii').get_fdata() == 0)


@pytest.mark.parametrize('multi_te', [False, True])
@pytest.mark.parametrize('estimator', ['RicianClipMax', 'RicianClipAbs'])
def test_combined_restoration_does_not_double_correct_sh(inputs, multi_te, estimator):
    from tests.test_rician_sh import protocol
    b, g = protocol()
    adc = .4*g[:, 0]**2 + .7*g[:, 1]**2 + 1.3*g[:, 2]**2
    signal = np.exp(-b*adc + .7*(b*adc)**2/6)
    data = np.broadcast_to(signal, (2, 2, 1, len(b))).copy()
    if multi_te:
        data = np.concatenate([data, .7*data], axis=-1)
        np.savetxt(inputs/'dwi.bval', [np.tile(b*1000, 2)])
        np.savetxt(inputs/'dwi.bvec', np.tile(g.T, (1, 2)))
        np.savetxt(inputs/'dwi.echotime', np.repeat([60, 90], len(b)))
    nib.save(nib.Nifti1Image(data, np.eye(4)), inputs/'dwi.nii')
    common = ['-SH', '-'+estimator, '-sigma', 'sigma.nii', '-mask', 'mask.nii', '-lmaxSH', '2']
    for name, extra in [('sh', []), ('combined', ['-DKI', '-db0restore', '-akc_outliers'])]:
        result = run_cli(inputs, name, common+extra)
        assert result.returncode == 0, result.stdout
    for name in ['S0.nii', 'S2.nii']:
        np.testing.assert_array_equal(nib.load(inputs/'sh'/name).get_fdata(),
                                      nib.load(inputs/'combined'/name).get_fdata())


@pytest.mark.parametrize('mode,multi_te', [('-SH', False), ('-SMI', True)])
def test_cli_combined_nonneg_does_not_constrain_tensor_fit(inputs, mode, multi_te):
    from tests.test_rician_sh import protocol
    if multi_te:
        image = nib.load(inputs/'dwi.nii')
        data = np.concatenate([image.get_fdata(), .8*image.get_fdata()], axis=-1)
        nib.save(nib.Nifti1Image(data, image.affine), inputs/'dwi.nii')
        b, g = protocol()
        np.savetxt(inputs/'dwi.bval', [np.tile(b*1000, 2)])
        np.savetxt(inputs/'dwi.bvec', np.tile(g.T, (1, 2)))
        np.savetxt(inputs/'dwi.echotime', np.repeat([60, 90], len(b)))
    common = [mode, '-DTI', '-DKI', '-RicianExp', '-sigma', 'sigma.nii', '-mask', 'mask.nii']
    if mode == '-SMI':
        common += ['-load_prior', 'prior.mat', '-lmax', '2']
    else:
        common += ['-lmaxSH', '0,2,4']
    for name, extra in [('unconstrained', []), ('nonnegative', ['-nonNeg'])]:
        result = run_cli(inputs, name, common+extra)
        assert result.returncode == 0, result.stdout
        assert ('with nonnegative signal constraints' in result.stdout) == bool(extra)
    tensors = [p for p in (inputs/'unconstrained').glob('*.nii')
               if '_dti' in p.name or '_dki' in p.name]
    assert len(tensors) >= 2
    for path in tensors:
        np.testing.assert_array_equal(nib.load(path).get_fdata(),
            nib.load(inputs/'nonnegative'/path.name).get_fdata())

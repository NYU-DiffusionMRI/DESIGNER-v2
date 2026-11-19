"""Clipped correction for tensor fits and explicit SH estimator precedence."""
import json

import nibabel as nib
import numpy as np
import pytest

from tests.test_rician_sh import protocol
from tests.test_tmi_rician_sh import inputs, run_cli


def compare_maps(first, second, names=None):
    names = names or [p.name for p in first.glob('*.nii')]
    assert names
    for name in names:
        a, b = nib.load(first/name), nib.load(second/name)
        np.testing.assert_allclose(a.affine, b.affine)
        assert np.isfinite(a.get_fdata()).all()
        assert np.isfinite(b.get_fdata()).all()
        np.testing.assert_allclose(a.get_fdata(), b.get_fdata(), rtol=2e-5, atol=2e-6)


@pytest.mark.parametrize('multi_te,gnc,settings,estimator', [
    (False,False,[],'Max'), (True,False,[],'Max'), (False,True,[],'Max'),
    (False,False,['-fit_constraints','0,1,0'],'Max'),
    (False,False,['-akc_outliers'],'Max'), (True,False,['-akc_outliers'],'Max'),
    (False,False,['-fit_smoothing','50'],'Max'), (True,False,['-fit_smoothing','50'],'Max'),
    (False,False,[],'Abs'), (True,False,[],'Abs'), (False,True,[],'Abs'),
    (True,False,['-akc_outliers','-fit_smoothing','50'],'Abs')])
def test_tensor_clipping_matches_explicitly_corrected_input(inputs, multi_te, gnc, settings, estimator):
    image = nib.load(inputs/'dwi.nii')
    data = image.get_fdata()
    if settings:
        b,g=protocol()
        adc=g*g@np.array([.4,.7,1.3])
        data=np.broadcast_to(np.exp(-b*adc+.7*(b*adc)**2/6),data.shape).copy()
        nib.save(nib.Nifti1Image(data,image.affine),inputs/'dwi.nii')
        nib.save(nib.Nifti1Image(np.full(data.shape[:3],.01),image.affine),inputs/'sigma.nii')
    if multi_te:
        data = np.concatenate([data, .7*data], axis=-1)
        nib.save(nib.Nifti1Image(data, image.affine), inputs/'dwi.nii')
        b, g = protocol()
        np.savetxt(inputs/'dwi.bval', [np.tile(b*1000, 2)])
        np.savetxt(inputs/'dwi.bvec', np.tile(g.T, (1, 2)))
        np.savetxt(inputs/'dwi.echotime', np.repeat([60, 90], len(b)))
    options = ['-DTI', '-DKI', '-WDKI', '-mask', 'mask.nii'] + settings
    if gnc:
        field = np.broadcast_to(np.diag([.95, 1.03, 1.01]).ravel(), (2, 2, 1, 9))
        nib.save(nib.Nifti1Image(field, image.affine), inputs/'gnc.nii')
        options += ['-gnc', 'gnc.nii', '-gnc_format', 'scanner']
    flagged = run_cli(inputs, 'flagged', options+['-RicianClip'+estimator, '-sigma', 'sigma.nii'])
    assert flagged.returncode == 0, flagged.stdout
    sigma = nib.load(inputs/'sigma.nii').get_fdata()
    difference = data**2 - sigma[..., None]**2
    corrected = np.sqrt(abs(difference) if estimator == 'Abs' else np.maximum(0, difference))
    if not settings and estimator == 'Max':
        assert np.any(corrected == 0)  # Exercises the tensor fitter's existing log floor.
    nib.save(nib.Nifti1Image(corrected, image.affine), inputs/'dwi.nii')
    explicit = run_cli(inputs, 'explicit', options)
    assert explicit.returncode == 0, explicit.stdout
    compare_maps(inputs/'flagged', inputs/'explicit')


def test_expectation_takes_precedence_for_both_tensors_and_sh(inputs):
    common = ['-DTI', '-DKI', '-SH', '-lmaxSH', '0,2,4', '-mask', 'mask.nii',
              '-sigma', 'sigma.nii']
    for name, flags in [('clip', ['-RicianClipMax']), ('absolute', ['-RicianClipAbs']),
                        ('rician', ['-RicianExp']), ('both', ['-RicianClipAbs', '-RicianExp'])]:
        result = run_cli(inputs, name, common+flags)
        assert result.returncode == 0, result.stdout
    tensor_names = [p.name for p in (inputs/'rician').glob('*.nii')
                    if p.name not in ['S0.nii', 'S2.nii', 'S4.nii']]
    compare_maps(inputs/'both', inputs/'rician', tensor_names)
    assert not np.allclose(nib.load(inputs/'rician/mk_dki.nii').get_fdata(),
                           nib.load(inputs/'clip/mk_dki.nii').get_fdata())
    compare_maps(inputs/'both', inputs/'rician', ['S0.nii', 'S2.nii', 'S4.nii'])
    standalone = run_cli(inputs, 'sh_only', ['-SH', '-RicianExp', '-lmaxSH', '0,2,4',
                                           '-mask', 'mask.nii', '-sigma', 'sigma.nii'])
    assert standalone.returncode == 0, standalone.stdout
    compare_maps(inputs/'rician', inputs/'sh_only', ['S0.nii', 'S2.nii', 'S4.nii'])
    assert json.loads((inputs/'both/SH_shells.json').read_text())['estimator'] == 'RicianExp'
    assert not np.allclose(nib.load(inputs/'clip/S0.nii').get_fdata(),
                           nib.load(inputs/'rician/S0.nii').get_fdata())
    assert not np.allclose(nib.load(inputs/'clip/S0.nii').get_fdata(),
                           nib.load(inputs/'absolute/S0.nii').get_fdata())
    assert json.loads((inputs/'absolute/SH_shells.json').read_text())['estimator'] == 'RicianClipAbs'


def test_rician_sh_also_takes_precedence_in_smi(inputs):
    common = ['-SMI', '-RicianExp', '-sigma', 'sigma.nii', '-mask', 'mask.nii',
              '-load_prior', 'prior.mat', '-lmax', '2']
    for name, flags in [('rician', []), ('both', ['-RicianClipAbs'])]:
        result = run_cli(inputs, name, common+flags)
        assert result.returncode == 0, result.stdout
        assert 'Rician expectation NLS' in result.stdout
    compare_maps(inputs/'both', inputs/'rician')


@pytest.mark.parametrize('multi_te', [False, True])
@pytest.mark.parametrize('estimator', ['Max', 'Abs'])
def test_clipping_is_not_repeated_after_b0_restoration(inputs, multi_te, estimator):
    b, g = protocol()
    diffusivity = .4*g[:, 0]**2 + .7*g[:, 1]**2 + 1.3*g[:, 2]**2
    signal = np.exp(-b*diffusivity + .7*(b*diffusivity)**2/6)
    data = np.arange(900., 1300., 100.).reshape(2, 2, 1, 1)*signal
    if multi_te:
        data = np.concatenate([data, .7*data], axis=-1)
        np.savetxt(inputs/'dwi.bval', [np.tile(b*1000, 2)])
        np.savetxt(inputs/'dwi.bvec', np.tile(g.T, (1, 2)))
        np.savetxt(inputs/'dwi.echotime', np.repeat([60, 90], len(b)))
    nib.save(nib.Nifti1Image(data, np.eye(4)), inputs/'dwi.nii')
    nib.save(nib.Nifti1Image(np.full((2, 2, 1), 20.), np.eye(4)), inputs/'sigma.nii')
    common = ['-DKI', '-db0restore', '-akc_outliers', '-mask', 'mask.nii']
    flagged = run_cli(inputs, 'flagged', common+['-RicianClip'+estimator, '-sigma', 'sigma.nii'])
    assert flagged.returncode == 0, flagged.stdout
    corrected = np.sqrt(abs(data**2-20**2) if estimator == 'Abs' else np.maximum(0, data**2-20**2))
    nib.save(nib.Nifti1Image(corrected, np.eye(4)), inputs/'dwi.nii')
    explicit = run_cli(inputs, 'explicit', common)
    assert explicit.returncode == 0, explicit.stdout
    compare_maps(inputs/'flagged', inputs/'explicit')


@pytest.mark.parametrize('options,error', [
    (['-DTI', '-RicianClipMax'], 'requires an explicit noise map'),
    (['-DTI', '-RicianClipAbs'], 'requires an explicit noise map'),
    (['-DKI', '-RicianExp'], 'requires an explicit noise map'),
    (['-DTI', '-RicianClipMax', '-sigma', 'sigma.nii', '-RicianSH_maxiter', '5'], 'with -SH/-SMI'),
    (['-SH', '-RicianClipMax', '-RicianClipAbs', '-sigma', 'sigma.nii'], 'Choose either'),
    (['-DKI', '-RicianClip', '-sigma', 'sigma.nii'], 'two explicit choices'),
])
def test_invalid_tensor_correction_options(inputs, options, error):
    result = run_cli(inputs, 'invalid', options)
    assert result.returncode != 0 and error in result.stdout

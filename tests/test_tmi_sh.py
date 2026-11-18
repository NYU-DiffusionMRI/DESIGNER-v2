"""Normal tmi entry point for standalone SH, without microstructure maps."""
import json

import nibabel as nib
import numpy as np
import pytest

from tests.test_rician_sh import protocol
from tests.test_tmi_rician_sh import inputs, run_cli


@pytest.mark.parametrize('rician', [False, True])
def test_single_shell_cli_saves_separate_invariants(inputs, rician):
    image = nib.load(inputs/'dwi.nii')
    nib.save(nib.Nifti1Image(image.get_fdata()[..., 6:96], image.affine), inputs/'dwi.nii')
    b, g = protocol()
    np.savetxt(inputs/'dwi.bval', [b[6:96]*1000])
    np.savetxt(inputs/'dwi.bvec', g[6:96].T)
    options = ['-SH', '-lmaxSH', '6', '-mask', 'mask.nii']
    if rician:
        options += ['-RicianSH', '-sigma', 'sigma.nii']
    result = run_cli(inputs, 'sh', options)
    assert result.returncode == 0, result.stdout
    assert 'SMI loading' not in result.stdout
    assert 'SMI model initialized' not in result.stdout
    assert 'nonnegative signal constraints' in result.stdout if rician else 'linear least squares' in result.stdout
    assert {p.name for p in (inputs/'sh').glob('*.nii')} == {'S0.nii', 'S2.nii', 'S4.nii', 'S6.nii'}
    for name in ['S0', 'S2', 'S4', 'S6']:
        output = nib.load(inputs/'sh'/f'{name}.nii')
        assert output.shape == (2, 2, 1, 1)
        assert np.isfinite(output.get_fdata()).all()
        np.testing.assert_allclose(output.affine, image.affine)
    metadata = json.loads((inputs/'sh/SH_shells.json').read_text())
    assert metadata['estimator'] == ('RicianSH' if rician else 'LLS')
    assert metadata['shells'] == [dict(volume=0, b_value_s_per_mm2=1000., b_shape=1.,
                                     echo_time_ms=0., n_measurements=90, lmax=6)]


def test_multi_echo_cli_preserves_every_shell_and_zero_order(inputs):
    image = nib.load(inputs/'dwi.nii')
    data = np.concatenate([image.get_fdata(), .7*image.get_fdata()], axis=-1)
    nib.save(nib.Nifti1Image(data, image.affine), inputs/'dwi.nii')
    b, g = protocol()
    np.savetxt(inputs/'dwi.bval', [np.tile(b*1000, 2)])
    np.savetxt(inputs/'dwi.bvec', np.tile(g.T, (1, 2)))
    np.savetxt(inputs/'dwi.echotime', np.repeat([60, 90], b.size))
    result = run_cli(inputs, 'sh', ['-SH', '-lmaxSH', '0'])
    assert result.returncode == 0, result.stdout
    assert {p.name for p in (inputs/'sh').glob('*.nii')} == {'S0.nii'}
    metadata = json.loads((inputs/'sh/SH_shells.json').read_text())
    output = nib.load(inputs/'sh/S0.nii').get_fdata()
    assert output.shape == (2, 2, 1, 6)
    for shell in metadata['shells']:
        ids = np.isclose(np.tile(b*1000, 2), shell['b_value_s_per_mm2']) & (
            np.repeat([60, 90], b.size) == shell['echo_time_ms'])
        np.testing.assert_allclose(output[..., shell['volume']], data[..., ids].mean(axis=-1), atol=1e-7)


@pytest.mark.parametrize('options,error', [
    (['-SH', '-SMI'], 'Select either -SH'),
    (['-SH', '-lmaxSH', '3'], 'nonnegative even integer'),
    (['-SH', '-lmaxSH', '-2'], 'nonnegative even integer'),
    (['-SH', '-lmaxSH', 'no'], 'nonnegative even integer'),
    (['-SH', '-lmaxSH', '0,2,3'], 'nonnegative even integer'),
    (['-SH', '-lmaxSH', '0,,4'], 'nonnegative even integer'),
    (['-SH', '-lmaxSH', '2,4'], 'expected 3, got 2'),
    (['-SH', '-lmax', '4'], 'use -lmaxSH for SH fitting orders'),
    (['-SH', '-lmax', '0,2,4'], 'use -lmaxSH for SH fitting orders'),
    (['-SH', '-lmax', '2', '-lmaxSH', '4'], 'use -lmaxSH for SH fitting orders'),
    (['-lmaxSH', '4'], '-lmaxSH requires -SH'),
    (['-SMI', '-lmaxSH', '0,2,4'], '-lmaxSH requires -SH'),
    (['-SMI', '-lmax', '0,2,4'], '-SMI -lmax accepts one invariant order'),
    (['-SH', '-RicianSH'], 'requires an explicit noise map'),
    (['-SH', '-RicianSH_maxiter', '2'], 'requires -RicianSH'),
    (['-SH', '-RicianSH', '-sigma', 'sigma.nii', '-RicianSH_maxiter', '0'], 'positive integer'),
])
def test_cli_rejects_invalid_standalone_options(inputs, options, error):
    result = run_cli(inputs, 'invalid', options)
    assert result.returncode != 0
    assert error in result.stdout


@pytest.mark.parametrize('rician', [False, True])
def test_cli_per_shell_orders_preserve_joint_shell_outputs(inputs, rician):
    options = ['-SH', '-lmaxSH', '0,2,4', '-mask', 'mask.nii']
    if rician:
        options += ['-RicianSH', '-sigma', 'sigma.nii']
    result = run_cli(inputs, 'mixed_orders', options)
    assert result.returncode == 0, result.stdout
    root = inputs/'mixed_orders'
    metadata = json.loads((root/'SH_shells.json').read_text())
    assert metadata['rank_one_processing']
    assert [s['lmax'] for s in metadata['shells']] == [0, 2, 4]
    assert {p.name for p in root.glob('*.nii')} == {'S0.nii', 'S2.nii', 'S4.nii'}
    s2 = nib.load(root/'S2.nii').get_fdata()
    s4 = nib.load(root/'S4.nii').get_fdata()
    assert s2.shape == s4.shape == (2, 2, 1, 3)
    assert np.all(s4[..., :2] == 0)
    assert np.any(s2[..., 1] > 0) and np.any(s4[..., 2] > 0)


@pytest.mark.parametrize('mode', ['-SH', '-SMI'])
def test_cli_iteration_cap_saves_aligned_diagnostics(inputs, mode):
    options = [mode, '-RicianSH', '-sigma', 'sigma.nii', '-mask', 'mask.nii',
               '-RicianSH_maxiter', '2', '-lmaxSH' if mode == '-SH' else '-lmax', '4']
    if mode == '-SMI':
        options += ['-load_prior', 'prior.mat']
    result = run_cli(inputs, 'capped', options)
    assert result.returncode == 0, result.stdout
    summary = json.loads((inputs/'capped/RicianSH_diagnostics.json').read_text())
    assert summary['constrained_maxiter'] == 2
    for name in ('RicianSH_early_stop', 'RicianSH_iterations', 'RicianSH_nonconverged'):
        image = nib.load(inputs/'capped'/f'{name}.nii')
        assert image.shape == (2, 2, 1, 3)
        np.testing.assert_allclose(image.affine, nib.load(inputs/'dwi.nii').affine)
        values = image.get_fdata()
        assert np.all(values[..., 0] == 0)
        assert np.all((values >= 0) & (values <= (2 if name.endswith('iterations') else 1)))
        if name in summary:
            assert np.count_nonzero(values, axis=(0,1,2)).tolist() == summary[name]

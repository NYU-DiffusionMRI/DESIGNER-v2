"""Installed tmi entry point: expected-mean tensors and option interactions."""
import json
import nibabel as nib
import numpy as np
import pytest
from scipy.special import hyp1f1
from tests.test_rician_sh import protocol
from tests.test_tmi_rician_sh import inputs,run_cli
from tests.test_tmi_rician_clip import compare_maps


@pytest.mark.parametrize('multi_te,gnc', [(False,False),(True,False),(False,True)])
def test_cli_recovers_known_tensor_from_expected_magnitudes(inputs,multi_te,gnc):
    b,g=protocol();sigma=60.;shape=(2,2,1)
    field=np.diag([.95,1.03,1.01]) if gnc else np.eye(3)
    q=g@field.T
    clean=1000*np.exp(-b*(q*q@np.array([.4,.7,1.5])))
    signals=[clean,.7*clean] if multi_te else [clean]
    observed=np.concatenate([sigma*np.sqrt(np.pi/2)*hyp1f1(-.5,1,-.5*(s/sigma)**2) for s in signals])
    nib.save(nib.Nifti1Image(np.broadcast_to(observed,shape+observed.shape).copy(),np.eye(4)),inputs/'dwi.nii')
    nib.save(nib.Nifti1Image(np.full(shape,sigma),np.eye(4)),inputs/'sigma.nii')
    flags=['-DTI','-DKI','-WDKI','-RicianExp','-sigma','sigma.nii','-mask','mask.nii',
           '-fit_constraints','1,1,1','-akc_outliers']
    if multi_te:
        np.savetxt(inputs/'dwi.bval',[np.tile(b*1000,2)])
        np.savetxt(inputs/'dwi.bvec',np.tile(g.T,(1,2)))
        np.savetxt(inputs/'dwi.echotime',np.repeat([60,90],len(b)))
    if gnc:
        nib.save(nib.Nifti1Image(np.broadcast_to(field.ravel(),shape+(9,)).copy(),np.eye(4)),inputs/'gnc.nii')
        flags+=['-gnc','gnc.nii','-gnc_format','scanner']
    result=run_cli(inputs,'expectation',flags)
    assert result.returncode==0,result.stdout
    assert 'ignoring -fit_constraints and -akc_outliers' in result.stdout
    suffixes=['_te60.0','_te90.0'] if multi_te else ['']
    for suffix in suffixes:
        for model in ['dti','dki','wdki']:
            root=inputs/'expectation';tag=model+suffix
            np.testing.assert_allclose(nib.load(root/f'md_{tag}.nii').get_fdata(),np.mean([.4,.7,1.5]),atol=2e-5)
            np.testing.assert_allclose(nib.load(root/f'ad_{tag}.nii').get_fdata(),1.5,atol=2e-5)
            if model!='dti':np.testing.assert_allclose(nib.load(root/f'mk_{tag}.nii').get_fdata(),0,atol=1e-4)
            if model!='wdki':
                info=json.loads((root/f'RicianExp_diagnostics_{tag}.json').read_text())
                assert info['constraints']=='none' and info['nonconverged']==0
    assert not list((inputs/'expectation').glob('*akc*'))


@pytest.mark.parametrize('flag,value', [('db0restore',None),('fit_smoothing','50'),('polyreg',None),('dt','unused.nii')])
def test_expectation_rejects_unimplemented_refinement_combinations(inputs,flag,value):
    options=['-DKI','-RicianExp','-sigma','sigma.nii','-'+flag]
    if value is not None:options.append(value)
    result=run_cli(inputs,'invalid',options)
    assert result.returncode!=0 and 'does not yet support -'+flag in result.stdout


def test_renamed_flag_has_no_legacy_tensor_clipping_behavior(inputs):
    result=run_cli(inputs,'invalid',['-DKI','-RicianSH','-sigma','sigma.nii'])
    assert result.returncode!=0 and 'renamed to -RicianExp' in result.stdout


def test_nonconvergence_is_saved_and_warned_without_aborting(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from unittest.mock import Mock
    from designer2 import tmi
    info = {key: np.array([value]) for key, value in dict(
        status=0, nfev=1200, cost=3., initial_cost=10., retried=True,
        nonconverged=True, nonpositive_D=True, signal_underflow=False).items()}
    save = Mock()
    monkeypatch.setattr(tmi, 'save_params', save)
    logger = Mock()
    tmi._save_tensor_rician_diagnostics(SimpleNamespace(rician_diagnostics=info),
        np.ones((1, 1, 1), bool), None, tmp_path, 'dki', logger)
    assert json.loads((tmp_path/'RicianExp_diagnostics_dki.json').read_text())['nonconverged'] == 1
    assert save.call_args.args[0]['RicianExp_nonconverged'][0, 0, 0] == 1
    logger.warning.assert_called_once()

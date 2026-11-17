"""Exercise the installed tmi entry point on small, generated magnitude images."""
import os
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest
from scipy.io import loadmat, savemat

from tests.test_rician_sh import protocol


@pytest.fixture
def inputs(tmp_path):
    b, g = protocol()
    signal = .6*np.exp(-b*2*g[:,2]**2) + .4*np.exp(-b*(.5+g[:,2]**2))
    rng = np.random.default_rng(12)
    shape = (2,2,1)
    noisy = np.hypot(signal + .12*rng.normal(size=shape+(b.size,)),
                     .12*rng.normal(size=shape+(b.size,)))
    for name, data in [('dwi', noisy), ('sigma', np.full(shape,.12)), ('mask', np.ones(shape))]:
        nib.save(nib.Nifti1Image(data, np.eye(4)), tmp_path/f'{name}.nii')
    np.savetxt(tmp_path/'dwi.bval', [b*1000])
    np.savetxt(tmp_path/'dwi.bvec', g.T)
    priors = loadmat(Path(__file__).parents[1]/'constant/SMI_Gaussian_noFWPrior.mat')
    prior = next(v for k,v in priors.items() if not k.startswith('__'))
    savemat(tmp_path/'prior.mat', {'prior':prior[:3000]})
    return tmp_path


def run_cli(root, name, options):
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', SMI_BLAS_THREADS='1')
    cmd = ['tmi', 'dwi.nii', name, '-n_cores', '1'] + options
    result = subprocess.run(cmd, cwd=root, env=env, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180)
    (root/f'{name}.log').write_text(result.stdout)
    return result


@pytest.mark.parametrize('multi_te', [False, True])
def test_cli_both_estimators(inputs, multi_te):
    if multi_te:
        dwi = nib.load(inputs/'dwi.nii')
        nib.save(nib.Nifti1Image(np.concatenate((dwi.get_fdata(), .8*dwi.get_fdata()), axis=-1),
                               dwi.affine), inputs/'dwi.nii')
        b,g=protocol()
        np.savetxt(inputs/'dwi.bval', [np.tile(b*1000,2)])
        np.savetxt(inputs/'dwi.bvec', np.tile(g.T,(1,2)))
        np.savetxt(inputs/'dwi.echotime', np.repeat([60,90], b.size))
    common = ['-SMI','-sigma','sigma.nii','-mask','mask.nii','-load_prior','prior.mat','-lmax','2']
    for name, extra in [('lls',[]),('rician',['-RicianSH'])]:
        result=run_cli(inputs,name,common+extra)
        assert result.returncode == 0, result.stdout
        assert ('Rician expectation NLS' if extra else 'linear least squares') in result.stdout
        if extra:
            assert 'nonnegative signal constraints' in result.stdout
    lls_files=sorted((inputs/'lls').glob('*.nii'))
    assert {p.name for p in lls_files} == {p.name for p in (inputs/'rician').glob('*.nii')}
    for p in lls_files:
        x=nib.load(p); y=nib.load(inputs/'rician'/p.name)
        assert x.shape == y.shape
        assert np.all(np.isfinite(y.get_fdata()))
        np.testing.assert_allclose(x.affine,y.affine)
    a=nib.load(inputs/'lls/rotinvs_smi.nii').get_fdata()
    b=nib.load(inputs/'rician/rotinvs_smi.nii').get_fdata()
    assert not np.allclose(a,b)
    np.testing.assert_array_equal(a[...,0],b[...,0])


@pytest.mark.parametrize('options, error', [
    (['-RicianSH'], 'requires -SMI'),
    (['-SMI','-RicianSH'], 'requires an explicit noise map'),
    (['-SMI','-RicianSH','-sigma','sigma.nii'], 'nonnegative')])
def test_cli_errors(inputs, options, error):
    if '-sigma' in options:
        nib.save(nib.Nifti1Image(np.full((2,2,1),-1.),np.eye(4)),inputs/'sigma.nii')
    result=run_cli(inputs,'invalid',options)
    assert result.returncode != 0
    assert error in result.stdout


def test_cli_rejects_misaligned_map(inputs):
    affine=np.eye(4); affine[0,3]=10
    nib.save(nib.Nifti1Image(np.ones((2,2,1)),affine),inputs/'sigma.nii')
    result=run_cli(inputs,'misaligned',['-SMI','-RicianSH','-sigma','sigma.nii'])
    assert result.returncode != 0
    assert 'aligned to the DWI grid' in result.stdout

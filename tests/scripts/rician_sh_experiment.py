"""Reproducible Standard Model / Rician noise experiment through the tmi CLI.

Run in the dev container: python tests/scripts/rician_sh_experiment.py OUTPUT
Requires the normal tmi installation, nibabel, scipy and matplotlib.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from scipy.special import eval_legendre


def sphere(n):
    z = 1 - 2*(np.arange(n)+.5)/n
    phi = np.arange(n)*np.pi*(3-np.sqrt(5))
    return np.column_stack((np.sqrt(1-z*z)*np.cos(phi), np.sqrt(1-z*z)*np.sin(phi), z))


def save(root, name, data):
    nib.save(nib.Nifti1Image(np.asarray(data, np.float32), np.eye(4)), root/name)


def run_tmi(root, image, name, extra=()):
    command = ['tmi', str(root/image), str(root/name), '-SMI', '-sigma', str(root/'sigma.nii'),
               '-mask', str(root/'mask.nii'), '-n_cores', '1', *extra]
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', SMI_BLAS_THREADS='1')
    start = time.monotonic()
    with (root/f'{name}.log').open('w') as log:
        result = subprocess.run(command, cwd=root, env=env, stdout=log, stderr=subprocess.STDOUT)
    record = {'command':command, 'seconds':time.monotonic()-start, 'returncode':result.returncode}
    (root/f'{name}_command.json').write_text(json.dumps(record, indent=2))
    if result.returncode:
        raise RuntimeError(f'tmi failed: see {root/name}.log')
    return record


def experiment(root, repeats=128, reuse=False):
    root.mkdir(parents=True, exist_ok=True)
    snrs = np.array([8., 15., 30., 60.])
    params = [(.45,1.8,1.5,.6,2.), (.6,2.2,1.7,.45,4.), (.75,2.6,2.,.3,8.)]
    b = np.repeat([0.,1.,2.,3.], [6,60,60,60])
    dirs = np.vstack((np.zeros((6,3)), np.tile(sphere(60),(3,1))))
    orientations = sphere(4096)
    dot2 = (dirs @ orientations.T)**2
    clean = np.zeros((len(params), b.size))
    truth = {key:[] for key in ['f','Da','DePar','DePerp','fw','p2','p4']}
    for i,(f,da,depar,deperp,kappa) in enumerate(params):
        weights = np.exp(kappa*(orientations[:,2]**2-1)); weights /= weights.sum()
        # Direct orientation quadrature of the Standard Model, independent of SH fitting.
        intra = np.exp(-b[:,None]*da*dot2)
        extra = np.exp(-b[:,None]*(deperp+(depar-deperp)*dot2))
        clean[i] = (f*intra+(1-f)*extra) @ weights
        values = [f,da,depar,deperp,0., weights@eval_legendre(2,orientations[:,2]),
                  weights@eval_legendre(4,orientations[:,2])]
        for key,val in zip(truth,values):
            truth[key].append(float(val))
    shape = (snrs.size,len(params),repeats)
    s0 = 1000.
    clean = np.broadcast_to(s0*clean[None,:,None,:], shape+(b.size,)).copy()
    sigma = np.broadcast_to(s0/snrs[:,None,None],shape).copy()
    rng = np.random.default_rng(20261117)
    noisy = np.hypot(clean+sigma[...,None]*rng.standard_normal(clean.shape),
                     sigma[...,None]*rng.standard_normal(clean.shape))
    corrected = np.sqrt(np.abs(noisy**2 - sigma[...,None]**2))
    metadata = {'seed':20261117, 'snr_definition':'noise-free b0 / complex-component sigma',
                'snrs':snrs.tolist(), 'repeats_per_tissue_per_snr':repeats,
                'tissues':params, 'truth':truth, 'directions_per_diffusion_shell':60,
                'b_s_per_mm2':[0,1000,2000,3000], 'odf':'Watson about z; 4096-point sphere quadrature',
                'S0':s0, 'training':'unmodified default SMI prior and seed; IAS,EAS',
                'noise':'sqrt((S+sigma*N1)^2+(sigma*N2)^2); independent standard Gaussian N1,N2',
                'corrected_lls':'sqrt(abs(magnitude**2 - sigma**2)) before the default LLS tmi fit'}
    if reuse and (root/'experiment.json').exists():
        previous=json.loads((root/'experiment.json').read_text())
        for key in ['seed','snrs','repeats_per_tissue_per_snr','tissues','S0']:
            if previous[key] != json.loads(json.dumps(metadata[key])):
                raise ValueError(f'Cannot reuse fits with changed {key}')
    for filename,data in [('clean.nii',clean),('noisy.nii',noisy),('corrected.nii',corrected),('sigma.nii',sigma),
                          ('mask.nii',np.ones(shape))]:
        save(root,filename,data)
    for stem in ['clean','noisy','corrected']:
        np.savetxt(root/f'{stem}.bval',[b*1000])
        np.savetxt(root/f'{stem}.bvec',dirs.T)
    for key,values in truth.items():
        save(root,f'truth_{key}.nii',np.broadcast_to(np.array(values)[None,:,None],shape))
    (root/'experiment.json').write_text(json.dumps(metadata,indent=2))
    for image,name,extra in [('clean.nii','clean_lls',[]),('noisy.nii','lls',[]),
                              ('noisy.nii','rician',['-RicianSH']),('corrected.nii','corrected_lls',[])]:
        record=root/f'{name}_command.json'
        if reuse and record.exists() and json.loads(record.read_text())['returncode']==0:
            print(f'Reusing completed {name}',flush=True)
        else:
            print(f'Running {name}',flush=True)
            run_tmi(root,image,name,extra)
    rows=[]
    for filename in sorted((root/'lls').glob('*.nii')):
        key=filename.name.removesuffix('_smi.nii')
        baseline=nib.load(root/'clean_lls'/filename.name).get_fdata()
        for method in ['lls','corrected_lls','rician']:
            fitted=nib.load(root/method/filename.name).get_fdata()
            for si,snr in enumerate(snrs):
                if key in truth:
                    targets=[('SM truth', np.broadcast_to(np.array(truth[key])[:,None],fitted[si].shape)),
                             ('noiseless tmi',baseline[si])]
                else:
                    targets=[('noiseless tmi',baseline[si])]
                for target,values in targets:
                    channels=range(fitted.shape[-1]) if fitted.ndim==4 else [None]
                    for channel in channels:
                        x=fitted[si] if channel is None else fitted[si,...,channel]
                        y=values if channel is None else values[...,channel]
                        error=x-y
                        rows.append({'map':filename.name,'channel':channel if channel is not None else '',
                                     'method':method,'snr':snr,'target':target,'n':x.size,
                                     'mean':float(x.mean()),'bias':float(error.mean()),
                                     'rmse':float(np.sqrt(np.mean(error**2))),
                                     'mae':float(np.mean(np.abs(error)))})
    with (root/'metrics.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    methods=[('lls','#d67b22','Uncorrected LLS'),
             ('corrected_lls','#22875e','Corrected LLS'),('rician','#2679b5','RicianSH')]
    for metric in ['bias','rmse']:
        fig,axes=plt.subplots(2,4,figsize=(15,7),constrained_layout=True)
        for ax,key in zip(axes.flat,truth):
            for method,color,label in methods:
                selected=[r for r in rows if r['map']==key+'_smi.nii' and r['method']==method and r['target']=='SM truth']
                ax.plot([r['snr'] for r in selected],[r[metric] for r in selected],'o-',color=color,label=label)
            if metric=='bias':ax.axhline(0,color='.5',lw=.7)
            ax.set_title(key);ax.set_xlabel('b0 SNR');ax.set_ylabel(metric.upper())
            ax.set_xscale('log');ax.minorticks_off();ax.set_xticks(snrs,labels=[str(int(v)) for v in snrs])
        axes.flat[0].legend();axes.flat[-1].axis('off')
        fig.suptitle(f'Standard Model truth: 3 tissues × {repeats} independent Rician realizations per SNR')
        fig.savefig(root/f'parameter_{metric}.png',dpi=160);plt.close(fig)
    selected=[r for r in rows if r['map']=='rotinvs_smi.nii' and r['channel'] in (1,2,3)]
    fig,axes=plt.subplots(1,3,figsize=(12,3.8),constrained_layout=True)
    for ax,ch in zip(axes,[1,2,3]):
        for method,color,label in methods:
            v=[r for r in selected if r['channel']==ch and r['method']==method]
            ax.plot([r['snr'] for r in v],[r['bias']/s0 for r in v],'o-',label=label,color=color)
        ax.axhline(0,color='.5',lw=.7);ax.set_title(f'Spherical mean, b={ch*1000}');ax.set_xlabel('b0 SNR');ax.set_ylabel('Bias / S0')
        ax.set_xscale('log');ax.minorticks_off();ax.set_xticks(snrs,labels=[str(int(v)) for v in snrs])
    axes[0].legend();fig.savefig(root/'spherical_mean_bias.png',dpi=160);plt.close(fig)
    print(f'Saved maps, metrics and figures to {root}',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('--repeats',type=int,default=128)
    parser.add_argument('--reuse-existing',action='store_true',help='Reuse completed fits for the same generated experiment; use only when estimator code is unchanged')
    args=parser.parse_args()
    experiment(args.output.resolve(),args.repeats,args.reuse_existing)

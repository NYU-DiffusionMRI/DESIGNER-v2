"""Compare every NIfTI map/channel emitted by two tmi runs.

Usage: python tests/scripts/compare_rician_sh_maps.py LLS RICIAN MASK OUTPUT
"""
import argparse
import csv
import html
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np


def mask_on_grid(mask_image, reference):
    """Match voxel-axis order without interpolation; reject different grids."""
    orientation=nib.orientations.ornt_transform(
        nib.orientations.io_orientation(mask_image.affine),
        nib.orientations.io_orientation(reference.affine))
    aligned=mask_image.as_reoriented(orientation)
    if aligned.shape!=reference.shape[:3] or not np.allclose(aligned.affine,reference.affine,atol=1e-4):
        raise ValueError('Mask grid differs from output maps after voxel-axis reorientation')
    return aligned


def compare(lls, rician, mask_file, output, label="LLS", description=""):
    output.mkdir(parents=True,exist_ok=True)
    left={p.name:p for p in lls.glob('*.nii*')}
    right={p.name:p for p in rician.glob('*.nii*')}
    if left.keys()!=right.keys():
        raise ValueError(f'Output maps differ: {left.keys() ^ right.keys()}')
    if not left:
        raise ValueError('No NIfTI output maps found')
    reference=nib.load(left[sorted(left)[0]])
    mask_image=mask_on_grid(nib.load(mask_file),reference)
    mask=mask_image.get_fdata()>0
    # MRtrix may write a different voxel-axis orientation than the input mask.
    # Reorient by exact flips/permutations above, never by resampling.
    planes=[(2,int(np.argmax(mask.sum(axis=(0,1)))),'axial'),
            (1,int(np.argmax(mask.sum(axis=(0,2)))),'coronal')]
    rows=[]; figures=[]; scalar=[]; differences=output/'differences';differences.mkdir(exist_ok=True)
    for name in sorted(left):
        aimg=nib.load(left[name]);bimg=nib.load(right[name])
        if aimg.shape!=bimg.shape or not np.allclose(aimg.affine,bimg.affine):
            raise ValueError(f'Geometry differs: {name}')
        if aimg.shape[:3]!=mask.shape or not np.allclose(aimg.affine,mask_image.affine,atol=1e-4):
            raise ValueError(f'Mask grid differs: {name}')
        a=aimg.get_fdata();b=bimg.get_fdata()
        delta=b-a
        nib.save(nib.Nifti1Image(delta,aimg.affine),differences/name)
        flat_a=a.reshape(mask.shape+(-1,));flat_b=b.reshape(mask.shape+(-1,))
        n=flat_a.shape[-1]
        for c in range(n):
            x=flat_a[...,c][mask];y=flat_b[...,c][mask]
            valid=np.isfinite(x)&np.isfinite(y);d=y[valid]-x[valid]
            rows.append({'map':name,'channel':c,'voxels':int(mask.sum()),'finite_pairs':int(valid.sum()),
                         'lls_mean':float(x[valid].mean()),'rician_mean':float(y[valid].mean()),
                         'lls_median':float(np.median(x[valid])),'rician_median':float(np.median(y[valid])),
                         'mean_delta':float(d.mean()),'median_delta':float(np.median(d)),
                         'mae':float(np.mean(np.abs(d))),'rms_delta':float(np.sqrt(np.mean(d*d))),
                         'p95_abs_delta':float(np.percentile(np.abs(d),95)),
                         'max_abs_delta':float(np.max(np.abs(d))),'exactly_equal':bool(np.array_equal(x,y))})
        if n==1:scalar.append((name,a,b))
        for axis,index,plane in planes:
            fig,axes=plt.subplots(n,3,figsize=(11,3*n),squeeze=False,constrained_layout=True)
            for c in range(n):
                x=flat_a[...,c];y=flat_b[...,c];d=y-x
                limits=np.percentile(np.r_[x[mask],y[mask]],[1,99])
                if limits[1]<=limits[0]:limits[1]=limits[0]+1
                extent=max(float(np.percentile(np.abs(d[mask]),99)),1e-12)
                selected_mask=np.take(mask,index,axis=axis).T
                for j,(image,cm,lo,hi,title) in enumerate([
                    (x,'viridis',*limits,label),(y,'viridis',*limits,'RicianSH'),
                    (d,'coolwarm',-extent,extent,f'RicianSH − {label}')]):
                    slice_=np.ma.masked_where(~selected_mask,np.take(image,index,axis=axis).T)
                    art=axes[c,j].imshow(slice_,origin='lower',cmap=cm,vmin=lo,vmax=hi)
                    axes[c,j].set_title(f'{title} · channel {c}');axes[c,j].axis('off')
                    fig.colorbar(art,ax=axes[c,j],shrink=.75)
            fig.suptitle(f'{name} · native {plane} slice {index} · display clipped to percentiles')
            filename=f'{name.split(".nii")[0]}_{plane}.png'
            fig.savefig(output/filename,dpi=120);plt.close(fig);figures.append(filename)
    with (output/'map_metrics.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (output/'map_metrics.json').write_text(json.dumps(rows,indent=2))
    # Compact overview of every scalar SMI map.
    fig,axes=plt.subplots(len(scalar),3,figsize=(10,2.7*len(scalar)),squeeze=False,constrained_layout=True)
    axis,index,_=planes[1];ms=np.take(mask,index,axis=axis).T
    for row,(name,a,b) in enumerate(scalar):
        limits=np.percentile(np.r_[a[mask],b[mask]],[1,99]);limits[1]=max(limits[1],limits[0]+1e-8)
        d=b-a;extent=max(float(np.percentile(np.abs(d[mask]),99)),1e-12)
        for col,(v,cm,lo,hi) in enumerate([(a,'viridis',*limits),(b,'viridis',*limits),(d,'coolwarm',-extent,extent)]):
            art=axes[row,col].imshow(np.ma.masked_where(~ms,np.take(v,index,axis=axis).T),origin='lower',cmap=cm,vmin=lo,vmax=hi)
            axes[row,col].axis('off');axes[row,col].set_title(name if col==0 else ['','RicianSH',f'RicianSH − {label}'][col])
            fig.colorbar(art,ax=axes[row,col],shrink=.8)
    fig.suptitle(f'Real-data comparison · {label} and RicianSH · identical pair scales')
    fig.savefig(output/'scalar_overview.png',dpi=140);plt.close(fig)
    intro=f"""<h1>Real-data RicianSH comparison</h1>
<p>{html.escape(description)}</p>
<p>Reference: {html.escape(label)}. Every output map and every 4D channel is included.
Delta means RicianSH minus the reference. Statistics use all mask voxels; images use native
voxel planes (mask reoriented by exact axis flips/permutations if needed), common 1st–99th percentile limits for each pair, and a separate symmetric
difference scale. Channel indices start at zero. Rotational invariants are ordered by
SH band, then shell. No ground truth is available for these real data.</p>
<p><a href="map_metrics.csv">All numerical metrics (CSV)</a> ·
<a href="../provenance.json">Input and noise provenance</a></p>"""
    table='<table><tr>'+''.join(f'<th>{k}</th>' for k in ['Map','Channel',label+' mean','RicianSH mean','Mean delta','RMS delta','Max |delta|'])+'</tr>'
    for r in rows:
        values=[r['map'],r['channel']]+[f'{r[k]:.6g}' for k in ['lls_mean','rician_mean','mean_delta','rms_delta','max_abs_delta']]
        table+='<tr>'+''.join(f'<td>{html.escape(str(v))}</td>' for v in values)+'</tr>'
    table+='</table>'
    body=intro+table+'<h2>All scalar maps</h2><img src="scalar_overview.png">'
    for figure in figures:body+=f'<details><summary>{html.escape(figure)}</summary><img loading="lazy" src="{figure}"></details>'
    (output/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>RicianSH map comparison</title><style>body{font:16px system-ui;max-width:1200px;margin:40px auto;padding:0 24px;color:#223}p{line-height:1.5}table{border-collapse:collapse;font-size:13px}td,th{padding:7px;border-bottom:1px solid #ddd;text-align:right}td:first-child{text-align:left}img{width:100%;height:auto}summary{cursor:pointer;padding:14px;background:#eef2f5;margin:10px 0}</style>'+body)
    print(f'Compared {len(left)} maps / {len(rows)} channels; saved to {output}',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ['lls','rician','mask','output']:parser.add_argument(key,type=Path)
    parser.add_argument('--label',default='LLS')
    parser.add_argument('--description',default='')
    args=parser.parse_args();compare(args.lls,args.rician,args.mask,args.output,args.label,args.description)

import lib.io as mio
import nibabel as nib
import numpy as np
import os

def parallel_outlier_smooth(inds, kernel, outlier_locations, dwi_norm, dwi, smoothlevel):
        
    k = kernel // 2
    x = inds[0]
    y = inds[1]
    z = inds[2]

    akcpatch = True
    while np.all(akcpatch) == True:
        xmin = 0 if x-k-1 < 0 else x-k-1
        xmax = dwi.shape[0] if x+k > dwi.shape[0] else x+k
        ymin = 0 if y-k-1 < 0 else y-k-1
        ymax = dwi.shape[1] if y+k > dwi.shape[1] else y+k
        zmin = 0 if z-k-1 < 0 else z-k-1
        zmax = dwi.shape[2] if z+k > dwi.shape[2] else z+k
    
        psize = (xmax - xmin) * (ymax - ymin) * (zmax - zmin)
        akcpatch = outlier_locations[xmin:xmax, ymin:ymax, zmin:zmax].flatten()
        k += 2


    ref = np.tile(np.reshape(dwi_norm[x,y,z,:],(1,dwi.shape[-1])),(psize,1))
    patch = np.reshape(dwi_norm[xmin:xmax,ymin:ymax,zmin:zmax,:],(psize, dwi.shape[-1]))
    patchorig = np.reshape(dwi[xmin:xmax,ymin:ymax,zmin:zmax,:],(psize, dwi.shape[-1]))
    intensities = np.sqrt(np.sum((patch-ref)**2, axis=1)) / dwi.shape[-1]
        
    min_idx = np.argsort(intensities)
    min_wgs = intensities[min_idx]
    wgs_max = min_wgs[-1]
    min_wgs[akcpatch] = wgs_max
    
    if not smoothlevel:
        goodidx = min_wgs <= np.median(min_wgs)
    else:
        goodidx = min_wgs <= np.percentile(min_wgs, smoothlevel)

    
    min_idx = min_idx[goodidx]
    min_wgs = min_wgs[goodidx]
    wgs_max = np.max(min_wgs)
    wgs_inv = wgs_max - min_wgs

    wgs_nrm = wgs_inv/np.sum(wgs_inv)
    wval = (patchorig[min_idx,:] * 
            (wgs_nrm[...,None] @ np.ones((1,dwi.shape[-1])))
            ).sum(axis=0)

    return wval

def refit_or_smooth(outlier_locations, dwi, mask=None, smoothlevel=None, n_cores=-3):
    from joblib import Parallel, delayed

    if mask is None:
        outinds = np.array(np.where(outlier_locations == 1))
    else:
        outinds = np.array(np.where(mask == 1))


    dwi_norm = abs(dwi) / np.amax(dwi, axis=(0,1,2))
    dwi_new = dwi.copy()
    kernel = 7
    # for i in range(len(outinds[0])):
    #     wval = parallel_outlier_smooth(outinds[:,i], kernel, outlier_locations, dwi_norm, dwi, smoothlevel)
    
    wval = (Parallel(n_jobs=n_cores, prefer='processes')
            (delayed(parallel_outlier_smooth)(
                outinds[:,i], kernel, outlier_locations, dwi_norm, dwi, smoothlevel
            ) for i in range(len(outinds[0]))))

    dwi_new[outinds[0,:],outinds[1,:],outinds[2,:],:] = np.array(wval)

    return dwi_new

#=======================================tmi black voxel======================================
def parallel_outlier_akc(inds, akc_mask_tmp, akc_dirs,rk,md,fa):
    import numpy as np
    k = 7 // 2
    x = inds[0]
    y = inds[1]
    z = inds[2]
    O=akc_mask_tmp[x,y,z]
    val_rk=rk[x,y,z]
    val_md=md[x,y,z]
    val_fa=fa[x,y,z]
    max_akc=np.max(akc_dirs[x,y,z,:])
    min_akc=np.min(akc_dirs[x,y,z,:])
    # if max_akc>10 or min_akc<0:
    #     O=1

    if val_rk<0.5 and val_md<2:
        O=1
    return int(O)

def akc_out(outlier_inds, akc_mask_tmp, akc_dirs,rk,md,fa, n_cores=-3):
    from joblib import Parallel, delayed
    import numpy as np

    new_akc_mask=akc_mask_tmp.copy()

    wval = (Parallel(n_jobs=n_cores, prefer='processes')
            (delayed(parallel_outlier_akc)(
                outlier_inds[:,i], akc_mask_tmp, akc_dirs,rk,md,fa
            ) for i in range(len(outlier_inds[0]))))
    
    wval=np.asarray(wval)
    wval=wval.astype(int)
    new_akc_mask[outlier_inds[0,:],outlier_inds[1,:],outlier_inds[2,:]] = wval[:]
    
    return new_akc_mask


def parallel_outlier_slope(inds, kernel, outlier_locations, bval, dwi_norm, dwi, fa,md, md_mask, smoothlevel):
    import numpy as np
    import warnings

    warnings.filterwarnings("ignore",
    message="Mean of empty slice",
    category=RuntimeWarning)

    x, y, z = inds
    k = kernel // 2
    # --- precomputed shells ---
    bval_rounded = np.round(bval, 2)
    nonb0shell = bval > 0.05
    lowbval = np.sort(np.unique(bval_rounded[bval_rounded < 1.2]))
    bx = lowbval[-1]
    bxx = lowbval[0]
    bxshell = bval_rounded == bx
    bxxshell = bval_rounded == bxx

    # --- grow patch ---
    while True:
        xmin = max(x - k - 1, 0)
        xmax = min(x + k, dwi.shape[0])
        ymin = max(y - k - 1, 0)
        ymax = min(y + k, dwi.shape[1])
        zmin = max(z - k - 1, 0)
        zmax = min(z + k, dwi.shape[2])

        akcpatch = outlier_locations[xmin:xmax, ymin:ymax, zmin:zmax]

        if not akcpatch.all():
            break
        k += 2

    # --- flatten patches ---
    psize = akcpatch.size

    fapatch = fa[xmin:xmax, ymin:ymax, zmin:zmax].ravel()
    mdpatch = md[xmin:xmax, ymin:ymax, zmin:zmax].ravel()
    csfpatch = md_mask[xmin:xmax, ymin:ymax, zmin:zmax].ravel()

    # --- similarity intensities ---
    patch = dwi_norm[xmin:xmax, ymin:ymax, zmin:zmax, nonb0shell]
    patch = patch.reshape(psize, -1)

    ref = dwi_norm[x, y, z, nonb0shell][None, :]
    diff = patch - ref
    intensities = np.sqrt((diff * diff).sum(axis=1)) / patch.shape[1]

    # --- FA / MD rejection mask ---
    omit = (
        (fapatch > fa[x,y,z] + fapatch.std()) |
        (fapatch < fa[x,y,z] - fapatch.std()) |
        (mdpatch > md[x,y,z] + mdpatch.std()) |
        (mdpatch < md[x,y,z] - mdpatch.std())
    )

    # --- rank + exclude ---
    min_wgs = intensities.copy()
    wgs_max = min_wgs.max()
    min_wgs[akcpatch.ravel()] = wgs_max
    min_wgs[csfpatch.ravel()] = wgs_max
    min_wgs[omit] = wgs_max

    if not smoothlevel:
        goodidx = min_wgs <= min_wgs.mean()
    else:
        thr = np.percentile(min_wgs, smoothlevel)
        goodidx = (min_wgs <= thr) & (min_wgs != wgs_max)

    # --- slopes ---
    bx_bv = np.log(dwi[x,y,z,bxshell].mean())
    bxx_bv = np.log(dwi[x,y,z,bxxshell])

    slope_bv = np.abs((bx_bv - bxx_bv) / (bx - bxx))

    patch_bx = np.log(dwi[xmin:xmax,ymin:ymax,zmin:zmax][:,:,:,bxshell]).reshape(psize, -1)
    patch_bxx = np.log(dwi[xmin:xmax,ymin:ymax,zmin:zmax][:,:,:,bxxshell]).reshape(psize, -1)

    bx_ok = patch_bx[goodidx].mean(axis=1)
    bxx_ok = patch_bxx[goodidx]

    slope_ok = np.abs((bx_ok[:,None] - bxx_ok) / (bx - bxx))

    mask = slope_ok > slope_bv
    masked = np.where(mask, slope_ok, np.nan)
    mean_vals = np.nanmean(masked, axis=0)
    slope_ok2 = np.where(np.isnan(mean_vals), slope_bv, mean_vals)

    wval_log = np.abs(slope_ok2) * bx + bx_bv
    return np.exp(wval_log).squeeze()

def b0restore_slope(outlier_locations, dwi, bval,k,perc,fa,md, mask=None, n_cores=-3):
    from joblib import Parallel, delayed
    import numpy as np
    from tqdm import tqdm

    if mask is None:
        outinds = np.array(np.where(outlier_locations == 1))
    else:
        outinds = np.array(np.where(mask == 1))

    dwi_norm = abs(dwi) / np.amax(dwi, axis=(0,1,2))
    dwi_new = dwi.copy()
    kernel = k
    
    csf_t=1.2
    md_mask = md.copy()
    md_mask[md>=csf_t] = int(1)
    md_mask[md<csf_t] = int(0)
    md_mask = md_mask.astype(bool)
    
    # print('=====================parallel processing=====================')
    # print(len(outinds[0]))
    wval = (Parallel(n_jobs=n_cores, prefer='processes')
            (delayed(parallel_outlier_slope)(
                outinds[:,i], kernel, outlier_locations, bval,dwi_norm, dwi,fa, md,md_mask, perc
            ) for i in tqdm(range(len(outinds[0])))))

    wval=np.asarray(wval)
    b0idx=np.where(bval<0.01)[0]
    # print(wval.shape) #[#outliers, #b0s]
    b0s=dwi[:,:,:,bval<0.01]
    for i in range(len(b0idx)):
        dwi_new[outinds[0,:],outinds[1,:],outinds[2,:],b0idx[i]] = wval[:,i]

    return dwi_new
#=======================================tmi black voxel======================================


def save_params(paramDict, niiex, model, outdir, format='nifti'):
    
    params = paramDict.keys()
    for key in params:
        if 'L' in key:
            outpath = os.path.join(r"{}".format(outdir), ("%s.mif" % (key)))
        else:
            outpath = os.path.join(r"{}".format(outdir), ("%s_%s.mif" % (key, model)))

        if not os.path.exists(outpath):
            vol = paramDict[key]
            ndims = vol.ndim

            mif = mio.Image(
                data = paramDict[key],
                vox = niiex.vox,
                transform = niiex.transform,
                grad = niiex.grad
            )

            if format == 'mif':
                mif.save(outpath)

            elif format == 'nifti':

                v = np.array(mif.vox)
                if ndims == 4:
                    v[-1] = 1
                elif ndims == 3:
                    v = np.append(v, 1)
                
                nii = nib.Nifti1Image(mif.data, mif.transform @ np.diag(v))
                nib.save(nii, outpath.replace('.mif', '.nii'))



def parallel_outlier_smooth(inds, kernel, outlier_locations, dwi_norm, dwi, smoothlevel):
    import numpy as np
    
    
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
    import numpy as np

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

def save_params(paramDict, niiex, model, outdir):
    from ants import from_numpy, image_write
    import os

    params = paramDict.keys()
    for key in params:
        outpath = os.path.join(outdir, ('%s_%s.nii' % (key, model)))
        vol = paramDict[key]
        ndims = vol.ndim

        out = from_numpy(
        paramDict[key], origin=niiex.origin[:ndims], spacing=niiex.spacing[:ndims], direction=niiex.direction[:ndims,:])
        image_write(out, outpath)


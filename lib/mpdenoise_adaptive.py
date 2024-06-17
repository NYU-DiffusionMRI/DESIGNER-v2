"""A Python class containing an implimentation of MPPCA denoising.
        
    Inputs are a 4D image with dimentions (X x Y x Z x N)
    
    Usage:
    import mpdenoise as mp
    imgdn, sigma, nparameters = mp.denoise(img, kernel=[7,7,7], M=60, shrinkage='frobenius')
    
    LICENCE
    Authors: Benjamin Ades-Aron (Benjamin.Ades-Aron@nyulangone.org)
    Copyright (c) 2016 New York University
    
    Permission is hereby granted, free of charge, to any non-commercial entity
    ('Recipient') obtaining a copy of this software and associated
    documentation files (the 'Software'), to the Software solely for
    non-commercial research, including the rights to use, copy and modify the
    Software, subject to the following conditions:
    
    1. The above copyright notice and this permission notice shall be
    included by Recipient in all copies or substantial portions of the
    Software.
    
    2. THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIESOF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
    NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BELIABLE FOR ANY CLAIM,
    DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
    OTHERWISE, ARISING FROM, OUT OF ORIN CONNECTION WITH THE SOFTWARE OR THE
    USE OR OTHER DEALINGS IN THE SOFTWARE.
    
    3. In no event shall NYU be liable for direct, indirect, special,
    incidental or consequential damages in connection with the Software.
    Recipient will defend, indemnify and hold NYU harmless from any claims or
    liability resulting from the use of the Software by recipient.

    4. Neither anything contained herein nor the delivery of the Software to
    recipient shall be deemed to grant the Recipient any right or licenses
    under any patents or patent application owned by NYU.

    5. The Software may only be used for non-commercial research and may not
    be used for clinical care.

    6. Any publication by Recipient of research involving the Software shall
    cite the references listed below.

    REFERENCES
    Veraart, J.; Fieremans, E. & Novikov, D.S. Diffusion MRI noise mapping
    using random matrix theory Magn. Res. Med., 2016, early view, doi:
    10.1002/mrm.26059
    """

import os, sys
import numpy as np
import multiprocessing
from joblib import Parallel, delayed
import scipy.linalg
import warnings

warnings.filterwarnings('ignore')

class MP(object):
    def __init__(self, dwi, kernel, patchtype, patchsize, shrinkage, algorithm, crop):
        if dwi.ndim > 4:
            self.coil = True
        else:
            self.coil = False

        if shrinkage is None:
            self.shrink = 'threshold'
        else:
            self.shrink = shrinkage

        if algorithm is None:
            self.algo = 'jespersen'
        else:
            self.algo = algorithm

        if kernel is None:
            nvols = dwi.shape[-1]
            p = np.arange(3, nvols, 2)
            pf = np.where(p**3 >= nvols)[0][0]
            defaultKernel = p[pf]
            kernel = np.array([defaultKernel,defaultKernel,defaultKernel])
        else:
            kernel = np.array(kernel)

        kernel = kernel+np.mod(kernel,2)-1
        self.kernel = kernel.astype(int)
        k = self.kernel // 2
        
        if patchtype is None:
            self.patchtype = 'box'
        else:
            self.patchtype = patchtype

        if patchsize is None:
            self.patchsize = np.prod(kernel)
        else:
            self.patchsize = patchsize

        if self.patchtype == 'box':
            self.pos_distances = None
        elif self.patchtype == 'adaptive':
            if self.patchsize >= np.prod(kernel):
                print('Warning: selecting sane default adaptive patch size')
                self.patchsize = np.floor(np.prod(kernel) - 0.2*np.prod(kernel))
            elif self.patchsize <= dwi.shape[-1]:
                print('Warning: selecting sane default adaptive patch size')
                self.patchsize = dwi.shape[-1] + 1
            if not isinstance(self.patchsize, int):
                self.patchsize = self.patchsize.astype(int)
            pi, pj, pk = np.where(np.ones((self.kernel)))
            patchcoords = np.vstack((pi,pj,pk)).T
            self.pos_distances = np.sum((patchcoords - k)**2, axis=1)

        if crop is None:
            self.crop = 0
        else:
            self.crop = crop
                    
        if self.coil:
            pwidth = (k[0],k[0]), (k[1],k[1]), (k[2],k[2]), (0,0), (0,0)
        else:
            pwidth = (k[0],k[0]), (k[1],k[1]), (k[2],k[2]), (0,0)
        self.pwidth = pwidth
        self.dwi = np.pad(dwi, pad_width=pwidth, mode='wrap')

        print('Denoising data with parameters:')
        print('kernel     = ' + str(self.kernel))
        print('patch type = ' + str(self.patchtype))
        print('patch size = ' + str(self.patchsize))
        print('shrinkage  = ' + str(self.shrink))
        print('algorithm  = ' + str(self.algo))
            
    def box_patch(self, dwi, coords):
        # extracts a patch of size kx x ky x kz from the padded input at specified coords
        k = self.kernel // 2
        X = dwi[coords[0]-k[0]:coords[0]+k[0]+1, coords[1]-k[1]:coords[1]+k[1]+1, coords[2]-k[2]:coords[2]+k[2]+1, ...]
        return X

    def normalize(self, im):
        im = ((im - np.min(im)) * (1/(np.max(im) - np.min(im)) * 1.0))
        return im

    def refine(self, X_tmp, coords):
        refval = X_tmp[np.prod(self.kernel)//2,...]
        if self.coil:
            int_distances = 1/(X_tmp.shape[-1]*X_tmp.shape[-2]) * np.sum((X_tmp - refval)**2, axis=(1,2))
        else:
            int_distances = 1/(X_tmp.shape[-1]) * np.sum((X_tmp - refval)**2, axis=1)

        wdists = (self.pos_distances * int_distances)
        iidx =np.argpartition(wdists, self.patchsize)[:self.patchsize]
        minind = np.argmin(self.pos_distances[iidx])

        # debug = True
        # if debug:
        #     print('sup')
        #     k = self.kernel // 2
        #     sx, sy, sz, N = self.dwi.shape
        #     mask = np.zeros((sx,sy,sz))
        #     mask[coords[0]-k[0]:coords[0]+k[0]+1, coords[1]-k[1]:coords[1]+k[1]+1, coords[2]-k[2]:coords[2]+k[2]+1] = 1
        #     iidx = bn.argpartition(wdists, self.patchsize)[:self.patchsize]
        #     pos_im = np.reshape(self.pos_distances,self.kernel)
        #     int_im = np.reshape(int_distances,self.kernel)
        #     w_im = np.reshape(wdists,self.kernel)
                
        #     pmask = mask.copy()
        #     pmask[coords[0]-k[0]:coords[0]+k[0]+1, coords[1]-k[1]:coords[1]+k[1]+1, coords[2]-k[2]:coords[2]+k[2]+1] = pos_im
        #     imask = mask.copy()
        #     imask[coords[0]-k[0]:coords[0]+k[0]+1, coords[1]-k[1]:coords[1]+k[1]+1, coords[2]-k[2]:coords[2]+k[2]+1] = int_im
        #     wmask = mask.copy()
        #     wmask[coords[0]-k[0]:coords[0]+k[0]+1, coords[1]-k[1]:coords[1]+k[1]+1, coords[2]-k[2]:coords[2]+k[2]+1] = w_im

        #     minw = np.partition(wdists, self.patchsize)[self.patchsize]
        #     omask = wmask.copy()
        #     omask[wmask>minw] = 1

        #     import matplotlib.pyplot as plt
        #     plt.imshow(np.squeeze(self.dwi[coords[0],:,:,55]))
        #     alphas = np.squeeze(mask[coords[0],:,:])
        #     plt.imshow(alphas, alpha=0.5*alphas)
        #     plt.show()

        #     plt.imshow(np.squeeze(self.dwi[coords[0],:,:,55]))
        #     pd = np.squeeze(pmask[coords[0],:,:])
        #     plt.imshow(1-pd, alpha=0.75*alphas)
        #     plt.show()

        #     plt.imshow(np.squeeze(self.dwi[coords[0],:,:,55]))
        #     id = np.squeeze(imask[coords[0],:,:])
        #     plt.imshow(1-id, alpha=0.75*alphas)
        #     plt.show()

        #     plt.imshow(np.squeeze(self.dwi[coords[0],:,:,55]))
        #     wd = np.squeeze(wmask[coords[0],:,:])
        #     plt.imshow(1-wd, alpha=0.75*alphas)
        #     plt.show()

        #     plt.imshow(np.squeeze(self.dwi[coords[0],:,:,55]))
        #     od = np.squeeze(omask[coords[0],:,:])
        #     plt.imshow(1-od,alpha=alphas*(1-od))
        #     plt.show()

        return iidx, minind

    def unpad(self, x, pad_width):
        slices = []
        for c in pad_width:
            e = None if c[1] == 0 else -c[1]
            slices.append(slice(c[0], e))
        return x[tuple(slices)]

    def padded_sampling(self, mask):
        # outputs a grid x, y, z of which coordinates to loop over when processing
        sx, sy, sz = mask.shape
        k = self.kernel // 2
        mask[:k[0],:,:] = 0
        mask[sx-k[0]:,:,:] = 0
        mask[:,:k[1],:] = 0
        mask[:,sy-k[1]:,:] = 0
        mask[:,:,:k[2]] = 0
        mask[:,:,sz-k[2]:] = 0
        self.mask = mask
        x, y, z = np.where(mask==1)
        return x.astype(int), y.astype(int), z.astype(int)

    def sampling(self, mask):
        # outputs a grid x, y, z of which coordinates to loop over when processing
        x, y, z = np.where(mask==1)
        return x.astype(int), y.astype(int), z.astype(int)

    def eig_shrink(self, vals, gamma):
        t = 1 + np.sqrt(gamma)
        s = np.zeros((vals.shape))
        x = vals[vals > t]
        s[vals > t] = np.sqrt((x**2 - gamma - 1)**2 - 4*gamma) / x
        return np.diag(s)
    
    def denoise(self, coords, max_retries=3):
        X = self.box_patch(self.dwi, coords)

        if self.coil:
            X = np.reshape(X,(np.prod(self.kernel), self.dwi.shape[-2], self.dwi.shape[-1]))
        else:
            X = np.reshape(X,(np.prod(self.kernel), self.dwi.shape[-1]))

        if self.patchtype == 'adaptive':
            Xn = self.normalize(X)
            nonlocalinds, minind = self.refine(Xn, coords)
            X = X[nonlocalinds,...]
        else:
            minind = np.prod(self.kernel)//2

        if self.coil:
            X = X.reshape(self.patchsize*self.dwi.shape[-2], self.dwi.shape[-1])
        
        M = X.shape[0]
        N = X.shape[1]
        Mp = np.min((M,N))
        Np = np.max((M,N))
        
        if M < N:  
            X = np.conj(X).T
        
        for retry in range(max_retries):
            try:
                u, vals, v = scipy.linalg.svd(X, full_matrices=False)
                break
            except np.linalg.LinAlgError:
                if retry == max_retries - 1:
                    raise
                print(f"SVD did not converge. Retrying {retry + 1}/{max_retries}...")
                X += np.random.normal(scale=1e-6, size=X.shape)
        # u,vals,v = scipy.linalg.svd(X, full_matrices=False)
        vals = (vals**2).astype('float32')

        order = np.argsort(vals)[::-1]
        u = u[:,order]
        v = v[:,order]
        vals = vals[order]

        tn = self.crop
        ptn = np.arange(0,Mp-tn).T
        p = np.arange(0,Mp).T

        csum = np.cumsum(vals[::-1])[::-1]
        if self.algo == 'veraart':
            sigmasq_1 = csum / ((Mp-p)*Np)
            range_mp = 4*np.sqrt((Mp-ptn)*(Np-tn))

        elif self.algo == 'cordero-grande':
            sigmasq_1 = csum / ((Mp-p)*(Np-p))
            range_mp = 4*np.sqrt((Mp-ptn)*(Np-ptn))

        elif self.algo == 'jespersen':
            sigmasq_1 = csum / ((Mp-p)*(Np-p))
            range_mp = 4*np.sqrt((Np-tn)*(Mp))

        range_data = vals[:Mp-tn] - vals[Mp-1]
        sigmasq_2 = range_data / range_mp

        t = np.where(sigmasq_2 < sigmasq_1[:Mp-tn])

        if t[0].any():
            t = t[0][0]
        else:
            t = Mp-tn-1
        
        sigma = np.sqrt(sigmasq_1[t])
        if (sigma == 0) or (not np.isfinite(sigma)):
            sigma = np.finfo(float).eps

        npars = t
        if self.shrink == 'threshold':
            vals[t:] = 0
            s = np.matrix(u) * np.diag(np.sqrt(vals)) * np.matrix(v)
        else:
            vals_norm = np.sqrt(vals)/(np.sqrt(Mp)*sigma)
            vals_frobnorm = (np.sqrt(Mp)*sigma) * self.eig_shrink(vals_norm, Np/Mp)
            s = np.matrix(u) * vals_frobnorm * np.matrix(v) 
        
        if M<N:
            s = np.conj(s).T

        # except:
        #     #print('Warning: no MP fit')
        #     s = X
        #     sigma = np.nan
        #     npars = np.nan
        #     if M<N:
        #         s = np.conj(s).T 

        if self.coil:
            signal = np.squeeze(s[minind::self.dwi.shape[-2], :])
        else:    
            signal = np.squeeze(s[minind, :]) 

        return signal, sigma, npars

    def process(self):
        sx, sy, sz, N = self.dwi.shape
        mask = np.ones((sx,sy,sz))

        x, y, z = self.padded_sampling(mask)
        xsize = int(x.size)
        coords = np.vstack((x,y,z))
        
        inputs = range(0, xsize)
        num_cores = multiprocessing.cpu_count()
        
        # # parallel
        signal, sigma, npars = zip(*Parallel(n_jobs=num_cores, prefer='processes')\
           (delayed(self.denoise)(coords[:,i]) for i in inputs))
        
    #   #  serial
        # k = self.kernel // 2
        # #for t in inputs:
        # crds = np.array([48+k[0], 38+k[1], 25+k[2]])
        # a, b, c = self.denoise(crds)
           
        # reconstruct original data matrix
        Sigma = np.zeros((sx, sy, sz))
        Npars = np.zeros((sx, sy, sz))
        Signal = np.zeros((sx, sy, sz, N), dtype=complex)
        for nn in inputs:
            Sigma[x[nn], y[nn], z[nn]] = sigma[nn]
            Npars[x[nn], y[nn], z[nn]] = npars[nn]
            Signal[x[nn], y[nn], z[nn], :] = signal[nn]

        Signal = self.unpad(Signal, self.pwidth)
        Npars = self.unpad(Npars, self.pwidth[:][:-1])
        Sigma = self.unpad(Sigma, self.pwidth[:][:-1])
        return Signal, Sigma, Npars

def denoise(img, kernel=None, patchtype=None, patchsize=None, shrinkage=None, algorithm=None, phase=None):

    import ants
    if phase is not None:
        mag = img.copy()
        img_phi = ants.image_read(phase).numpy()
        nii = ants.image_read(phase)

        minphi = np.min(img_phi)
        maxphi = np.max(img_phi)
        
        if (maxphi - minphi) > 2 * np.pi:
            print('rescaling phase from -pi to pi')
            phi = (img_phi - minphi) / (maxphi - minphi) * 2 * np.pi - np.pi
        
        kernel_phase = [15, 15, 1]
        
        print('phase denoising - square 2D patching')
        img = mag * np.exp(1j*phi)
        mp_phase = MP(img, kernel_phase, patchtype='box', patchsize=None, shrinkage='threshold', algorithm='cordero-grande', crop=phi.shape[-1]//2)
        img_dn1, sigma1, npars1 = mp_phase.process()
        phi_dn = np.angle(img_dn1)

        # out = ants.from_numpy(npars1, origin=nii.origin[:-1], spacing=nii.spacing[:-1], direction=nii.direction[:-1,:])
        # ants.image_write(out, 'phase_npars.nii')
        # img_phase = np.angle(img * np.exp(-1j*phi_dn))
        # out = ants.from_numpy(img_phase, origin=nii.origin, spacing=nii.spacing, direction=nii.direction)
        # ants.image_write(out, 'phase_dn.nii')
        out = ants.from_numpy(phi_dn, origin=nii.origin, spacing=nii.spacing, direction=nii.direction)
        ants.image_write(out, 'phase_dn.nii')
        # import pdb; pdb.set_trace()

        print('magnitude denoising - adaptive patching')
        img_np = np.real(img*np.exp(-1j*phi_dn))
        mp = MP(img_np, kernel, patchtype='adaptive', patchsize=None, shrinkage='frob', algorithm='jespersen', crop=0)
        Signal, Sigma, Npars = mp.process()
        
    else:
        zeroinds = np.where(img==0)
        img[zeroinds] = np.finfo(img.dtype).eps

        mp = MP(img.copy(), kernel, patchtype='adaptive', patchsize=None, shrinkage='frob', algorithm='cordero-grande', crop=0)
        Signal, Sigma, Npars = mp.process()
        Signal[zeroinds] = 0

    return abs(Signal), Sigma, Npars

# if __name__ == "__main__":
#     denoise()
    

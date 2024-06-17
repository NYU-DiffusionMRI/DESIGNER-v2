"""A Python class containing an implementation of MPPCA denoising for magnitude and/or phase images.
        
    Inputs are a 4D image with dimensions (X x Y x Z x N)
    
    Usage:
    import mpcomplex as mp
    imgdn, sigma, nparameters = mp.denoise(img, kernel=[5, 5, 5], shrinkage='frob', algorithm=veraart)
    
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
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
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

import gc
import numpy as np
from scipy.linalg import svd
from numpy_groupies import aggregate

class MP(object):
    def __init__(self, dwi, kernel, step, shrinkage, algorithm, crop, n_cores):
        from skimage.util.shape import view_as_windows

        self.n_cores = n_cores

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
            default_kernel = p[pf]
            kernel = np.array([default_kernel,default_kernel,default_kernel])
        else:
            kernel = np.array(kernel)

        # kernel = kernel+np.mod(kernel,2)-1
        self.kernel = kernel.astype(int)
        k = self.kernel // 2
        self.patchsize = np.prod(self.kernel)
                    
        pwidth = (k[0],k[0]), (k[1],k[1]), (k[2],k[2]), (0,0)
        self.pwidth = pwidth
        self.dwi = np.pad(dwi, pad_width=pwidth, mode='wrap')
        
        if step is None:
            self.step = self.kernel // 2 + 1
        else:
            self.step = step

        if crop is None:
            self.crop = 0
        elif crop is True:
            self.crop = self.dwi.shape[-1]//2
        else:
            self.crop = crop

        self.imsize = self.dwi.shape
        self.dwi = self.dwi.reshape(-1, self.imsize[-1])

        t = np.arange(self.dwi.shape[0]).reshape(self.imsize[:3])
        self.temp = view_as_windows(t, self.kernel, step=self.step).reshape(-1, np.prod(self.kernel))

        print('Denoising data with parameters:')
        print('kernel      = ' + str(self.kernel))
        print('step        = ' + str(self.step))
        print('shrinkage   = ' + str(self.shrink))
        print('algorithm   = ' + str(self.algo))
        gc.collect()
            

    def unpad(self, x, pad_width):
        slices = []
        for c in pad_width:
            e = None if c[1] == 0 else -c[1]
            slices.append(slice(c[0], e))
        return x[tuple(slices)]

    def eig_shrink(self, vals, gamma):
        """
        eigenvalue shrinkage
        """
        t = 1 + np.sqrt(gamma)
        s = np.zeros((vals.shape))
        x = vals[vals > t]
        s[vals > t] = np.sqrt((x**2 - gamma - 1)**2 - 4*gamma) / x
        return np.diag(s)
    
    def denoise(self, X):
        """
        Main function to denoise the data matrix X
        """

        M = X.shape[0]
        N = X.shape[1]
        Mp = np.min((M,N))
        Np = np.max((M,N))
        
        if M < N:  
            X = np.conj(X).T
     
        u,vals,v = svd(X, full_matrices=False)
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

        return s, sigma, npars #, nonlocalinds

    def patchav(self, wp, signal, temp):
        """
        weighted patch averaging
        """
        return aggregate(temp, wp*signal, func='sum', size=self.dwi.shape[0], dtype=self.dwi.dtype
            ).reshape(self.imsize[:3])

    def get_weights(self, temp):
        from skimage.util.shape import view_as_windows

        i,j,k = np.unravel_index(temp, self.imsize[:3])
        distance = (i-np.mean(i, axis=1, keepdims=True))**2 + \
                     (j-np.mean(j, axis=1, keepdims=True))**2 + \
                     (k-np.mean(k, axis=1, keepdims=True))**2
        count = np.histogram(temp, np.arange(self.dwi.shape[0]+1))[0]

        den = view_as_windows(self.patchav(
            1, distance.flatten(), temp.flatten()
            ), self.kernel, self.step).reshape(
            -1, np.prod(self.kernel))
        cbl = view_as_windows(count.reshape(self.imsize[:3]
            ), self.kernel, self.step).reshape(
            -1, np.prod(self.kernel))
    
        num = den - distance
        wp = num / (den * (cbl-1) + np.finfo(float).eps)
        wp[cbl==1] = 1 
        return wp.flatten()

    def im_reconstruct(self, wp, image):
        from joblib import Parallel, delayed
        image = np.array(image)

        if image.ndim == 1:
            image = np.tile(image, (self.patchsize, 1)).T
            rec_img = self.patchav(wp, image.flatten(), self.temp.flatten())
            rec_img = self.unpad(rec_img, self.pwidth[:][:-1])
        elif image.ndim == 3:
            # Process in chunks to manage memory usage
            S = Parallel(n_jobs=-3, prefer='threads')\
            (delayed(self.patchav)(
                wp, image[:,:,i].flatten(), self.temp.flatten()
                ) for i in range(self.imsize[-1]))
            rec_img = np.array(S).transpose(1,2,3,0)
            rec_img = self.unpad(rec_img, self.pwidth)

        gc.collect()
        return rec_img

    def process(self):
        from joblib import Parallel, delayed
        
        num_patches = self.temp.shape[0]
        num_vols = self.dwi.shape[1]

        # Preallocating arrays for signal, sigma, and npars
        signal = np.empty((num_patches, num_vols), dtype=complex)
        sigma = np.empty(num_patches, dtype=np.float32)
        npars = np.empty(num_patches, dtype=np.int32)

        
        results = Parallel(n_jobs=self.n_cores, prefer='processes')\
            (delayed(self.denoise)(self.dwi[self.temp[i,:],:]) for i in range(num_patches))
        
        signal, sigma, npars = zip(*results)


        print('...patch averaging...')
        wp = self.get_weights(self.temp)
        signal = self.im_reconstruct(wp, signal)
        sigma = self.im_reconstruct(wp, sigma)
        npars = self.im_reconstruct(wp, npars)

        gc.collect()

        return signal, sigma, npars

def denoise(img, kernel=None, step=None, shrinkage=None, algorithm=None, crop=0, phase=None, n_cores=-1):
    ''' 
    denoising of complex data is implimented as a 2-pass procedure
    first denoise the phase using a large 2D patch,
    then denoise the magnitude using user input arguments
    '''
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
        step_phase = [kernel_phase[0]//2-1, kernel_phase[1]//2-1, 1]
        
        print('phase denoising')
        img = mag*np.exp(1j*phi)
        mp_phase = MP(img, kernel_phase, step_phase, 'threshold', 'cordero-grande', phi.shape[-1]//2, n_cores)
        img_dn1, sigma1, npars1 = mp_phase.process()
        phi_dn = np.angle(img_dn1)
        
        # out = ants.from_numpy(abs(npars1), origin=nii.origin[:-1], spacing=nii.spacing[:-1], direction=nii.direction[:-1,:])
        # ants.image_write(out, 'phase_npars.nii')        
        # img_phase = np.angle(img * np.exp(-1j*phi_dn))
        # out = ants.from_numpy(img_phase, origin=nii.origin, spacing=nii.spacing, direction=nii.direction)
        # ants.image_write(out, 'phase_dn.nii')
        out = ants.from_numpy(phi_dn, origin=nii.origin, spacing=nii.spacing, direction=nii.direction)
        ants.image_write(out, 'phase_dn.nii')
        # import pdb; pdb.set_trace()

        print('magnitude denoising')
        img_np = np.real(img*np.exp(-1j*phi_dn))
        mp = MP(img_np, kernel, step, shrinkage, algorithm, crop, n_cores)
        Signal, Sigma, Npars = mp.process()
    else:
        zeroinds = np.where(img==0)
        img[zeroinds] = np.finfo(img.dtype).eps

        mp = MP(img, kernel, step, shrinkage, algorithm, crop, n_cores)
        Signal, Sigma, Npars = mp.process()
        Signal[zeroinds] = 0



    return Signal, Sigma, Npars

def generate_parser():
    import argparse
    parser = argparse.ArgumentParser(
        description='denoise from command line'
        )
    parser.add_argument(
        'input',
        type=str,
        help='input path to 4D magnitude nifti image'
    )
    parser.add_argument(
        'output',
        type=str,
        help="output filename"
    )
    parser.add_argument(
        '-noisemap',
        type=str,
        nargs='?',
        help="output noise map filename"
    )
    parser.add_argument(
        '-nparsmap',
        type=str,
        nargs='?',
        help="output number of parameters map filename"
    )
    parser.add_argument(
        '-phase',
        type=str,
        nargs='?',
        help='input path to 4D phase nifti image'
    )
    parser.add_argument(
        '-extent',
        type=str,
        nargs='?',
        default='5,5,5',
        help='kernel extent, specify as 3 comma separated values or as a single scalar. (default=5,5,5)'
    )
    parser.add_argument(
        '-shrinkage',
        type=str,
        nargs='?',
        default='frob',
        help='type of shrinkage. Threshold or frob. (default=frob)'
    )
    parser.add_argument(
        '-algorithm',
        type=str,
        nargs='?',
        default='jespersen',
        help='MP cutoff algorithm. veraart or coredero-grande or jespersen. (default=jespersen)'
    )

    return parser

def getargs(args):
    parser = generate_parser()
    args = parser.parse_args(args)
    return args

def main():
    '''
    if called from the commmand line, import ants for image io,
    next load magnitude nifti and if available, phase nifti,
    next call the denoising function, applying user arguments,
    next save outputs
    
    '''
    import sys
    import ants

    args = getargs(sys.argv[1:])

    img_mag_ants = ants.image_read(args.input)
    img_mag = img_mag_ants.numpy()
    if args.phase:
        img_phi = ants.image_read(args.phase).numpy()
    else:
        img_phi = None

    if args.extent:
        extent = np.array([int(i) for i in args.extent.split(',')])
    else:
        extent = args.extent
   
    (Signal, Sigma, Npars) = denoise(img_mag, kernel=extent, step=extent//2,  shrinkage=args.shrinkage, algorithm=args.algorithm, crop=0, phase=img_phi)
    
    Signal_ants = ants.from_numpy(abs(Signal), origin=img_mag_ants.origin, spacing=img_mag_ants.spacing, direction=img_mag_ants.direction)
    ants.image_write(Signal_ants, args.output)

    if args.noisemap:
        Sigma_ants = ants.from_numpy(Sigma, origin=img_mag_ants.origin[:-1], spacing=img_mag_ants.spacing[:-1], direction=img_mag_ants.direction[:-1,:])
        ants.image_write(Sigma_ants, args.noisemap)

    if args.nparsmap:
        Npars_ants = ants.from_numpy(Npars, origin=img_mag_ants.origin[:-1], spacing=img_mag_ants.spacing[:-1], direction=img_mag_ants.direction[:-1,:])
        ants.image_write(Npars_ants, args.nparsmap)


if __name__ == "__main__":
    main()
    

#!/usr/bin/env python
# -*- coding: utf-8 -*-
# -------------------------------------------------------------------------
# Copyright (c) 2022, Rafael Neto Henriques, Benjamin Ades-Aron
#
# Developers : Rafael Neto Henriques (rafaelnh21@gmail.com)
#              Benjamin Ades-Aron (benjamin.ades-aron@nyulangone.org)
# -------------------------------------------------------------------------
# Adapted implementation of the gibbs removal procedure suggested by:
#
# Kellner E, Dhital B, Kiselev VG, Reisert M. Gibbs-ringing artifact removal
# based on local subvoxel-shifts. Magn Reson Med. 2015 doi: 10.1002/mrm.26054.
#
# Lee, H.H., Novikov, D.S., Fieremans, E., 2021. Removal of partial 
# Fourier‐induced Gibbs (RPG) ringing artifacts in MRI. 
# Magn Reson in Med. 86(5), pp.2733-2750.
#
# Full adaption of the original fully-sampled code is described in chapter 3
# of thesis:
#
# Neto Henriques, R., 2018. Advanced Methods for Diffusion MRI Data
# Analysis and their Application to the Healthy Ageing Brain
# (Doctoral thesis). https://doi.org/10.17863/CAM.29356
# -------------------------------------------------------------------------

import numpy as np
from tqdm import tqdm
from functools import partial
from multiprocessing import Pool
import scipy.fft as fft
from joblib import Parallel, delayed

class RPG(object):

    def _image_tv(self, x, axis=0, n_points=3):
        """ 
        Computes total variation (TV) of matrix x accross a given axis and
        along two directions.

        Parameters
        ----------
        x : 2D ndarray
            matrix x
        axis : int (0 or 1)
            Axis which TV will be calculated. Default a is set to 0.
        n_points : int
            Number of points to be included in TV calculation.

        Returns
        -------
        ptv : 2D ndarray
            Total variation calculated from the right neighbours of each point
        ntv : 2D ndarray
            Total variation calculated from the left neighbours of each point
        """
        xs = x.copy() if axis else x.T.copy()

        # Add copies of the data so that data extreme points are also analysed
        xs = np.concatenate((xs[:, (-n_points-1):], xs, xs[:, 0:(n_points+1)]),
                            axis=1)

        ptv = np.absolute(xs[:, (n_points+1):(-n_points-1)] -
                        xs[:, (n_points+2):(-n_points)])
        ntv = np.absolute(xs[:, (n_points+1):(-n_points-1)] -
                        xs[:, (n_points):(-n_points-2)])
        for n in range(1, n_points):
            ptv = ptv + np.absolute(xs[:, (n_points+1+n):(-n_points-1+n)] -
                                    xs[:, (n_points+2+n):(-n_points+n)])
            ntv = ntv + np.absolute(xs[:, (n_points+1-n):(-n_points-1-n)] -
                                    xs[:, (n_points-n):(-n_points-2-n)])

        if axis:
            return ptv, ntv
        else:
            return ptv.T, ntv.T

    def nn_resample(self, img, shape):
        """
        Computes the nearest neighbor resampling of a 2D image to a new grid shape

        Parameters
        ----------
        img : 2D ndarray
        shape : tuple of ints 
            Dimensions of the new grid img will be sampled to

        Returns
        -------
        img : 2D ndarray
            New image sampled to updated grid shape
        """
        def per_axis(in_sz, out_sz):
            ratio = 0.5 * in_sz / out_sz
            return np.round(np.linspace(ratio - 0.5, in_sz - ratio - 0.5, num=out_sz)).astype(int)

        return img[per_axis(img.shape[0], shape[0])[:, None],
                per_axis(img.shape[1], shape[1])]


    def _gibbs_removal_1d(self, x, axis=0, n_points=3):
        """ 
        Suppresses Gibbs ringing along a given axis using fourier sub-shifts.

        Parameters
        ----------
        x : 2D ndarray
            Matrix x.
        axis : int (0 or 1)
            Axis in which Gibbs oscillations will be suppressed. Default is set
            to 0
        n_points : int, optional
            Number of neighbours to access local TV (see note). Default is set to
            3.

        Returns
        -------
        xc : 2D ndarray
            Matrix with suppressed Gibbs oscillations along the given axis.

        Note
        ----
        This function suppresses the effects of Gibbs oscillations based on the
        analysis of local total variation (TV). Although artefact correction is
        done based on two adjanced points for each voxel, total variation should be
        accessed in a larger range of neighbours. The number of the neighbours to
        be considered in TV calculation can be adjusted using parameter n_points.
        """
        ssamp = np.linspace(0.01, .6, num=41)

        xs = x.copy() if axis else x.T.copy()
        
        # TV for shift zero (baseline)
        tvr, tvl = self._image_tv(xs, axis=1, n_points=n_points)
        tvp = np.minimum(tvr, tvl)
        tvn = tvp.copy()

        # Find optimal shift for gibbs removal
        isp = xs.copy()
        isn = xs.copy()
        sp = np.zeros(xs.shape)
        sn = np.zeros(xs.shape)
        N = xs.shape[1]
        c = fft.fft(xs, axis=1, norm='ortho', n=N)
        k = fft.fftfreq(N, 1 / (2.0j * np.pi))
        for s in ssamp:
            ks = k * s
            # Access positive shift for given s
            img_p = abs(fft.ifft(c * np.exp(ks), axis=1, norm='ortho', n=N))
            tvsr, tvsl = self._image_tv(img_p, axis=1, n_points=n_points)
            tvs_p = np.minimum(tvsr, tvsl)

            # Access negative shift for given s
            img_n = abs(fft.ifft(c * np.exp(-ks), axis=1, norm='ortho', n=N))
            tvsr, tvsl = self._image_tv(img_n, axis=1, n_points=n_points)
            tvs_n = np.minimum(tvsr, tvsl)

            # Update positive shift params
            isp[tvp > tvs_p] = img_p[tvp > tvs_p]
            sp[tvp > tvs_p] = s
            tvp[tvp > tvs_p] = tvs_p[tvp > tvs_p]

            # Update negative shift params
            isn[tvn > tvs_n] = img_n[tvn > tvs_n]
            sn[tvn > tvs_n] = s
            tvn[tvn > tvs_n] = tvs_n[tvn > tvs_n]

        # check non-zero sub-voxel shifts
        idx = np.nonzero(sp + sn)

        # use positive and negative optimal sub-voxel shifts to interpolate to
        # original grid points

        xs[idx] = (isp[idx] - isn[idx])/(sp[idx] + sn[idx])*sn[idx] + isn[idx]

        # import scipy.interpolate as sci
        # import matplotlib.pyplot as plt
        # plt.imshow(img_p); plt.colorbar(); plt.show(); 
        # # #xs_ = sci.interp2d(idx[0], idx[1], vals[idx], kind='cubic')
        # import pdb;
        # pdb.set_trace()

        return xs if axis else xs.T


    def _weights(self, shape):
        """ 
        Computes the weights necessary to combine two images processed by
        the 1D Gibbs removal procedure along two different axes [1]_.

        Parameters
        ----------
        shape : tuple
            shape of the image

        Returns
        -------
        G0 : 2D ndarray
            Weights for the image corrected along axis 0.
        G1 : 2D ndarray
            Weights for the image corrected along axis 1.

        References
        ----------
        .. [1] Kellner E, Dhital B, Kiselev VG, Reisert M. Gibbs-ringing artifact
            removal based on local subvoxel-shifts. Magn Reson Med. 2016
            doi: 10.1002/mrm.26054.
        """

        G0 = np.zeros(shape)
        G1 = np.zeros(shape)
        k0 = np.linspace(-np.pi, np.pi, num=shape[0])
        k1 = np.linspace(-np.pi, np.pi, num=shape[1])

        # Middle points
        K1, K0 = np.meshgrid(k1[1:-1], k0[1:-1])
        cosk0 = 1.0 + np.cos(K0)
        cosk1 = 1.0 + np.cos(K1)
        G1[1:-1, 1:-1] = cosk0 / (cosk0+cosk1)
        G0[1:-1, 1:-1] = cosk1 / (cosk0+cosk1)

        # Boundaries
        G1[1:-1, 0] = G1[1:-1, -1] = 1
        G1[0, 0] = G1[-1, -1] = G1[0, -1] = G1[-1, 0] = 1/2
        G0[0, 1:-1] = G0[-1, 1:-1] = 1
        G0[0, 0] = G0[-1, -1] = G0[0, -1] = G0[-1, 0] = 1/2

        return G0, G1


    def _gibbs_removal_2d(self, image, n_points=3, G0=None, G1=None):
        """ 
        Suppress Gibbs ringing of a 2D image.

        Parameters
        ----------
        image : 2D ndarray
            Matrix cotaining the 2D image.
        n_points : int, optional
            Number of neighbours to access local TV (see note). Default is
            set to 3.
        G0 : 2D ndarray, optional.
            Weights for the image corrected along axis 0. If not given, the
            function estimates them using function :func:`_weights`
        G1 : 2D ndarray
            Weights for the image corrected along axis 1. If not given, the
            function estimates them using function :func:`_weights`

        Returns
        -------
        imagec : 2D ndarray
            Matrix with Gibbs oscillations reduced along axis a.

        References
        ----------
        Please cite the following articles
        .. [1] Neto Henriques, R., 2018. Advanced Methods for Diffusion MRI Data
            Analysis and their Application to the Healthy Ageing Brain
            (Doctoral thesis). https://doi.org/10.17863/CAM.29356
        .. [2] Kellner E, Dhital B, Kiselev VG, Reisert M. Gibbs-ringing artifact
            removal based on local subvoxel-shifts. Magn Reson Med. 2016
            doi: 10.1002/mrm.26054.
        """
        if np.any(G0) is None or np.any(G1) is None:
            G0, G1 = self._weights(image.shape)

        img_c1 = self._gibbs_removal_1d(image, axis=1, n_points=n_points)
        img_c0 = self._gibbs_removal_1d(image, axis=0, n_points=n_points)

        C1 = fft.fft2(img_c1, norm='ortho')
        C0 = fft.fft2(img_c0, norm='ortho')
        imagec = abs(fft.ifft2(fft.fftshift(C1)*G1 + fft.fftshift(C0)*G0, norm='ortho'))
       
        return imagec

    def _gibbs_removal_2d_ff(self, image, n_points=3, G0=None, G1=None):
        if np.any(G0) is None or np.any(G1) is None:
            G0, G1 = self._weights(image.shape)

        img_f1 = abs(fft.ifft2(fft.fftshift(fft.fft2(image, norm='ortho'))*G1, norm='ortho'))
        img_f0 = abs(fft.ifft2(fft.fftshift(fft.fft2(image, norm='ortho'))*G0, norm='ortho'))

        img_c1 = self._gibbs_removal_1d(img_f1, axis=1, n_points=n_points)
        img_c0 = self._gibbs_removal_1d(img_f0, axis=0, n_points=n_points)

        imagec = abs(img_c1 + img_c0)

        return imagec

    def _gibbs_removal_2d_y(self, image, n_points=3):
        """ 
        Suppress Gibbs ringing of a 2D image along a single dimension.
        For use with partial fourier factors 5/8 and 7/8.

        Parameters
        ----------
        image : 2D ndarray
            Matrix cotaining the 2D image.
        n_points : int, optional
            Number of neighbours to access local TV (see note). Default is
            set to 3.

        Returns
        -------
        imagec : 2D ndarray
            Matrix with Gibbs oscillations reduced only in 1 direction.

        References
        ----------
        Please cite the following articles
        .. [1] Neto Henriques, R., 2018. Advanced Methods for Diffusion MRI Data
            Analysis and their Application to the Healthy Ageing Brain
            (Doctoral thesis). https://doi.org/10.17863/CAM.29356
        .. [2] Kellner E, Dhital B, Kiselev VG, Reisert M. Gibbs-ringing artifact
            removal based on local subvoxel-shifts. Magn Reson Med. 2016
            doi: 10.1002/mrm.26054.
        .. [3] Lee, H.H., Novikov, D.S., Fieremans, E., 2021. Removal of partial 
            Fourier-induced Gibbs (RPG) ringing artifacts in MRI. 
            Magn Reson in Med. 86(5), pp.2733-2750.
        """

        # remove gibbs along axis 1
        imagec = self._gibbs_removal_1d(image, axis=0, n_points=n_points)

        return imagec

    def _gibbs_removal_2d_y_2(self, image, n_points=3, G0=None, G1=None):
        """ 
        Suppress Gibbs ringing of a 2D image with 6/8 partial fourier.

        Parameters
        ----------
        image : 2D ndarray
            Matrix cotaining the 2D image.
        n_points : int, optional
            Number of neighbours to access local TV (see note). Default is
            set to 3.
        G0 : 2D ndarray, optional.
            Weights for the image corrected along axis 0. If not given, the
            function estimates them using function :func:`_weights`
        G1 : 2D ndarray
            Weights for the image corrected along axis 1. If not given, the
            function estimates them using function :func:`_weights`

        Returns
        -------
        imagec : 2D ndarray
            Matrix with Gibbs oscillations reduced along axis a.

        References
        ----------
        Please cite the following articles
        .. [1] Neto Henriques, R., 2018. Advanced Methods for Diffusion MRI Data
            Analysis and their Application to the Healthy Ageing Brain
            (Doctoral thesis). https://doi.org/10.17863/CAM.29356
        .. [2] Kellner E, Dhital B, Kiselev VG, Reisert M. Gibbs-ringing artifact
            removal based on local subvoxel-shifts. Magn Reson Med. 2016
            doi: 10.1002/mrm.26054.
        .. [3] Lee, H.H., Novikov, D.S., Fieremans, E., 2021. Removal of partial 
            Fourier-induced Gibbs (RPG) ringing artifacts in MRI. 
            Magn Reson in Med. 86(5), pp.2733-2750.
        """
        if np.any(G0) is None or np.any(G1) is None:
            G0, G1 = self._weights(image.shape)

        img_c1 = self._gibbs_removal_1d(image, axis=1, n_points=n_points)
        img_c0 = self._gibbs_removal_1d(image, axis=0, n_points=n_points)

        img_c00 = img_c0[::2,:]
        img_c01 = img_c0[1::2,:]

        img_c10 = self._gibbs_removal_1d(img_c00, axis=0, n_points=n_points)
        img_c11 = self._gibbs_removal_1d(img_c01, axis=0, n_points=n_points)
        img_c = np.zeros_like(image)

        img_c[::2,:] = img_c10
        img_c[1::2,:] = img_c11

        C1 = fft.fft2(img_c)
        C0 = fft.fft2(img_c1)
        imagec = abs(fft.ifft2(fft.fftshift(C1)*G1 + fft.fftshift(C0)*G0))

        return imagec

    def _parallel_pf_gibbs_removal(self, vol, shap, pf_fact, n_points, G0, G1, testff):
        """ 
        Wrapper function for parallel processing
        """

        if pf_fact == 5/8:
            scale = 4
            for sc in range(scale):
                vol[sc::scale,:] = self._gibbs_removal_2d(vol[sc::scale,:], n_points=n_points)
            vol = self._gibbs_removal_2d_y(vol, n_points=n_points)

        elif pf_fact == 6/8:
            vol = self._gibbs_removal_2d_y_2(vol, n_points=n_points)
            #vol = self._gibbs_removal_2d(vol, n_points=n_points, G0=G0, G1=G1)
            # from hong-hsi's supplemental material, the following will have better perfamnce
            # but oversmooth:
            # scale = 2
            # for sc in range(scale):
            #     vol[sc::scale,:] = _gibbs_removal_2d(vol[sc::scale,:], n_points=n_points)
        
        elif pf_fact == 7/8:
            scale = 4
            vol = self.nn_resample(vol, (shap[1]*3, shap[2]))
            for sc in range(scale):
                vol[sc::scale,:] = self._gibbs_removal_2d(vol[sc::scale,:], n_points=n_points)
            vol = self.nn_resample(vol, (shap[1], shap[2]))
            vol = self._gibbs_removal_2d_y(vol, n_points=n_points)

        
        elif pf_fact == 1 and testff:
            vol = self._gibbs_removal_2d_ff(vol, n_points=n_points, G0=G0, G1=G1)
        elif pf_fact == 1 and not testff:
            vol = self._gibbs_removal_2d(vol, n_points=n_points, G0=G0, G1=G1)

        return vol

    def gibbs_removal(self, vol, slice_axis=2, n_points=3, pe_dim=1, pf_fact=1, nproc=1, testff=False):
        """ 
        Suppresses Gibbs ringing artefacts of images volumes.

        Parameters
        ----------
        vol : ndarray ([X, Y]), ([X, Y, Z]) or ([X, Y, Z, g])
            Matrix containing one volume (3D) or multiple (4D) volumes of images.
        slice_axis : int (0, 1, or 2)
            Data axis corresponding to the number of acquired slices. Default is
            set to the third axis
        n_points : int, optional
            Number of neighbour points to access local TV (see note). Default is
            set to 3.

        Returns
        -------
        vol : ndarray ([X, Y]), ([X, Y, Z]) or ([X, Y, Z, g])
            Matrix containing one volume (3D) or multiple (4D) volumes of corrected
            images.

        Notes
        -----
        For 4D matrix last element should always correspond to the number of
        diffusion gradient directions.

        References
        ----------
        Please cite the following articles
        .. [1] Neto Henriques, R., 2018. Advanced Methods for Diffusion MRI Data
            Analysis and their Application to the Healthy Ageing Brain
            (Doctoral thesis). https://doi.org/10.17863/CAM.29356
        .. [2] Kellner E, Dhital B, Kiselev VG, Reisert M. Gibbs-ringing artifact
            removal based on local subvoxel-shifts. Magn Reson Med. 2016
            doi: 10.1002/mrm.26054.
        .. [3] Lee, H.H., Novikov, D.S., Fieremans, E., 2021. Removal of partial 
            Fourier-induced Gibbs (RPG) ringing artifacts in MRI. 
            Magn Reson in Med. 86(5), pp.2733-2750.
        """
        nd = vol.ndim
        
        # check the axis corresponding to different slices
        # 1) This axis cannot be larger than 2
        if slice_axis > 2:
            raise ValueError("Different slices have to be organized along" +
                            "one of the 3 first matrix dimensions")
        
        # check matrix dimension
        elif nd == 2:
            vol = vol[None,...]
        elif nd == 3:
            vol = np.moveaxis(vol, slice_axis, 0)
        elif nd == 4:
            vol = np.moveaxis(vol, (slice_axis, 3), (0,1))

        if nd == 4:
            inishap = vol.shape
            vol = vol.reshape((inishap[0] * inishap[1], inishap[2],  inishap[3]))
        
        if nd > 4:
            raise ValueError("Data have to be a 4D, 3D or 2D matrix")
        elif nd < 2:
            raise ValueError("Data is not an image")

        # Produce weigthing functions for 2D Gibbs removal
        shap = vol.shape
        G0, G1 = self._weights(shap[1:])

        # check that phase encoding axis is either 0 or 1
        if pe_dim > 2:
            raise ValueError("Phase encoding axis must be either 0 (LR) or 1 (AP)")

        # if phase incoding direction is along axis 0, swap first two axes
        if pe_dim == 1:
            vol = np.swapaxes(vol, 1, 2)

        pf_opts = [5/8, 6/8, 7/8, 1]
        if not pf_fact in pf_opts:
            raise ValueError("RPG degibbs only supports partial fourier factors of 5/8, 6/8, 7/8, and 1")

        vol = vol.copy()

        inputs = tqdm(range(vol.shape[0]))

        vol = (Parallel(n_jobs=nproc, prefer='processes')
            (delayed(self._parallel_pf_gibbs_removal)(
                vol[i,...], shap=shap, pf_fact=pf_fact, n_points=n_points, G0=G0, G1=G1, testff=testff
            ) for i in inputs))
        vol = np.array(vol)

        # pool = Pool(nproc)
        # partial_func = partial(
        #     self._parallel_pf_gibbs_removal, 
        #     shap=shap, pf_fact=pf_fact, n_points=n_points, G0=G0, G1=G1, testff=testff
        # )
        # vol[:, :, :] = pool.map(partial_func, vol)
        # pool.close()
        # pool.join()
    
        #vol[:, :, :] = self._parallel_pf_gibbs_removal(vol[0, :, :], shap=shap, pf_fact=pf_fact, n_points=n_points, G0=G0, G1=G1, testff=testff)

        # Reshape data to original format
        if pe_dim == 1:
            vol = np.swapaxes(vol, 1, 2)
        if nd == 2:
            vol = vol.squeeze()
        if nd == 3:
            vol = np.moveaxis(vol, 0, slice_axis)
        if nd == 4:
            vol = vol.reshape(inishap)
            vol = np.moveaxis(vol, (0, 1), (slice_axis, 3))

        return vol

def gibbs_removal(image, slice_axis=2, n_points=3, pe_dim=1, pf_fact=1, nproc=1, testff=False):
    rpg = RPG()
    return rpg.gibbs_removal(image, slice_axis, n_points, pe_dim, pf_fact, nproc, testff)


def fft2_mri(img):
    n = img.shape[0]*img.shape[1]
    return np.fft.fftshift(np.fft.fft2(np.fft.fftshift(img))) / 2

def ifft2_mri(img):
    n = img.shape[0]*img.shape[1]
    return np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(img))) / 2

def make_gibbs_phantom(img, size, scale, pf, sigma, pe_dir=1):
    nps = np.ceil(size/scale)
    pk = fft2_mri(img)
    scale = (1-scale)/2
    pk = pk[round(img.shape[0]*scale):round(img.shape[0]*(1-scale)), 
            round(img.shape[1]*scale):round(img.shape[1]*(1-scale))]
    pk = pk[:size,:size]
    pk = pk + sigma * np.random.randn(size,size) + 1j * sigma * np.random.randn(size,size)
    p = ifft2_mri(pk)
    if pe_dir == 1:
        pk[:,round(size*pf):] = 0
    elif pe_dir == 0:
        pk[round(size*pf):,:] = 0
    p_pf = ifft2_mri(pk)
    return p, p_pf

if __name__ == '__main__':
    from phantominator import shepp_logan
    import matplotlib.pyplot as plt
    import scipy.io as sio

    sz = 128
    scale = 0.5
    pf = 6/8
    sigma = 0.02
    phantom = np.flip(shepp_logan(256))

    phantom_gibbs, phantom_gibbs_pf = make_gibbs_phantom(phantom, sz, scale, pf, sigma)
    sio.savemat('/Users/benaron/Desktop/gibbs_phantom.mat', {'phantom': phantom_gibbs})

    testff = True
    g1 = gibbs_removal(abs(phantom_gibbs), testff=testff, n_points=10)
    
    testff = False
    g2 = gibbs_removal(abs(phantom_gibbs), testff=testff, n_points=10)

    mat = sio.loadmat('/Volumes/Research/cbi05data/data1/Hamster/Ben/gibbs_phantom_hhdg.mat')
    g3 = mat['P_pf_dg']

    plt.subplot(3,3,1)
    plt.imshow(abs(phantom_gibbs), vmin=0, vmax=1)
    plt.title('gibbs phantom (no pf)'); plt.axis('off')
    plt.subplot(3,3,2)
    plt.imshow(abs(g1), vmin=0, vmax=1)
    plt.title('gibbs corrected (filter before unring1d)'); plt.axis('off')
    plt.subplot(3,3,3)
    plt.imshow(abs(g1)-abs(phantom_gibbs), vmin=-.25, vmax=.25)
    plt.title('residual'); plt.axis('off')
    plt.subplot(3,3,4)
    plt.imshow(abs(phantom_gibbs), vmin=0, vmax=1)
    plt.title('gibbs phantom (no pf)'); plt.axis('off')
    plt.subplot(3,3,5)
    plt.imshow(abs(g2), vmin=0, vmax=1)
    plt.title('gibbs corrected (filter after unring1d)'); plt.axis('off')
    plt.subplot(3,3,6)
    plt.imshow(abs(g2)-abs(phantom_gibbs), vmin=-.25, vmax=.25)
    plt.title('residual'); plt.axis('off')
    plt.subplot(3,3,7)
    plt.imshow(abs(phantom_gibbs), vmin=0, vmax=1)
    plt.title('gibbs phantom (no pf)'); plt.axis('off')
    plt.subplot(3,3,8)
    plt.imshow(g3, vmin=0, vmax=1)
    plt.title('gibbs corrected (Lee matlab code)'); plt.axis('off')
    plt.subplot(3,3,9)
    plt.imshow(abs(g3)-abs(phantom_gibbs), vmin=-.25, vmax=.25)
    plt.title('residual'); plt.axis('off')
    plt.show()




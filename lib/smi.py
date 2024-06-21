"""
SMI: Standard Model Imaging class
"""
import os
import sys
from math import comb
import numpy as np
import itertools

import scipy.ndimage as snd
import scipy.cluster.hierarchy as sch
import scipy.special as scp
import scipy.linalg as scl
import scipy.io as sio
import scipy.optimize as sco

import warnings

class SMI(object):

    def __init__(self, bval, bvec, beta=None, echo_time=None, merge_distance=None, cs_phase=1, flag_fit_fodf=0, 
        flag_rectify_fodf=0, compartments=None, n_levels=10, l_max=None, rotinv_lmax=None, 
        noise_bias=None, training_bounds=None, training_prior=None, n_training=1e5, 
        l_max_training = None, seed=42):
        """
        Setting some default values and initialization required for class functions
        """

        if not cs_phase:
            self.cs_phase = 1
        else:
            self.cs_phase = cs_phase

        self.merge_distance = merge_distance
        self.flag_fit_fodf = flag_fit_fodf
        self.flag_rectify_fodf = flag_rectify_fodf
        self.n_levels = n_levels
        self.l_max = l_max
        self.project_forward_SH = False
        
        self.rotinv_lmax = rotinv_lmax
        self.n_training = int(n_training)

        self.set_compartments(compartments)

        if not l_max_training:
            self.l_max_training = 6
        else:
            self.l_max_training = l_max_training
        
        if training_bounds is not None and training_prior is not None:
            self.prior = training_prior
            if not self.fit_T2:
                self.prior = np.hstack(
                    (self.prior[:,:4], 100 * np.ones((self.prior.shape[0], 2)), self.prior[:,5:])
                    )  
        else:
            if training_bounds is None:
                self.lb_training = [0.05, 1, 1, 0.1,   0,  50,  50, 0.05]
                self.ub_training = [0.95, 3, 3, 1.2, 0.5, 150, 120, 0.99]
            else:
                self.lb_training = training_bounds[0,:]
                self.ub_training = training_bounds[1,:]
            self.flag_get_priors = True
            self.flag_Depar_Deperp = False

        if not noise_bias:
            self.flag_rician_bias = False
        elif noise_bias == 'rician':
            self.flag_rician_bias = True

        if seed is not None:
            self.seed = seed
            np.random.seed(self.seed)
        else:
            self.seed = None

        # set up required inputs
        self.set_bvals_bvecs(bval, bvec)
        self.set_bshape(beta)
        self.set_echotime(echo_time)

    def set_mask(self, mask):
        """
        Function to add a brain mask to the class if the mask was input
        If a mask was not input the full image volume is processed.
        """
        if mask is not None:
            self.mask = mask
        else:
            self.mask = np.ones(self.dwishape[:3])

    def set_sigma(self, sigma, dwi):
        """
        Function to add noise level sigma to the class if sigma was input
        If sigma was not input, it is esimtated from b0 data.
        """
        if sigma is None:
            if self.echo_time is None:
                b0_for_sigma_flag = self.b < 0.1
            else:
                b0_for_sigma_flag = (
                    (self.b < 0.1) & (self.echo_time == np.min(self.echo_time))
                    )
        else:
            self.sigma = sigma
            b0_for_sigma_flag = None

        if b0_for_sigma_flag is not None:
            std_dev = dwi[:,:,:,b0_for_sigma_flag].std(axis=-1)
            self.sigma = snd.gaussian_filter(std_dev, sigma=1)
            warnings.warn('sigma was not an input (not recommended)')
        
    def set_bvals_bvecs(self, bval, bvec):
        """
        Function to add gradient table to the class. Gradient magnitude and direction
        are required inputs.
        """
        if not bval.any() or not bvec.any():
            raise Exception('b and dirs are compulsory arguments for SMI to run')
        else:
            self.b = bval
            self.dirs = bvec

            if np.max(self.b) < 500:
                self.flag_microstructure_units = 1
            else:
                self.flag_microstructure_units = 0
                

    def set_bshape(self, beta):
        if beta is None:
            self.beta = np.zeros(self.b.shape) + 1
        else:
            self.beta = beta
        
        try:
            if self.flag_microstructure_units:
                self.beta[self.b < 0.05] = 1
            else:
                self.beta[self.b < 50] = 1
        except:
            raise Exception('mismatch between input bshape and gradients, please check')

        if not self.merge_distance:
            if self.flag_microstructure_units:
                self.merge_distance = 0.1
            else:
                self.merge_distance = 100


    def set_echotime(self, echo_time):
        if echo_time is None:
            self.echo_time = np.zeros(self.b.shape)
        else:
            self.echo_time = echo_time

        if (echo_time is None) or (len(np.unique(echo_time)) == 1):
            self.fit_T2 = 0
        else:
            self.fit_T2 = 1
        

    def set_lmax_rotinv_smi(self):
        """
        Function to set a default l-max for rotational invariant and smi fits.
        """
        if not self.l_max:
            self.l_max = self.get_default_lmax()

        if not self.rotinv_lmax and np.max(self.l_max) >= 4:
            self.rotinv_lmax = 4
        elif not self.rotinv_lmax:
            self.rotinv_lmax = 2

    def set_compartments(self, compartments):
        """
        Function to choose which compartment fractions should be estimated.
        """
        if not compartments:
            compartments = ['EAS', 'IAS', 'FW']
            
        self.flag_compartments = [0,0,0,0]
        if 'IAS' in compartments:
            self.flag_compartments[0] = 1
        if 'EAS' in compartments:
            self.flag_compartments[1] = 1
        if 'FW' in compartments:
            self.flag_compartments[2] = 1
        if 'DOT' in compartments:
            self.flag_compartments[3] = 1
            raise Exception('Current version does not support a DOT compartment (update should be ready soon)')

    def set_priors(self):
        """
        Function to set up priors for SM fitting
        """
        if self.flag_get_priors:
            self.prior = self.get_uniformly_distributed_SM_prior()

        if (self.rotinv_lmax == 2 and self.prior.shape[1] < 8) \
            or (self.rotinv_lmax == 4 and self.prior.shape[1] < 9) \
            or (self.rotinv_lmax == 6 and self.prior.shape[1] < 10):
            raise Exception('Inconsistency between desired Lmax for ML fitting and prior distribution')


    def l2norm(self, x):
        """
        Computes the L2 norm of x
        """
        return np.sqrt(np.sum(x**2, axis=0))

    def get_uniformly_distributed_SM_prior(self):
        """
        function to identify lower and uppoer bounds on SM parameters and set up
        prior for each parameter
        Inputs:
        -------
        training bounds: (inherited from self)

        Outputs:
        --------
        Prior: ND array
        """
        # water fractions

        f = (np.random.uniform(size=(5 * self.n_training)) * 
            (self.ub_training[0] - self.lb_training[0]) + self.lb_training[0]
            )
        f_fw = (np.random.uniform(size=(5 * self.n_training)) * 
            (self.ub_training[4] - self.lb_training[4]) + self.lb_training[4]
            )
        keep = 1 > (f + f_fw)
        f = f[keep]
        f = f[:self.n_training]
        f_fw = f_fw[keep]
        f_fw = f_fw[:self.n_training]
        f_extra = 1 - f - f_fw

        # diffusivities
        Da = (np.random.uniform(size=self.n_training) * 
            (self.ub_training[1] - self.lb_training[1]) + self.lb_training[1]
            )
        if self.flag_Depar_Deperp:
            Depar = (np.random.uniform(size=(5 * self.n_training)) * 
                (self.ub_training[2] - self.lb_training[2]) + self.lb_training[2]
                )
            Deperp = (np.random.uniform(size=(5 * self.n_training)) * 
                (self.ub_training[3] - self.lb_training[3]) + self.lb_training[3]
                )
            keep = Depar > Deperp
            Depar = Depar[keep]
            Deperp = Deperp[keep]
        else:
            Depar = (np.random.uniform(size=self.n_training) * 
                (self.ub_training[2] - self.lb_training[2]) + self.lb_training[2]
                )
            Deperp = (np.random.uniform(size=self.n_training) * 
                (self.ub_training[3] - self.lb_training[3]) + self.lb_training[3]
                )

        # T2 values
        T2a = (np.random.uniform(size=self.n_training) * 
            (self.ub_training[5] - self.lb_training[5]) + self.lb_training[5]
            )
        T2e = (np.random.uniform(size=self.n_training) * 
            (self.ub_training[6] - self.lb_training[6]) + self.lb_training[6]
            )

        # p2m
        p2 = (np.random.uniform(size=self.n_training) * 
            (self.ub_training[7] - self.lb_training[7]) + self.lb_training[7]
            )
        px = p2.copy()
        dims = 2 * np.arange(0, self.l_max_training + 2, 2)
        plm = np.array([]).reshape(0, self.n_training)
        for i in range(len(dims) - 1):
            dim = dims[i + 1] + 1
            x = np.random.standard_normal((dim, self.n_training))
            initial_norm = self.l2norm(x)
            x_normalised = x / np.tile(initial_norm, (dim, 1))
            x_newnorm = x_normalised * np.tile(px.T, (dim, 1))
            plm = np.vstack((plm, x_newnorm))
            px = np.random.uniform(size=self.n_training) * px * 0.9
        
        prior = np.vstack((f, Da, Depar, Deperp, f_fw, T2a, T2e, p2))

        # rotinvs
        dims = 2 * np.arange(2, self.rotinv_lmax + 2, 2)
        for i in range(len(dims) - 1):
            dim = dims[i]
            dim_next = dims[i + 1]
            px = self.l2norm(plm[dim:(dim+dim_next+1), :])
            prior = np.vstack((prior, px))

        return prior.T

    def vectorize(self, image, mask):
        """
        - If the input is 1D or 2D: unpatch it to 3D or 4D using a mask.
        - If the input is 3D or 4D: vectorize to 2D using a mask.
        """
        if mask is None:
            mask = np.ones((image.shape[0], image.shape[1], image.shape[2]))

        if image.ndim == 1:
            image = np.expand_dims(image, axis=0)

        if image.ndim == 2:
            n = image.shape[0]
            s = np.zeros((mask.shape[0], mask.shape[1], mask.shape[2], n))
            for i in range(0, n):
                dummy = np.zeros(mask.shape)
                dummy[mask.astype(bool)] = image[i,:]
                s[:,:,:,i] = dummy

        if image.ndim == 3:
            image = np.expand_dims(image, axis=-1)

        if image.ndim == 4:
            s = np.zeros((image.shape[-1], np.sum(mask).astype(int)))
            for i in range(0, image.shape[-1]):
                tmp = image[:,:,:,i]
                s[i,:] = tmp[mask.astype(bool)]
             
        return np.squeeze(s)

    def cart_to_sph(self, x, y, z):
        """
        Function to switch from cartesian to spherical coordinates
        Inputs:
        -------
        x, y, z

        Outputs:
        --------
        theta, phi, r
        """
        hxy = np.hypot(x, y)
        r = np.hypot(hxy, z)
        phi = np.arctan2(z, hxy)
        theta = np.arctan2(y, x)
        return theta, phi, r

    def get_even_sh(self, dirs, l_max):
        """
        if CS_phase=1, then the definition uses the Condon-Shortley phase factor
        of (-1)^m. Default is CS_phase=0 (so this factor is ommited)
            
        By: Santiago Coelho
        """

        def legendre(n, X) :
            res = []
            for m in range(int(n+1)):
                res.append(scp.lpmv(m, n, X))
            return np.array(res)

        if dirs.shape[1] != 3:
            dirs = dirs.T
        
        n_meas = dirs.shape[0]
        phi, theta, _ = self.cart_to_sph(dirs[:,0], dirs[:,1], dirs[:,2])
        theta = np.pi / 2 - theta

        l = np.arange(0, np.max(l_max) + 2, 2, dtype=int)
        l_all = np.array([])
        m_all = np.array([])
        for i in range(len(l)):
            l_all = np.hstack((l_all, l[i] * np.ones(2 * l[i] + 1)))
            m_all = np.hstack((m_all, np.arange(-l[i], l[i] + 1)))

        k_lm = (np.sqrt((2 * l_all + 1) / (4 * np.pi) * 
                scp.factorial(l_all - abs(m_all)) / scp.factorial(l_all + abs(m_all)))
                )

        extra_factor = np.ones(k_lm.shape)
        extra_factor[m_all != 0] = np.sqrt(2)
        if self.cs_phase:
            extra_factor = extra_factor * (-1) ** m_all
        
        p_l_in_cos_theta = np.zeros((len(l_all), n_meas))
        phi_term = np.zeros((len(l_all), n_meas))
        id_which_pl = np.zeros(len(l_all))
        for i in range(len(l_all)):
            all_pls = legendre(l_all[i], np.cos(theta))
            p_l_in_cos_theta[i,:] = all_pls[abs(m_all[i]).astype(int), :]
            id_which_pl[i] = abs(m_all[i])
            if m_all[i] > 0:
                phi_term[i,:] = np.cos(m_all[i] * phi)
            elif m_all[i] == 0:
                phi_term[i,:] = 1
            elif m_all[i] < 0:
                phi_term[i,:] = np.sin(-m_all[i] * phi)

        y_lm = np.tile(extra_factor, (n_meas, 1)).T * np.tile(k_lm, (n_meas, 1)).T * \
            phi_term * p_l_in_cos_theta

        return y_lm.T


    def fit2D4D_LLS_RealSphHarm_wSorting_norm_var(self, dwi, mask, rank1sl=True):
        """
        This function performs the linear fit of Slm on each shell using real spherical
        harmonics as defined in https://cs.dartmouth.edu/~wjarosz/publications/dissertation/appendixB.pdf
        
        The Condon-Shortley phase is used in this definition
        
        Inputs: 
        -------
        DWI: assumed to be either 2D [Nmeasurements x Nvoxels] or 4D [Nx x Ny x Nz x Nmeasurements]
        bval, bvec, beta, TE indicate the acquisition
        If beta or TE are constant throughout the acquisition, these can be left empty
        Lmax is the maximum order considered for the fitting of Slms
        (this can be a different order in all clusters but input order is the same as given by Group_dwi_in_shells_b_beta_TE)
        
        Outputs: 
        --------
        Slm: size [N_SH_coeffs x Nvoxels x Nclusters]
        Sl:  size [Lmax/2+1 x Nvoxels x Nclusters]

        The order of the Rotational Invariants is given by Group_dwi_in_shells_b_beta_TE
        and it is provided in table_4D_sorted
    
        By: Santiago Coelho        
        """

        sz_dwi = dwi.shape
        if dwi.ndim == 4:
            flag_4D = 1
        elif dwi.ndim == 2:
            flag_4D = 0 
        else:
            raise Exception('DWI must be a 2D or 4D array')

        if flag_4D:
            dwi_2D = self.vectorize(dwi, mask)
        else:
            dwi_2D = dwi

        n_voxels = dwi_2D.shape[1]
        n_meas = dwi_2D.shape[0]
       
        if len(self.b) != n_meas:
            raise Exception('bval should be a vector with the same length as the number of DWI')
        
        bb = self.table_4d
        l_max = self.l_max
        n_shells = bb.shape[1]

        if np.isscalar(l_max):
            l_max = l_max * np.ones(n_shells)

        l_max[bb[0,:] < 0.025] = 0
        l_max[abs(bb[1,:]) < 0.025] = 0

        n_sh_coeffs = (1 + l_max * (l_max + 3) / 2).astype(int)
        if self.dirs.shape[1] == n_meas:
            self.dirs = self.dirs.T

        l = np.arange(0, np.max(l_max) + 2, 2, dtype=int)
        l_all = np.array([])
        m_all = np.array([])
        for i in range(len(l)):
            l_all = np.hstack((l_all, l[i] * np.ones(2 * l[i] + 1)))
            m_all = np.hstack((m_all, np.arange(-l[i], l[i] + 1)))
        
        n_l_all = np.sqrt((2 * l_all + 1) * (4 * np.pi))
        n_l_unique = np.sqrt((2 * l + 1) * (4 * np.pi))

        y_lm_matrix = self.get_even_sh(self.dirs, np.max(l_max))

        s_l_clusters_all = np.zeros((int(np.max(l_max) / 2 + 1), n_voxels, n_shells))
        s_lm_clusters_all = np.zeros((int(np.max(n_sh_coeffs)), n_voxels, n_shells))
        id_measurements = np.arange(n_meas)
        ids_current_cluster_all = []
        for i in range(n_shells):
            ids_current_cluster = (
                id_measurements[(abs(bb[0,i] - self.b) < self.merge_distance) &
                (abs(bb[1,i] - self.beta) < self.merge_distance) &
                (abs(bb[3,i] - self.echo_time) < 1e-1)]
            )
            if abs(len(ids_current_cluster) - bb[2,i]) > 0.1:
                raise Exception('count of elements in current cluster failed, check shell table and B inputs')
    
            if (abs(bb[1,i]) > 0.045) & (bb[0,i] > 0.05):
                ids_current_lm = np.arange(n_sh_coeffs[i])
                y_lm_current_cluster = y_lm_matrix[ids_current_cluster,:]
                s_lm_current_cluster, _,_,_ = scl.lstsq(
                    y_lm_current_cluster[:, ids_current_lm], 
                    dwi_2D[ids_current_cluster,:]
                    )
                s_l_current_cluster = np.zeros((int(l_max[i]/2+1), n_voxels))
                s_l_current_cluster[0,:] = abs(s_lm_current_cluster[0,:])
                for j in range(2, int(l_max[i]) + 2, 2):
                    ids_currentl_m = ((1/2*(j+1)*(j+2)-j) + np.arange(-j,j + 1) - 1).astype(int)
                    s_l_current_cluster[int(j/2),:] = np.sqrt(
                        np.sum(
                        s_lm_current_cluster[ids_currentl_m,:]**2, axis=0)
                        )      
            else:
                s_lm_current_cluster = np.zeros((int(n_sh_coeffs[i]), n_voxels))
                s_l_current_cluster = np.zeros((int(l_max[i]/2+1), n_voxels))
                s_lm_current_cluster[0,:] = (np.sqrt(4*np.pi) * 
                    np.mean(dwi_2D[ids_current_cluster,:], axis=0)
                    )
                s_l_current_cluster[0,:] = abs(s_lm_current_cluster[0,:])
               
            s_lm_clusters_all[:n_sh_coeffs[i],:,i] = (
                s_lm_current_cluster / 
                np.tile(n_l_all[:n_sh_coeffs[i]], 
                (n_voxels, 1)).T
                )
            s_l_clusters_all[:int(l_max[i]/2+1),:,i] = (
                s_l_current_cluster / 
                np.tile(n_l_unique[:int(l_max[i]/2+1)], 
                (n_voxels, 1)).T
                )
            ids_current_cluster_all.append(
                np.vstack((ids_current_cluster * 0 + i, ids_current_cluster))
                )

        if rank1sl:
            sl_dn = s_l_clusters_all.copy()
            slm_dn = s_lm_clusters_all.copy()
            id_shells = np.arange(n_shells)
            for ll in range(2, int(np.max(l_max)) + 2, 2):
                if np.sum(l_max >= ll) > 1:
                    id_useful_shells = id_shells[l_max>=ll]
                    id_current_band = (l_all == ll)
                    slm_voxels_raw = s_lm_clusters_all[id_current_band, :, :][:, :, id_useful_shells]
                    slm_voxels_dn = self.low_rank_denoising(slm_voxels_raw.transpose(1,0,2), 1).transpose(1,0,2)
                    slm_dn[id_current_band, :, :][:, :, id_useful_shells]
                    sl_dn[int(ll/2),:,id_useful_shells] = np.sqrt(np.sum(slm_voxels_dn**2, axis=0)).T
                  
            s_lm_clusters_all = slm_dn.copy()
            s_l_clusters_all = sl_dn.copy()            

        if flag_4D:
            sl = np.zeros((sz_dwi[0], sz_dwi[1], sz_dwi[2], int(np.max(l_max)/2+1), n_shells))
            slm = np.zeros((sz_dwi[0], sz_dwi[1], sz_dwi[2], np.max(n_sh_coeffs), n_shells))
            for i in range(n_shells):
                try:
                    sl[:,:,:,:,i] = self.vectorize(s_l_clusters_all[:,:,i], self.mask)
                    slm[:,:,:,:,i] = self.vectorize(s_lm_clusters_all[:,:,i], self.mask)
                except:
                    sl[:,:,:,:,i] = self.vectorize(s_l_clusters_all[:,:,i], self.mask)[...,None]
                    slm[:,:,:,:,i] = self.vectorize(s_lm_clusters_all[:,:,i], self.mask)[...,None]
        else:
            sl = s_l_clusters_all
            slm = s_lm_clusters_all
        
        ids = np.hstack(ids_current_cluster_all)

        if self.project_forward_SH:

            group_clusters_all = ids[0, :]
            ids_clusters_all = ids[1, :]
            n_clusters = len(np.unique(group_clusters_all))
            signal = np.zeros_like(dwi_2D)
            for ii in range(n_clusters):
                ids_current_cluster = ids_clusters_all[group_clusters_all == ii]
                s_lm_current_cluster = s_lm_clusters_all[:,:,ii]
                if np.sqrt(np.sum(s_lm_current_cluster[1:,:]**2)) < 1e-3:
                    signal[ids_current_cluster,:] = n_l_all[0] * (y_lm_matrix[[ids_current_cluster],0].T @ s_lm_current_cluster[[0],:])
                else:
                    signal[ids_current_cluster,:] = y_lm_matrix[ids_current_cluster,:] @ (np.tile(n_l_all, (n_voxels,1)).T * s_lm_current_cluster)
            self.signal_4d = self.vectorize(signal, self.mask)

        num_elements = self.rotinv_lmax // 2 + 1
        slices = [np.squeeze(sl[:, :, :, i, :]) for i in range(num_elements)]
        sl = np.concatenate(slices, axis=3)
    
        return sl.reshape(sz_dwi[0],sz_dwi[1],sz_dwi[2],-1)
        
    def low_rank_denoising(self, X, p):
        u,s,v = np.linalg.svd(X, full_matrices=False)
        s_dn = np.zeros((s.shape[0], s.shape[1], s.shape[1]))
        diag_inds = np.diag_indices(X.shape[2])
        s_dn[:,diag_inds[0][:p],diag_inds[1][:p]] = s[:,:p]
        u_s = np.einsum('ijk,ikl->ijl', u, s_dn)
        return np.einsum('ijk,ikl->ijl', u_s, v)
        

    def group_dwi_in_shells_b_beta_te(self):
        """
        Function to group b-shells, echo times, and shell-linearity into clusters
        Inputs
        ------
        b-values:  (inherited from self)
        beta:      (inherited from self)
        echo time: (inherited from self)

        Output
        ------
        table_4d: sorted array of b-values, b-shapes, cluster sizes, and echo times
        """
        
        # create an array of all b-values, b-shapes, and echo times
        clusters = np.vstack((self.b, self.beta, self.echo_time)).T
        # cluster acp parameter data based on euclidean distance
        id_clusters = sch.fclusterdata(clusters, self.merge_distance, criterion='distance')
        
        sorted_ids_int = []
        shell_ids_int = []
        id_original_sortparam = np.arange(len(id_clusters))
        table_4d = np.zeros((4, len(np.unique(id_clusters))))

        # loop over cluster indices and create an array decribing each cluster
        for idx, i in enumerate(np.unique(id_clusters)):
            table_4d[:,idx] = np.vstack((
                np.mean(self.b[id_clusters==i]), 
                np.mean(self.beta[id_clusters==i]), 
                np.sum(id_clusters==i),
                np.mean(self.echo_time[id_clusters==i]))
                ).squeeze()
            sorted_ids_int.append(id_original_sortparam[id_clusters==i])
            shell_ids_int.append(i * np.ones(np.sum(id_clusters==i)))
        sorted_ids_int = np.hstack(sorted_ids_int)
        shell_ids_int = np.hstack(shell_ids_int)
        
        # ensure b-shell units are on the order of 1
        table_aux = table_4d.copy()
        table_aux[1,:] = np.round(table_aux[1,:], 2)
        if self.flag_microstructure_units:
            table_aux[0,:] = np.round(table_aux[0,:], 2)
        else:
            table_aux[0,:] = np.round(table_aux[0,:], -1)

        # sort shell table from smallest to largest shells
        sorted_inds = np.lexsort((table_aux[3, :], table_aux[1, :], table_aux[0, :]))
        table_aux_sorted = table_aux[:, sorted_inds]
        
        sorted_ids = []
        shell_ids = np.array(sorted_ids_int * 0)
        for i in range(table_aux.shape[1]):
            aux = sorted_ids_int[shell_ids_int == sorted_inds[i]]
            sorted_ids.append(aux)
            shell_ids[aux] = i

        return table_aux_sorted
        

    def get_default_lmax(self):
        """
        Function to automatically calculate the best L-max for spherical harmonic fitting
        """
        table_4d = self.group_dwi_in_shells_b_beta_te()
        self.table_4d = table_4d

        b = table_4d[0,:]
        beta = table_4d[1,:]
        n_dirs = table_4d[2,:]
        l_max = np.zeros(table_4d.shape[1])

        if max(b) < 100:
            b_micro_units = b
        else:
            b_micro_units = b / 1000

        for i in range(len(b_micro_units)):
            if abs(beta[i]) < 0.1 or b_micro_units[i] < 0.1:
                l_max[i] = 0
            else:
                if b_micro_units[i] < 1.2:
                    l_max[i] = 2
                elif b_micro_units[i] <= 2.5:
                    l_max[i] = 4
                elif b_micro_units[i] < 6.5:
                    l_max[i] = 6
                else:
                    l_max[i] = 8
                n_free_param = l_max[i] * (l_max[i] + 3) / 2 + 1
                while n_free_param >= n_dirs[i] / 1.3:
                    l_max[i] = l_max[i] - 2
                    n_free_param = l_max[i] * (l_max[i] + 3) / 2 + 1
        
        return l_max

    def compute_extended_moments(self, moments, degree):
        """
        This function computes the extended matrix of moments for polynomial
        interpolation in multiple dimensions.

        Moments is [NxD] where N is the number of samples and D the dimension of
        the problem (the number of different moments)

        e.g. in 1D this extended matrix equals the Vandermonde matrix
        Moments_Extended = [1 X X^2 ... X^degree]

        Inputs
        ------
        moments: N, M ND array
        degree : int

        Outputs
        -------
        extended moments: N, M*degree ND array
        """

        n = moments.shape[0]
        n_moments = moments.shape[1]

        # stack ones onto the list of extended moments
        x0 = np.ones(n)
        if degree > 1:
            x1 = moments.copy()
            x = np.hstack((x0[...,None], x1))
        else:
            x = x0

        # loop up to input degree, computing the high order moments for each degree
        if degree >= 2:
            for current_degree in range(2, degree + 1):
                combinations = np.array(list(
                    itertools.combinations(
                        np.arange(n_moments + current_degree - 1), current_degree)
                ))

                n_total_combs = combinations.shape[0]
                flags = np.zeros((n_total_combs, n_moments + current_degree - 1), dtype=bool)
                sums = np.zeros((n_total_combs, n_moments +current_degree - 1))
                which = np.zeros((n_total_combs, current_degree), dtype=int)
                for i in range(n_total_combs):
                    flags[i, combinations[i,:]] = True
                    sums[i,:] = np.cumsum(~flags[i,:])
                    which[i,:] = sums[i, combinations[i,:]]
                current_x = np.ones((n, n_total_combs))
                
                for i in range(n_total_combs):
                    for j in range(current_degree):
                        current_id = which[i,j]
                        current_x[:,i] = current_x[:,i] * moments[:,current_id]
                x = np.hstack((x, current_x))

        return x

    def rotinv_kell_wfw_b_beta_te_numerical(self, ell, shells, kernel):
        """
        Inputs
        ------
        kernel: [f,Da,Depar,Deperp,f_w,T2a,T2e]
            (ND array with size: [Nvoxels x 5] or [Nvoxels x 7])
        ell: l-th order Legendre polynomial (int)

        outputs
        -------
        k for order l (ND array)

        This function evaluates the SM kernel projected on the ell-th order
        Legendre polynomial. The integral is computed numerically
        using Gaussian quadrature with 200 points (more than enough accuracy).
        
        b is the b-value of the shell in um/ms2 x=[f,Da,Depar-Deperp,Deperp,f_extra];
        
        By: Santiago Coelho
        
        
        """
        n_voxels = kernel.shape[0]
        if not np.isscalar(ell):
            raise Exception(' Only scalar L values supported')
        b = shells[0,:].squeeze()
        beta = shells[1,:].squeeze()
        te = shells[-1,:].squeeze()

        f = kernel[:,[0]]
        da = kernel[:,[1]]
        depar = kernel[:,[2]]
        deperp = kernel[:,[3]]
        deltae = depar - deperp
        fw = kernel[:,[4]]
        f_extra = 1-f-fw

        if kernel.shape[1] == 7:
            t2a = kernel[:,[5]]
            t2e = kernel[:,[6]]
        else:
            t2a = 1000 * np.ones_like(f)
            t2e = 1000 * np.ones_like(f)

        d_fw = 3
        t2_fw = 1500
        n = 200

        
        dwd = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        mat = sio.loadmat(os.path.join(dwd,'constant','constant_values.mat'))

        xx = mat['xx']
        w = mat['w']

        x = np.tile(xx, (1, len(b)))
        bval = np.tile(b, (n, 1))
        bshape = np.tile(beta, (n, 1))
        te_rep = np.tile(te, (n, 1))

        if ell == 0:
            pell_xi = 0 * x + 1
        elif ell == 2:
            pell_xi = 3/2 * x**2 - 1/2
        elif ell == 4:
            pell_xi = (35 * x**4 - 30 * x**2 + 3) / 8
        elif ell == 6:
            pell_xi = (231 * x**6 - 315 * x**4 + 105 * x**2 - 5) / 16
        elif ell == 8:
            pell_xi = (6435 * x**8 - 12012 * x**6 + 6930 * x**4 - 1260 * x**2 + 35) / 128
        
        kell = np.zeros((n_voxels, len(b)))
        for i in range(n_voxels):
            kell[i,:] = w @ ((f[i] * np.exp(-bval * bshape * da[i] * x**2 - bval * da[i]/3 * (1-bshape) - te_rep/t2a[i]) + 
                f_extra[i] * np.exp(-bval * bshape * deltae[i] * x**2 - bval * deltae[i]/3 * (1-bshape) - bval * deperp[i] - te_rep/t2e[i]) +
                fw[i] * np.exp(-bval * d_fw - te_rep/t2_fw)) * pell_xi)
        
        abs_beta = abs(beta)
       
        if ell > 0:
            kell[:, b < 1e-6] = 0
            kell[:, abs_beta < 1e-6] = 0
        else:
            kell[:, b < 1e-6] = (f * np.exp(-te[b < 1e-6] / t2a) + 
                f_extra * np.exp(-te[b < 1e-6] / t2e) + 
                (1 - f - f_extra) * np.exp(-te[b < 1e-6] / t2_fw)
                )

            if (abs_beta < 1e-6).any():
                kell[:, abs_beta < 1e-6] = (f * np.exp(-te[abs_beta < 1e-6, None] / t2a) *
                    np.exp(-b[abs_beta < 1e-6, None] * da / 3) + 
                    f_extra * np.exp(-te[abs_beta < 1e-6, None] / t2e) * 
                    np.exp(-b[abs_beta < 1e-6, None] * deltae / 3 - b[abs_beta < 1e-6, None] * deperp) + 
                    fw * np.exp(-te[abs_beta < 1e-6, None] / t2_fw) * 
                    np.exp(-b[abs_beta < 1e-6, None] * d_fw)
                    )

        return kell

    def generate_sm_wfw_b_beta_te_ws0_training_data(self):
        """
        Function to generate k tensor training data from a prior
        Inputs:
        -------
        prior: (inherited from self)
        Outputs:
        --------
        kernel_params: NDarray
        """

        f = self.prior[:,0]

        n_training = len(f)
        if not self.flag_compartments[0]:
            f = np.zeros((n_training, 1))
        else:
            f = self.prior[:,[0]]

        if not self.flag_compartments[2]:
            f_fw = np.zeros((n_training, 1))
        else:
            f_fw = self.prior[:,[5]]

        if not self.flag_compartments[1]:
            f = 1 - f_fw

        s0 = np.ones((len(f), 1))
        f_extra = 1 - f - f_fw
        
        kernel_params = np.hstack((s0, f, self.prior[:,1:4], f_fw, self.prior[:, 5:]))
        shells = self.table_4d[[0,1,3],:]

        k0_all = self.rotinv_kell_wfw_b_beta_te_numerical(0, shells, kernel_params[:,1:8])
        if self.rotinv_lmax == 2:
            k2_all = self.rotinv_kell_wfw_b_beta_te_numerical(2, shells, kernel_params[:,1:8])
            rot_invs = s0 * np.hstack((k0_all, 
                kernel_params[:,[8]] * abs(k2_all)))
        elif self.rotinv_lmax == 4:
            k2_all = self.rotinv_kell_wfw_b_beta_te_numerical(2, shells, kernel_params[:,1:8])
            k4_all = self.rotinv_kell_wfw_b_beta_te_numerical(4, shells, kernel_params[:,1:8])
            rot_invs = s0 * np.hstack((k0_all, 
                kernel_params[:,[8]] * abs(k2_all), 
                kernel_params[:,[9]] * abs(k4_all)))
        elif self.rotinv_lmax == 6:
            k2_all = self.rotinv_kell_wfw_b_beta_te_numerical(2, shells, kernel_params[:,1:8])
            k4_all = self.rotinv_kell_wfw_b_beta_te_numerical(4, shells, kernel_params[:,1:8])
            k6_all = self.rotinv_kell_wfw_b_beta_te_numerical(6, shells, kernel_params[:,1:8])
            rot_invs = s0 * np.hstack((k0_all, 
                kernel_params[:,[8]] * abs(k2_all), 
                kernel_params[:,[9]] * abs(k4_all),
                kernel_params[:,[10]] * abs(k6_all)))

        return rot_invs

    def standard_model_mlfit_rot_invs(self, rot_invs, sigma_norm_limits):
        """
        Standard model polynomial regression training and fitting

        Inputs:
        -------
        rot_invs: rotational invariants of the input diffusion image (NDarray)
        sigma_norm_limits: maximum and minimum values sigma/s0 can take

        Outputs:
        -------
        standard model kernel: 4D NDarray
        """

        if rot_invs.ndim == 4:
            flag_4d = 1
        elif rot_invs.ndim == 2:
            flag_4d = 0
        else:
            raise Exception('RotInvs must be a 2D or 4D array')

        if flag_4d:
            rot_invs_normalized = self.vectorize(rot_invs, self.mask)
            sigma_normalized = self.vectorize(self.sigma, self.mask)
        else:
            rot_invs_normalized = rot_invs
            sigma_normalized = self.sigma
        
        s0_lowest_te = rot_invs_normalized[0,:].copy()
        
        np.divide(
            rot_invs_normalized, s0_lowest_te, out=rot_invs_normalized, where=s0_lowest_te != 0
            )
        rot_invs_normalized = rot_invs_normalized.T
        np.divide(
            sigma_normalized, s0_lowest_te, out=sigma_normalized, where=s0_lowest_te != 0
            )
        
        shells = self.table_4d[0,:]
        beta = self.table_4d[1,:]
        n_dirs = self.table_4d[2,:]
        
        keep_non_zero_s0 = np.ones(self.table_4d.shape[1], dtype=bool)
        if self.rotinv_lmax == 0:
            keep_rot_invs_kernel = keep_non_zero_s0.copy()
        elif self.rotinv_lmax == 2:
            keep_non_zero_s2 = (abs(beta) > 0.2) & (shells > 0.1) & (self.l_max >=2)
            keep_rot_invs_kernel = np.hstack((
                keep_non_zero_s0, keep_non_zero_s2)
                )
        elif self.rotinv_lmax == 4:
            keep_non_zero_s2 = (abs(beta) > 0.2) & (shells > 0.1) & (self.l_max >=2)
            keep_non_zero_s4 = (abs(beta) > 0.2) & (shells > 0.1) & (self.l_max >= 4)
            keep_rot_invs_kernel = np.hstack((
                keep_non_zero_s0, keep_non_zero_s2, keep_non_zero_s4)
                )
        elif self.rotinv_lmax == 6:
            keep_non_zero_s2 = (abs(beta) > 0.2) & (shells > 0.1) & (self.l_max >=2)
            keep_non_zero_s4 = (abs(beta) > 0.2) & (shells > 0.1) & (self.l_max >= 4)
            keep_non_zero_s6 = (abs(beta) > 0.2) & (shells > 0.1) & (self.l_max >= 6)
            keep_rot_invs_kernel = np.hstack((
                keep_non_zero_s0, keep_non_zero_s2, keep_non_zero_s4, keep_non_zero_s6)
                )

        sigma_noise_norm_levels_edges = np.linspace(
            sigma_norm_limits[0], sigma_norm_limits[1], self.n_levels + 1
            )
        sigma_noise_norm_levels_ids = np.digitize(
            sigma_normalized, sigma_noise_norm_levels_edges
            )
    
        sigma_noise_norm_levels_ids[sigma_normalized < sigma_noise_norm_levels_edges[0]] = 0
        sigma_noise_norm_levels_ids[sigma_normalized >= sigma_noise_norm_levels_edges[-1]] = self.n_levels
        sigma_noise_norm_levels_mean = 1/2 * (
           sigma_noise_norm_levels_edges[1:] + sigma_noise_norm_levels_edges[:-1]
           )

        degree_kernel = 3
        x_fit_norm = self.compute_extended_moments(
            rot_invs_normalized[:, keep_rot_invs_kernel], degree_kernel
            )
        
        n_voxels_masked = rot_invs_normalized.shape[0]
        f_ml_fit = np.zeros(n_voxels_masked)
        da_ml_fit = np.zeros(n_voxels_masked)
        depar_ml_fit = np.zeros(n_voxels_masked)
        deperp_ml_fit = np.zeros(n_voxels_masked)
        f_fw_ml_fit = np.zeros(n_voxels_masked)
        t2a_ml_fit = np.zeros(n_voxels_masked)
        t2e_ml_fit = np.zeros(n_voxels_masked)

        f = self.prior[:,0]
        da = self.prior[:,1]
        depar = self.prior[:,2]
        deperp = self.prior[:,3]
        f_fw = self.prior[:,4]
        t2a = self.prior[:,5]
        t2e = self.prior[:,6]
        n_training = len(f)
        if not self.flag_compartments[0]:
            f = np.zeros(n_training)
        if not self.flag_compartments[2]:
            f_fw = np.zeros(n_training)
        if not self.flag_compartments[1]:
            f = 1 - f_fw

        rotinvs_train = self.generate_sm_wfw_b_beta_te_ws0_training_data()

        if self.rotinv_lmax == 0:
            p2_ml_fit = np.zeros(n_voxels_masked)
            sigma_ndirs_factor = np.sqrt(n_dirs)
            rotinvs_train = rotinvs_train[:,:-1//2]
        if self.rotinv_lmax == 2:
            p2_ml_fit = np.zeros(n_voxels_masked)
            sigma_ndirs_factor = np.sqrt(np.hstack((n_dirs, n_dirs * 5)))
        elif self.rotinv_lmax == 4:
            p2_ml_fit = np.zeros(n_voxels_masked)
            p4_ml_fit = np.zeros(n_voxels_masked)
            sigma_ndirs_factor = np.sqrt(np.hstack((n_dirs, n_dirs * 5, n_dirs * 9)))
        elif self.rotinv_lmax == 6:
            p2_ml_fit = np.zeros(n_voxels_masked)
            p4_ml_fit = np.zeros(n_voxels_masked)
            p6_ml_fit = np.zeros(n_voxels_masked)
            sigma_ndirs_factor = np.sqrt(np.hstack((n_dirs, n_dirs * 5, n_dirs * 9, n_dirs * 13)))
        
        rotinvs_train_norm = rotinvs_train / rotinvs_train[:,[0]]

        for i in range(1, len(sigma_noise_norm_levels_edges)):
            flag_current_noise_level = (sigma_noise_norm_levels_ids == i)

            if not np.any(flag_current_noise_level):
                continue
            
            sigma_rotinvs_training = sigma_noise_norm_levels_mean[i-1] / sigma_ndirs_factor
            meas_rotinvs_train = (rotinvs_train_norm + 
                sigma_rotinvs_training * np.random.standard_normal(size=rotinvs_train_norm.shape)
                )
            
            x_train = self.compute_extended_moments(
                meas_rotinvs_train[:, keep_rot_invs_kernel], degree=degree_kernel)
        
            pinv_x = scl.pinv(x_train)
            # pinv_x = np.linalg.pinv(x_train)
            coeffs_f = pinv_x @ f
            coeffs_da = pinv_x @ da
            coeffs_depar = pinv_x @ depar
            coeffs_deperp = pinv_x @ deperp
            coeffs_f_fw = pinv_x @ f_fw
            coeffs_t2a = pinv_x @ t2a
            coeffs_t2e = pinv_x @ t2e
            
            f_ml_fit[flag_current_noise_level] = (
                x_fit_norm[flag_current_noise_level, :] @ coeffs_f
                )
            da_ml_fit[flag_current_noise_level] = (
                x_fit_norm[flag_current_noise_level, :] @ coeffs_da
                )            
            depar_ml_fit[flag_current_noise_level] = (
                x_fit_norm[flag_current_noise_level, :] @ coeffs_depar
                )
            deperp_ml_fit[flag_current_noise_level] = (
                x_fit_norm[flag_current_noise_level, :] @ coeffs_deperp
                )
            f_fw_ml_fit[flag_current_noise_level] = (
                x_fit_norm[flag_current_noise_level, :] @ coeffs_f_fw
                )
            t2a_ml_fit[flag_current_noise_level] = (
                x_fit_norm[flag_current_noise_level, :] @ coeffs_t2a
                )
            t2e_ml_fit[flag_current_noise_level] = (
                x_fit_norm[flag_current_noise_level, :] @ coeffs_t2e
                )

            if self.rotinv_lmax >= 2:
                p2 = self.prior[:, 7]
                coeffs_p2 = pinv_x @ p2
                p2_ml_fit[flag_current_noise_level] = (
                    x_fit_norm[flag_current_noise_level, :] @ coeffs_p2
                    )
                p2_ml_fit[p2_ml_fit < 0] = 0 
                p2_ml_fit[p2_ml_fit > 1] = 1
            if self.rotinv_lmax >= 4:
                p4 = self.prior[:, 8]
                coeffs_p4 = pinv_x @ p4
                p4_ml_fit[flag_current_noise_level] = (
                    x_fit_norm[flag_current_noise_level, :] @ coeffs_p4
                    )
                p4_ml_fit[p4_ml_fit < 0] = 0 
                p4_ml_fit[p4_ml_fit > 1] = 1
            if self.rotinv_lmax >=6:
                p6 = self.prior[:, 9]
                coeffs_p6 = pinv_x @ p6
                p6_ml_fit[flag_current_noise_level] = (
                    x_fit_norm[flag_current_noise_level, :] @ coeffs_p6
                    )
                p6_ml_fit[p6_ml_fit < 0] = 0 
                p6_ml_fit[p6_ml_fit > 1] = 1

        f_ml_fit[f_ml_fit < 0] = 0
        f_ml_fit[f_ml_fit > 1] = 1
        da_ml_fit[da_ml_fit < 0] = 0
        da_ml_fit[da_ml_fit > 3] = 3
        depar_ml_fit[depar_ml_fit < 0] = 0
        depar_ml_fit[depar_ml_fit > 3] = 3
        deperp_ml_fit[deperp_ml_fit < 0] = 0
        deperp_ml_fit[deperp_ml_fit > 3] = 3
        f_fw_ml_fit[f_fw_ml_fit < 0] = 0
        f_fw_ml_fit[f_fw_ml_fit > 1] = 1
        t2a_ml_fit[t2a_ml_fit < 30] = 30
        t2a_ml_fit[t2a_ml_fit > 150] = 150
        t2e_ml_fit[t2e_ml_fit < 30] = 30
        t2e_ml_fit[t2e_ml_fit > 150] = 150 
    
        if self.rotinv_lmax == 2:
            pn_ml_fit = p2_ml_fit
        elif self.rotinv_lmax == 4:
            pn_ml_fit = np.vstack((p2_ml_fit, p4_ml_fit))
        elif self.rotinv_lmax >= 6:
            pn_ml_fit = np.vstack((p2_ml_fit, p4_ml_fit, p6_ml_fit))
        
        if not self.fit_T2:
            kernel = np.vstack(
                (f_ml_fit, da_ml_fit, depar_ml_fit, deperp_ml_fit, f_fw_ml_fit, pn_ml_fit)
                )
        else:
            kernel = np.vstack(
                (f_ml_fit, da_ml_fit, depar_ml_fit, deperp_ml_fit, f_fw_ml_fit, t2a_ml_fit, t2e_ml_fit, pn_ml_fit)
                )

        if flag_4d:
            kernel = self.vectorize(kernel, self.mask)

        return kernel

    def get_plm_from_s_and_kernel(self, dwi, kernel):
        """
        Function to get plms from SMI kernel and dwi
        Inputs:
        -------

        Outputs:
        --------
        """

        lmax = int(np.max(self.l_max))
        l_all = np.repeat(np.arange(0, lmax + 2, 2), 2 * np.arange(0, lmax + 2, 2) + 1)
        n_l = np.sqrt((2 * l_all + 1) * 4 * np.pi)
        y_lm_matrix = self.get_even_sh(self.dirs, lmax)
        if lmax == 0:
            plm = np.nan
            warnings.warn('plm will not be estimated since only S0 is estimated')
        else:
            dwi_norm = self.vectorize(dwi, self.mask)
            kernel = self.vectorize(kernel, self.mask)

            f = kernel[0,:]
            da = kernel[1,:]
            depar = kernel[2,:]
            deperp = kernel[3,:]
            fw = kernel[4,:]

            n_dwi = len(self.b)
            n_voxels = len(f)

            if len(np.unique(self.echo_time)) > 1:
                t2a = kernel[5,:]
                t2e = kernel[6,:]
            else:
                t2a = np.ones_like(f)
                t2e = np.ones_like(f)

            x = np.vstack((f, da, depar, deperp, fw, t2a, t2e)).T
            shells = np.vstack((self.b, self.beta, self.echo_time))

            kell = np.zeros((n_voxels, n_dwi, int(lmax / 2 + 1)))
            kell[:,:,0] = self.rotinv_kell_wfw_b_beta_te_numerical(0, shells, x)
            if lmax >= 2:
                kell[:,:,1] = self.rotinv_kell_wfw_b_beta_te_numerical(2, shells, x)
            if lmax >= 4:
                kell[:,:,2] = self.rotinv_kell_wfw_b_beta_te_numerical(4, shells, x)
            if lmax >= 6:
                kell[:,:,3] = self.rotinv_kell_wfw_b_beta_te_numerical(6, shells, x)
            if lmax >= 8:
                kell[:,:,4] = self.rotinv_kell_wfw_b_beta_te_numerical(8, shells, x)
            kell = kell.transpose(2, 1, 0)

            plm = np.zeros((int(lmax * (lmax + 3) / 2 + 1), n_voxels))
            nlm = (2 * np.arange(0, lmax + 2, 2) + 1).astype(int)
            for i in range(n_voxels):
                kl_times_ylm = np.repeat(kell[:,:,i], nlm, axis=0).T * (y_lm_matrix * n_l)
                plm[:,i] = scl.lstsq(kl_times_ylm, dwi_norm[:,i])[0]

            plm = self.vectorize(plm, self.mask)
            norm_ell = 4 * np.pi * (2 * np.arange(2, 10, 2) + 1)

            if lmax == 2:
                pl = np.sqrt(np.sum(plm[:,:,:,:5]**2, axis=3) / norm_ell[0])
            elif lmax == 4:
                pl = np.stack((np.sqrt(np.sum(plm[:,:,:,:5]**2, axis=3) / norm_ell[0]),
                            np.sqrt(np.sum(plm[:,:,:,5:14]**2, axis=3) / norm_ell[1])
                ), axis=3)
            elif lmax == 6:
                pl = np.stack((np.sqrt(np.sum(plm[:,:,:,:5]**2, axis=3) / norm_ell[0]),
                            np.sqrt(np.sum(plm[:,:,:,5:14]**2, axis=3) / norm_ell[1]),
                            np.sqrt(np.sum(plm[:,:,:,14:27]**2, axis=3) / norm_ell[2])
                ), axis=3)
            elif lmax == 8:
                pl = np.stack((np.sqrt(np.sum(plm[:,:,:,:5]**2, axis=3) / norm_ell[0]),
                            np.sqrt(np.sum(plm[:,:,:,5:14]**2, axis=3) / norm_ell[1]),
                            np.sqrt(np.sum(plm[:,:,:,14:27]**2, axis=3) / norm_ell[2]),
                            np.sqrt(np.sum(plm[:,:,:,27:44]**2, axis=3) / norm_ell[3])
                ), axis=3)

        return plm, pl

    def compute_eps_odf_rectification(self, plm):
        """
        Function to compute epsilon and rectify a previously computed set of plms

        Inputs:
        -------

        Outputs:
        -------
        """

        n_odfs = plm.shape[1]
        dwd = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        mat = sio.loadmat(os.path.join(dwd,'constant','constant_values.mat'))
        dirs = mat['dirs']
        leb_weights = mat['leb_weights']

        n_lm = plm.shape[0]
        lmax = np.sqrt(2 * n_lm + 9/4) - 3/2
        y_lm_n = self.get_even_sh(dirs, lmax)
        epsilon_all = np.zeros(n_odfs)
        plm_rect_all = np.zeros((n_lm, n_odfs))
        
        for i in range(n_odfs):
            objfun = (lambda epsilon, leb_weights, y_lm_n, plm: 
                abs(leb_weights @ (
                abs(y_lm_n @ np.hstack((1/np.sqrt(4*np.pi), plm)) - epsilon) - 
                y_lm_n @ np.hstack((1/np.sqrt(4*np.pi), plm)) - epsilon)
                ))
            epsilon_all[i] = sco.fmin(func=objfun, x0=0.1, ftol=1e-5, maxiter=1000, args=(leb_weights, y_lm_n, plm[:,i]))

            rect_odf = 1/2 * (abs(y_lm_n @ np.hstack((1/np.sqrt(4*np.pi), plm[:,i])) - epsilon_all[i]) + 
                y_lm_n @ np.hstack((1/np.sqrt(4*np.pi), plm[:,i])) - epsilon_all[i])

            plm_rect = leb_weights @ (rect_odf[...,None] * y_lm_n)
            plm_rect = plm_rect[:,1:]
            plm_rect_all[:,i] = plm_rect
        
        plm_rect = plm_rect_all
        if lmax == 2:
            pl_rect = np.sqrt(np.sum(plm_rect[:5,:]**2, axis=3)).T
        elif lmax == 4:
            pl_rect = np.stack((np.sqrt(np.sum(plm_rect[:5,:]**2, axis=3)).T,
                        np.sqrt(np.sum(plm_rect[5:14,:]**2, axis=3)).T
            ), axis=0)
        elif lmax == 6:
            pl_rect = np.stack((np.sqrt(np.sum(plm_rect[:5,:]**2, axis=3)).T,
                        np.sqrt(np.sum(plm_rect[5:14,:]**2, axis=3)).T,
                        np.sqrt(np.sum(plm_rect[14:27,:]**2, axis=3)).T
            ), axis=0)
        elif lmax == 8:
            pl_rect = np.stack((np.sqrt(np.sum(plm_rect[:5,:]**2, axis=3)).T,
                        np.sqrt(np.sum(plm_rect[5:14,:]**2, axis=3)).T,
                        np.sqrt(np.sum(plm_rect[14:27,:]**2, axis=3)).T,
                        np.sqrt(np.sum(plm_rect[27:44,:]**2, axis=3)).T
            ), axis=0)

        return epsilon_all.T

    def fit(self, dwi, mask=None, sigma=None):
        """
        Main fitting function for SMI class
        Inputs: 
        -------
        dwi (4d NDarray)
        bval (1d NDarray)
        bvec (2d NDarray)
        optional (mask: 3d NDarray, sigma: 3d NDarray)

        Outputs:
        --------
        kernel (4d NDarray)
        RotInvs (4d NDarray)

        """

        self.dwishape = dwi.shape    
        self.set_mask(mask)

        self.set_sigma(sigma, dwi)

        # set l_max for smi and rotinv
        self.set_lmax_rotinv_smi()

        # set SMI priors
        self.set_priors()

        dwi = dwi.astype(np.float32)
        # correct for rician bias if the flag is on
        if self.flag_rician_bias:
            dwi = np.sqrt(abs(dwi**2 - sigma[...,None]**2))

        # spherical harmonic fit
        rot_invs = self.fit2D4D_LLS_RealSphHarm_wSorting_norm_var(dwi, self.mask)
        
        output = {}
        # standard model fit
        kernel = self.standard_model_mlfit_rot_invs(rot_invs, [0, 0.2])
        output['f'] = kernel[:,:,:,0]
        output['Da'] = kernel[:,:,:,1]
        output['DePar'] = kernel[:,:,:,2]
        output['DePerp'] = kernel[:,:,:,3]
        output['fw'] = kernel[:,:,:,4]
        if self.fit_T2:
            output['T2a'] = kernel[:,:,:,5]
            output['T2e'] = kernel[:,:,:,6]
            output['p2'] = kernel[:,:,:,7]
            if self.rotinv_lmax == 4:
                output['p4'] = kernel[:,:,:,8]
            if self.rotinv_lmax == 6:
                output['p6'] = kernel[:,:,:,9]
        else:
            output['p2'] = kernel[:,:,:,5]
            if self.rotinv_lmax == 4:
                output['p4'] = kernel[:,:,:,6]
            if self.rotinv_lmax == 6:
                output['p6'] = kernel[:,:,:,7]
        output['rotinvs'] = rot_invs


        if self.flag_fit_fodf:
            s0 = rot_invs[:,:,:,0]
            dwi_norm = np.divide(dwi, s0, where=s0 != 0)
            (plm, pl) = self.get_plm_from_s_and_kernel(dwi_norm, kernel)
            output['plm'] = plm
            output['pl'] = pl

        if self.flag_rectify_fodf and self.flag_fit_fodf:
            plm = self.vectorize(plm, self.mask)
            epsilon = self.compute_eps_odf_rectification(plm[1:, :])
            epsilon = self.vectorize(epsilon, mask)
            output['epsilon'] = epsilon

        return output
        
if __name__ == "__main__":
    sys.exit()
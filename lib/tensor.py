from joblib import Parallel, delayed
import numpy as np
from tqdm import tqdm
from lib.mpunits import vectorize
import os

import cvxpy as cp
import scipy.linalg as scl
import scipy.io as sio
from scipy.ndimage import zoom, gaussian_filter

import multiprocessing
import warnings
warnings.filterwarnings("ignore")

class TensorFitting(object):

    def __init__(self, grad, n_cores):
        self.grad = grad
        self.n_cores = n_cores

    def create_tensor_order(self, order):
        """
        Function to create a set of indices describing the order of a 6 parameter or 21 parameter diffusion tensor
        for a diffusion tensor, the default order goes as xx, xy, xz, yy, yz, zz, and so on
        """
        if order == 2:
            cnt = np.array([1, 2, 2, 1, 2, 1], dtype=int)
            ind = np.array(([1, 1], [1, 2], [1, 3], [2, 2], [2, 3], [3, 3])) - 1
        if order == 4:
            cnt = np.array([1, 4, 4, 6, 12, 6, 4, 12, 12, 4, 1, 4, 6, 4, 1], dtype=int)
            ind = np.array(([1,1,1,1],[1,1,1,2],[1,1,1,3],[1,1,2,2],[1,1,2,3],[1,1,3,3],\
                [1,2,2,2],[1,2,2,3],[1,2,3,3],[1,3,3,3],[2,2,2,2],[2,2,2,3],[2,2,3,3],[2,3,3,3],[3,3,3,3])) - 1
        return cnt, ind

    def diffusion_coeff(self, dt, dir):
        """
        Compute the apparent diffusion coefficient
        """
        # compute ADC
        dcnt, dind = self.create_tensor_order(2)
        ndir = dir.shape[0]
        bD = np.tile(dcnt,(ndir, 1)) * dir[:,dind[:, 0]] * dir[:,dind[:, 1]]
        adc = bD @ dt
        return adc


    def kurtosis_coeff(self, dt, dir):
        """
        Function to Compute the apparent kurtosis coefficient
        """
        # compute AKC
        wcnt, wind = self.create_tensor_order(4)
        ndir = dir.shape[0]
        T = np.tile(wcnt,(ndir, 1)) * dir[:,wind[:, 0]] * dir[:,wind[:, 1]] * dir[:,wind[:, 2]] * dir[:,wind[:, 3]]
        akc = T @ dt[6:]

        adc = self.diffusion_coeff(dt[:6], dir)
        md = np.sum(dt[np.array([0,3,5])], 0)/3
        akc = (akc * md[np.newaxis]**2) / (adc**2)
        return akc

    def w_kurtosis_coeff(self, dt, dir):
        """
        Function to compute W form of kurtosis coefficient
        """
        # compute AKC
        wcnt, wind = self.create_tensor_order(4)
        T = np.prod(np.reshape(dir[:,wind],(-1,15,4)), axis=2) @ np.diag(wcnt)
        akc = T @ dt[6:]
        return akc

    def fibonacci_sphere(self, samples=1, randomize=True):
        """
        Function to generate N "samples" evenly spaced points on a sphere
        """
        import random
        rnd = 1
        if randomize:
            rnd = random.random() * samples
        points = []
        offset = 2/samples
        increment = np.pi * (3. - np.sqrt(5.))
        for i in range(samples):
            y = ((i * offset) - 1) + (offset / 2)
            r = np.sqrt(1 - pow(y,2))
            phi = ((i + rnd) % samples) * increment
            x = np.cos(phi) * r
            z = np.sin(phi) * r
            points.append([x,y,z])
        return points

    def radial_sampling(self, dir, n):
        """
        Function to output the radial component of a set of evenly distribted direction vectors
        """
        # get the radial component of a set of directions
        dt = 2*np.pi/n
        theta = np.arange(0,2*np.pi,dt)
        dirs = np.vstack((np.cos(theta), np.sin(theta), 0*theta))
        
        v = np.hstack((-dir[1], dir[0], 0))
        s = np.sqrt(np.sum(v**2))
        c = dir[2]
        V = np.array([[0, -v[2], v[1]],[v[2], 0, -v[0]],[-v[1], v[0], 0]])
        R = np.eye(3) + V + np.matmul(V,V) * (1-c)/(s**2)
        dirs = np.matmul(R, dirs)
        return dirs 

    def wlls(self, shat, dwi, b, C=None):
        """
        Function to perform pseudoinversion on raw diffusion data
        """
        # compute a wlls fit using weights from inital fit shat
        w = np.diag(shat)
        try:
            if C is None:
                dt = scl.pinv(w @ b, check_finite=True) @ (w @ np.log(dwi))
            else:
                x = cp.Variable((22,1))
                A = w @ b
                B = w @ np.log(dwi[...,None])
                objective = cp.Minimize(cp.sum_squares(A @ x - B))
                constraints = [C @ x >= 0]
                prob = cp.Problem(objective, constraints)
                result = prob.solve()
                dt = x.value.squeeze()
                # dt,_ = opt.nnls(w @ b, w @ np.log(dwi))
                # res = opt.lsq_linear(w @ b, w @ np.log(dwi), (-5, 1), method='trf', tol=1e-12, max_iter=220000)
                # dt = res.x
        except:
            dt = np.zeros((b.shape[1]))

        return dt  

    def dti_tensor_params(self, nn):
        """
        Function to compute orthogonal projections of diffusion from a spherical tensor representation
        """
        # compute dti tensor eigenvalues and eigenvectors and sort them
        values, vectors = np.linalg.eig(nn)
        idx = np.argsort(-values)
        values = -np.sort(-values)
        vectors = vectors[:, idx]
        return values, vectors  

    def dki_tensor_params(self, v1, dt, fit_w=False):
        """
        Function to compute orthogonal projections of kurtosis from a spherical tensor representation
        """
        # kurtosis tensor parameters use average directional
        # statistics to approximate ak and rk
        dirs = np.vstack((v1, -v1))
        if fit_w:
            akc = self.w_kurtosis_coeff(dt, dirs)
            ak = np.mean(akc)
            dirs = self.radial_sampling(v1, 256).T
            akc = self.w_kurtosis_coeff(dt, dirs)
            rk = np.mean(akc)
        else:
            akc = self.kurtosis_coeff(dt, dirs)
            ak = np.mean(akc)
            dirs = self.radial_sampling(v1, 256).T
            akc = self.kurtosis_coeff(dt, dirs)
            rk = np.mean(akc)
        return ak, rk    

    def extract_parameters(self, dt, b, mask, extract_dti, extract_dki, fit_w=False):
        """
        Computes parameters MD, AD, RD, FA and Trace of the diffusion tensor, along with color encoded FA
        Computes parameters MK AK and RK of the kurtosis tensor
        """
        # extract all tensor parameters from dt
        # num_cores = multiprocessing.cpu_count()

        DT = np.reshape(np.concatenate(
            (dt[0,:], dt[1,:], dt[2,:], dt[1,:], dt[3,:], dt[4,:], dt[2,:], dt[4,:], dt[5,:])
            ), (3, 3, dt.shape[1]))
        
        # get the trace
        rdwi = np.exp(np.matmul(b[:,1:], dt))
        B = np.round(self.grad[:,-1])
        uB = np.unique(B)
        trace = np.zeros((dt.shape[1], uB.shape[0]))
        for ib in range(0, uB.shape[0]): 
            t = np.where(B == uB[ib])[0]
            trace[:,ib] = np.mean(rdwi[t,:], axis=0)

        nvox = dt.shape[1]
        inputs = range(0, nvox)
        values, vectors = zip(*Parallel(n_jobs=self.n_cores, prefer='processes')\
            (delayed(self.dti_tensor_params)(DT[:,:,i]) for i in inputs))        
        values = np.reshape(np.abs(values), (nvox, 3))
        vectors = np.reshape(vectors, (nvox, 3, 3))

        l1 = vectorize(values[:,0], mask)
        l2 = vectorize(values[:,1], mask)
        l3 = vectorize(values[:,2], mask)
        v1 = vectorize(vectors[:,:,0].T, mask)

        md = (l1+l2+l3)/3
        rd = (l2+l3)/2
        ad = l1
        fa = np.sqrt(1/2)*np.sqrt((l1-l2)**2+(l2-l3)**2+(l3-l1)**2)/np.sqrt(l1**2+l2**2+l3**2)
        trace = vectorize(trace.T, mask)
        fe = np.abs(np.stack((fa*v1[:,:,:,0], fa*v1[:,:,:,1], fa*v1[:,:,:,2]), axis=3))
        
        parameters = {}
        if extract_dti:
            parameters['md'] = md
            parameters['rd'] = rd
            parameters['ad'] = ad
            parameters['fa'] = fa
            parameters['fe'] = fe
            parameters['trace'] = trace

        if extract_dki:
            #dirs = np.array(self.fibonacci_sphere(256, True))
            dwd = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            mat = sio.loadmat(os.path.join(dwd,'constant','dirs256.mat'))
            dirs = mat['dirs']
            
            
            if fit_w:
                akc = self.w_kurtosis_coeff(dt, dirs)
                mk = np.mean(akc, axis=0)
                ak, rk = zip(*Parallel(n_jobs=self.n_cores, prefer='processes')\
                    (delayed(self.dki_tensor_params)(vectors[i,:,0], dt[:,i], fit_w=True) for i in inputs))
                ak = np.reshape(ak, (nvox))
                rk = np.reshape(rk, (nvox))
                ak = vectorize(ak, mask) * md**2 / ad**2
                rk = vectorize(rk, mask) * md**2 / rd**2
            else:
                akc = self.kurtosis_coeff(dt, dirs)
                mk = np.mean(akc, axis=0)
                ak, rk = zip(*Parallel(n_jobs=self.n_cores, prefer='processes')\
                   (delayed(self.dki_tensor_params)(vectors[i,:,0], dt[:,i], fit_w=False) for i in inputs))
                ak = np.reshape(ak, (nvox))
                rk = np.reshape(rk, (nvox))
                ak = vectorize(ak, mask)
                rk = vectorize(rk, mask)
            mk = vectorize(mk, mask)

            parameters['mk'] = mk
            parameters['ak'] = ak
            parameters['rk'] = rk
        
        return parameters

    def dki_fit(self, dwi, mask, constraints=None):
        """ 
        Diffusion parameter estimation using weighted linear least squares fitting
        Outputs the 6 parameter diffusion tensor and 21 parameter kurtosis tensor in dt[diffusion,kurtosis]
        Outputs S0 the true quantitative mean signal at zero gradient strength
        the gradient tensor b
        """
        # run the fit
        order = np.floor(np.log(np.abs(np.max(self.grad[:,-1])+1))/np.log(10))
        if order >= 2:
            self.grad[:, -1] = self.grad[:, -1]/1000

        dwi.astype(np.double)
        dwi[dwi <= 0] = np.finfo(np.double).eps
        dwi[~np.isfinite(dwi)] = np.finfo(np.double).eps

        self.grad = self.grad.astype(np.double)
        normgrad = np.sqrt(np.sum(self.grad[:,:3]**2, 1))
        normgrad[normgrad == 0] = 1

        self.grad[:,:3] = self.grad[:,:3]/np.tile(normgrad, (3,1)).T
        self.grad[np.isnan(self.grad)] = 0

        dcnt, dind = self.create_tensor_order(2)
        wcnt, wind = self.create_tensor_order(4)

        ndwis = dwi.shape[-1]
        bs = np.ones((ndwis, 1))
        bD = np.tile(dcnt,(ndwis, 1))*self.grad[:,dind[:, 0]]*self.grad[:,dind[:, 1]]
        bW = np.tile(wcnt,(ndwis, 1))*self.grad[:,wind[:, 0]]*self.grad[:,wind[:, 1]]*self.grad[:,wind[:, 2]]*self.grad[:,wind[:, 3]]
        b = np.concatenate((bs, (np.tile(-self.grad[:,-1], (6,1)).T * bD), np.squeeze(1/6 * np.tile(self.grad[:,-1], (15,1)).T**2) * bW), 1)

        dwi_ = vectorize(dwi, mask)
        init = np.matmul(np.linalg.pinv(b), np.log(dwi_))
        shat = np.exp(np.matmul(b, init))

        if np.any(constraints):
            dwd = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            mat = sio.loadmat(os.path.join(dwd,'constant','dirs_constrained.mat'))
            dir = mat['dir']

            ndir = dir.shape[0]
            C = []
            if constraints[0] > 0:
                c0 = np.hstack((np.zeros((ndir,1)), 
                    np.tile(dcnt, (ndir, 1)) * dir[:, dind[:,0]] * dir[:,dind[:,1]], 
                    np.zeros((ndir, 15))))
                C.append(c0)
            if constraints[1] > 0:
                c1 = np.hstack((np.zeros((ndir,7)), 
                    np.tile(wcnt, (ndir, 1)) * dir[:, wind[:,0]] * dir[:,wind[:,1]] * dir[:,wind[:,2]] * dir[:,wind[:,3]]))
                C.append(c1)
            if constraints[2] > 0:
                c2 = np.hstack((np.zeros((ndir,1)), 
                    3/np.max(self.grad[:,3]) * np.tile(dcnt, (ndir, 1)) * dir[:, dind[:,0]] * dir[:,dind[:,1]], 
                    -np.tile(wcnt, (ndir, 1)) * dir[:,wind[:,0]] * dir[:,wind[:,1]] * dir[:,wind[:,2]] * dir[:,wind[:,3]]))
                C.append(c2)
            C = np.reshape(C, (ndir*np.sum(constraints), -1))
        else:
            C = None

        inputs = tqdm(range(0, dwi_.shape[1]))
        # num_cores = multiprocessing.cpu_count()

        # for i in inputs:
        #     self.wlls(shat[:,i], dwi_[:,i], b, C)
        
        V = 0
        if V == 1:
            # vectorized fit is slow
            import timeit
            starting_time = timeit.default_timer()
            shat_vw = shat.T[:, :, np.newaxis] * np.eye(b.shape[0])
            wb = np.einsum('ijk,kl->ijl',shat_vw, b)
            pinvwb = np.linalg.pinv(wb)
            wlogdwi = np.einsum('ijk,ji->ij',shat_vw, np.log(dwi_))
            dt = np.einsum('ijk,ik->ji',pinvwb, wlogdwi)
            print("Time difference :", timeit.default_timer() - starting_time)
        else:
            dt = Parallel(n_jobs=self.n_cores,prefer='processes')\
                (delayed(self.wlls)(shat[:,i], dwi_[:,i], b, C) for i in inputs)
            dt = np.reshape(dt, (dwi_.shape[1], b.shape[1])).T
        

        s0 = np.exp(dt[0,:])
        dt = dt[1:,:]
        D_apprSq = 1/(np.sum(dt[(0,3,5),:], axis=0)/3)**2
        dt[6:,:] = dt[6:,:]*np.tile(D_apprSq, (15,1))
        return dt, s0, b

    def dti_fit(self, dwi, mask):

        # run the fit
        order = np.floor(np.log(np.abs(np.max(self.grad[:,-1])+1))/np.log(10))
        if order >= 2:
            self.grad[:, -1] = self.grad[:, -1]/1000

        dwi.astype(np.double)
        dwi[dwi <= 0] = np.finfo(np.double).eps
        dwi[~np.isfinite(dwi)] = np.finfo(np.double).eps

        self.grad = self.grad.astype(np.double)
        normgrad = np.sqrt(np.sum(self.grad[:,:3]**2, 1))
        normgrad[normgrad == 0] = 1

        self.grad[:,:3] = self.grad[:,:3]/np.tile(normgrad, (3,1)).T
        self.grad[np.isnan(self.grad)] = 0

        dcnt, dind = self.create_tensor_order(2)

        ndwis = dwi.shape[-1]
        bs = np.ones((ndwis, 1))
        bD = np.tile(dcnt,(ndwis, 1))*self.grad[:,dind[:, 0]]*self.grad[:,dind[:, 1]]
        b = np.concatenate((bs, (np.tile(-self.grad[:,-1], (6,1)).T * bD)),1)

        dwi_ = vectorize(dwi, mask)
        init = np.matmul(np.linalg.pinv(b), np.log(dwi_))
        shat = np.exp(np.matmul(b, init))

        inputs = tqdm(range(0, dwi_.shape[1]))
        # num_cores = multiprocessing.cpu_count()
        
        dt = Parallel(n_jobs=self.n_cores, prefer='processes')\
            (delayed(self.wlls)(shat[:,i], dwi_[:,i], b) for i in inputs)
        
        dt = np.reshape(dt, (dwi_.shape[1], b.shape[1])).T

        s0 = np.exp(dt[0,:])
        dt = dt[1:,:]
        return dt, s0, b

    def compute_outliers(self, dt, dir):
        akc = self.kurtosis_coeff(dt, dir)
        return np.any((akc < -1) | (akc > 10), axis=0)

    def outlierdetection(self, dt, mask, dir):

        nvxls = dt.shape[1]
        akc_mask = np.zeros((nvxls))
        nblocks = 200

        try:
            N = 10000
            akc_mask = Parallel(n_jobs=self.n_cores, prefer='threads')\
                (delayed(self.compute_outliers)
                (dt, dir[int(N/nblocks*(i-1)):int(N/nblocks*i),:]) for i in range(1, nblocks + 1)
                )
        except:
            N = 1000
            akc_mask = Parallel(n_jobs=8)\
                (delayed(self.compute_outliers)
                (dt, dir[int(N/nblocks*(i-1)):int(N/nblocks*i),:]) for i in range(1, nblocks + 1)
                )

        return np.sum(akc_mask, axis=0)
    
    
    def identify_outliers(self, dt, percentiles):
        low_thresh = np.percentile(dt, percentiles[0], method='median_unbiased')
        high_thresh = np.percentile(dt, percentiles[1], method='median_unbiased')
        return (dt < low_thresh) | (dt > high_thresh)
    
    def identify_outliers_med(self, dt, threshold):
        d = np.abs(dt - np.median(dt))
        mdev = np.median(d)
        s = d/mdev if mdev else np.zeros(len(d))
        return s > threshold
    
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
        import itertools

        n = moments.shape[0]
        n_moments = moments.shape[1]

        # stack ones onto the list of extended moments
        x0 = np.ones(n)
        if degree >= 1:
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
    
    def train_bayes_fit(self, dwi, dt, s0, b, mask, flag_dti='False'):
        """
        This function trains and evaluates a polynomial regression model to update
        a wlls fit without outlier voxels. 

        Inputs
        ------
        dwi: ND array 
        dt : wlls evaluated diffusion tensor Nx21 tensor
        s0 : first order element of dt
        b  : gradient matrix

        Outputs
        -------
        extended moments: dt fit from regression
        """

        SNR = 100
        sigma = 1 / SNR

        maxb = 3
        # grad_orig = self.grad
        # grad_keep = self.grad[:,3] < maxb
        # dwi = dwi[..., grad_keep]
        
        dwi[dwi<=0] = np.finfo(dwi.dtype).eps
        D_apprSq = 1/(np.sum(dt[(0,3,5),:], axis=0)/3)**2
        if not flag_dti == 'True':
            dt[6:,:] = dt[6:,:] / np.tile(D_apprSq, (15,1))
        dt = np.vstack((np.log(s0), dt))

        # first identify and remove outliers in dt
        if flag_dti == 'True':
            outlier_range = (1, 99)
            trace = [1,4,6]
        else:
            outlier_range = (5, 95)
            trace = [1,4,6,7,10,12,17,19,21]
          

        outlier_inds = np.zeros((len(trace), dt.shape[1]))
        for i in range(len(trace)):
            outlier_inds[i,:] = dt[trace[i], :] < -0.2
        non_negative_trace = (1 - outlier_inds).all(axis=0)
        dt_ = dt[:, non_negative_trace]

        outlier_inds = np.zeros_like(dt_)
        for i in range(dt.shape[0]):
            outlier_inds[i,:] = self.identify_outliers(dt_[i,:], (3, 99))
        non_outlier_voxels = (1 - outlier_inds).all(axis=0)

        # interpolate non outliers in dt for more training data
        train_size = 400000
        dt_init = dt_[:, non_outlier_voxels]
        scale_fact = train_size / dt_init.shape[1]

        #scale_fact = int(np.ceil(.5e6 / dt_init.shape[1]))
        dt_interp_init = zoom(dt_init, (1, scale_fact))
           
        # # set up training (zeroing diagnoal W elements)
        # train_size = 400000
        # id_non_outliers = np.random.permutation(train_size)
        # dt_train = np.zeros((dt.shape[0], train_size))
        # #inds_exclude = np.array([9,10,12,14,15,16,17,19,21]) - 1
        # for id in range(dt.shape[0]):
        #     #if np.any(id == inds_exclude):
        #     #    dt_train[id,:] = 0
        #     #else:
        #     dt_train[id,:] = np.real((dt_interp_init[id, id_non_outliers]) * 1.0)


        dt_train = dt_interp_init.copy()
        dt_train[0,:] = -1 * np.random.exponential(0.707, (1,train_size))
        #dt_train[0,:] = np.random.uniform(0,1,(1, train_size)) * 2
        s_training = np.exp(b @ dt_train).T
        s_training += np.abs(sigma *
                        (np.random.standard_normal(s_training.shape) + 
                        1j*np.random.standard_normal(s_training.shape)))
        

        # compress the training data for speed
        npars = 15
        x = np.log(s_training)
        #u,s,v = np.linalg.svd((x - x.mean()) / x.std(), full_matrices=False)
        #reduced_dim = v[:,npars:]
        data_training = x #@ reduced_dim

        # polynomial regression
        order_poly_fit = 1
        x_moments = self.compute_extended_moments(data_training, order_poly_fit)
        pinv_x = np.linalg.pinv(x_moments)
        
        coeffs_dt = np.zeros((x_moments.shape[1], dt.shape[0]))
        for i in range(dt.shape[0]):
            coeffs_dt[:,i] = pinv_x @ dt_train[i,:].T

        # apply the trained coefficients to the original data
        data_testing = np.log(vectorize(dwi, mask) / s0).T #@ reduced_dim
        x_moments_test = self.compute_extended_moments(data_testing, order_poly_fit)
        dt_poly = np.zeros_like(dt)
        #dt_poly_train = np.zeros_like(dt_train)
        for i in range(dt.shape[0]):
            kn = x_moments_test @ coeffs_dt[:,i]
            kn[kn < np.min(dt_train[i,:])] = np.min(dt_train[i,:])
            kn[kn > np.max(dt_train[i,:])] = np.max(dt_train[i,:])
            dt_poly[i,:] = kn
            #dt_poly_train[i,:] = x_moments @ coeffs_dt[:,i]

        # import matplotlib.pyplot as plt
        # plt.plot(dt_train[1,:],dt_poly_train[1,:],'o'); 
        # plt.xlim([-3,3])
        # plt.ylim([-3,3])
        # # fig, ax = plt.subplots(7,1)
        # # for i in range(dt.shape[0]):
        # #     ax[i].plot(dt_train[i,:], dt_poly_train[i,:],'o')
            
        # plt.show()
        # import pdb; pdb.set_trace()

        dt_poly[~np.isfinite(dt_poly)] = 0
        s0_poly = dt_poly[0,:]
        dt_poly = dt_poly[1:,:]
        D_apprSq = 1/(np.sum(dt_poly[(0,3,5),:], axis=0)/3)**2
        if not flag_dti == 'True':
            dt_poly[6:,:] = dt_poly[6:,:] * np.tile(D_apprSq, (15,1))

        return dt_poly
    
    def create_dw_tensors(self, dt, order):
        """
        Function to take vectorized diffusion tensor (no s0) and generate 3x3xN D tensor
        and 3x3x3x3xN W tensor
        """

        if order == 2:
            DD = np.reshape(np.concatenate(
                (dt[0,:], dt[1,:], dt[2,:], 
                dt[1,:], dt[3,:], dt[4,:], 
                dt[2,:], dt[4,:], dt[5,:])
                ) ,(3, 3, dt.shape[1])) #.transpose(2,0,1)
            
            WT = None
        
        if order == 4:
            DD = np.reshape(np.concatenate(
                (dt[0,:], dt[1,:], dt[2,:], 
                dt[1,:], dt[3,:], dt[4,:], 
                dt[2,:], dt[4,:], dt[5,:])
                ) ,(3, 3, dt.shape[1])) #.transpose(2,0,1)
            
            wt = dt[6:,:].copy()
            WT = np.reshape(np.concatenate(
                (wt[0,:], wt[1,:], wt[2,:], wt[1,:], wt[3,:], wt[4,:], wt[2,:], wt[4,:], wt[5,:],
                wt[1,:], wt[3,:], wt[4,:], wt[3,:], wt[6,:], wt[7,:], wt[4,:], wt[7,:], wt[8,:],
                wt[2,:], wt[4,:], wt[5,:], wt[4,:], wt[7,:], wt[8,:], wt[5,:], wt[8,:], wt[5,:],
                wt[1,:], wt[3,:], wt[4,:], wt[3,:], wt[5,:], wt[7,:], wt[4,:], wt[7,:], wt[8,:],
                wt[3,:], wt[6,:], wt[4,:], wt[6,:], wt[10,:], wt[11,:], wt[7,:], wt[11,:], wt[12,:],
                wt[4,:], wt[7,:], wt[8,:], wt[7,:], wt[11,:], wt[12,:], wt[8,:], wt[12,:], wt[13,:],
                wt[2,:], wt[4,:], wt[5,:], wt[4,:], wt[7,:], wt[8,:], wt[5,:], wt[8,:], wt[9,:],
                wt[4,:], wt[7,:], wt[8,:], wt[7,:], wt[11,:], wt[12,:], wt[8,:], wt[12,:], wt[13,:],
                wt[5,:], wt[8,:], wt[9,:], wt[8,:], wt[12,:], wt[13,:], wt[9,:], wt[13,:], wt[14,:])
                ) ,(3, 3, 3, 3, dt.shape[1])) #.transpose(4,0,1,2,3)
            
        return DD, WT 

    
    def train_rotated_bayes_fit(self, dwi, dt, s0, b, mask, flag_dti='False'):
        """
        This function trains and evaluates a polynomial regression model to update
        a wlls fit without outlier voxels. 

        Inputs
        ------
        dwi: ND array 
        dt : wlls evaluated diffusion tensor Nx21 tensor
        s0 : first order element of dt
        b  : gradient matrix

        Outputs
        -------
        extended moments: dt fit from regression
        """

        SNR = 100
        sigma = 1 / SNR

        maxb = 3
        # grad_orig = self.grad
        # grad_keep = self.grad[:,3] < maxb
        # dwi = dwi[..., grad_keep]
        
        np.maximum(dwi, np.finfo(dwi.dtype).eps, out=dwi)

        D_apprSq = 1/(np.sum(dt[(0,3,5),:], axis=0)/3)**2
        if not flag_dti == 'True':
            dt[6:,:] /= np.tile(D_apprSq, (15,1))

        # first identify and remove outliers in dt
        if flag_dti == 'True':
            outlier_range = (1, 99)
            all_tinds = np.arange(6)
            trace = [0,3,5]
            nontrace = np.delete(all_tinds, trace)
            DD, WT = self.create_dw_tensors(dt, order=2)
        else:
            outlier_range = (5, 95)
            all_tinds = np.arange(21)
            trace = [0,3,5,6,9,11,16,18,20]
            nontrace = np.delete(all_tinds, trace)
            DD, WT = self.create_dw_tensors(dt, order=4)
        
        n_brain_voxels = dt.shape[1]
        n_coeffs = dt.shape[0]
        n_rotations = 10

        H = np.random.standard_normal((n_rotations,3,3))
        H = (H + H.transpose(0,2,1))/2
        (Evalues, Evec) = np.linalg.eig(H)
        Evec[0,:,:] = np.eye(3)
        dtx = np.zeros((n_coeffs, n_brain_voxels, n_rotations))
        #Dx = np.einsum('ilj,klm,inj->iknm', Evec, DD, Evec)

        print('rotations!!')
        for rot in range(n_rotations):
            R = Evec[rot,:,:]
  
            #Dx = (R @ (DD @ R.T)).transpose(1,2,0)
            Dx = np.einsum('ab...,ac...,bd...->cd...', 
                 DD, R, R)
            dtx[:6, :, rot] = np.vstack((Dx[0,0,:], Dx[0,1,:], Dx[0,2,:], Dx[1,1,:], Dx[1,2,:], Dx[2,2,:]))

            if not (flag_dti == 'True'):
                #Wx = (R @ (WT @ R.T)).transpose(1,2,3,4,0).reshape(9,9,n_brain_voxels)
                Wx = np.einsum('abcd...,ae,bf,cg,dh->efgh...', 
                    WT, R, R, R, R
                    ).reshape(9,9,n_brain_voxels, order='F')
                dtx[6:, :, rot] = np.vstack((Wx[0,0,:], Wx[0,1,:], Wx[0,2,:], Wx[0,4,:], Wx[0,5,:], 
                                             Wx[0,8,:], Wx[1,4,:], Wx[1,5,:], Wx[1,8,:], Wx[6,8,:],
                                             Wx[4,4,:], Wx[4,5,:], Wx[4,8,:], Wx[5,8,:], Wx[8,8,:])) 
                
        n_rotated_voxels = n_brain_voxels * n_rotations
        dtx = dtx.reshape(n_coeffs, n_rotated_voxels)        
        
        outlier_inds = np.zeros((len(trace), n_rotated_voxels))
        for i in range(len(trace)):
            outlier_inds[i,:] = (dtx[trace[i], :] < -0.1) | (dtx[trace[i], :] > 4)
        non_negative_trace = (1 - outlier_inds).all(axis=0)
        dt_ = dtx[:, non_negative_trace]

        outlier_inds = np.zeros((len(nontrace), np.sum(non_negative_trace)))
        for i in range(len(nontrace)):
            outlier_inds[i,:] = self.identify_outliers(dt_[nontrace[i],:], (2.5, 97.5))
        non_outlier_voxels = (1 - outlier_inds).all(axis=0)

        # interpolate non outliers in dt for more training data
        train_size = 400000 * n_rotations
        dt_init = dt_[:, non_outlier_voxels]
        scale_fact = train_size / dt_init.shape[1]

        #scale_fact = int(np.ceil(.5e6 / dt_init.shape[1]))
        dt_interp_init = zoom(dt_init, (1, scale_fact))
           
        # # set up training (zeroing diagnoal W elements)
        # train_size = 400000
        # id_non_outliers = np.random.permutation(train_size)
        # dt_train = np.zeros((dt.shape[0], train_size))
        # #inds_exclude = np.array([9,10,12,14,15,16,17,19,21]) - 1
        # for id in range(dt.shape[0]):
        #     #if np.any(id == inds_exclude):
        #     #    dt_train[id,:] = 0
        #     #else:
        #     dt_train[id,:] = np.real((dt_interp_init[id, id_non_outliers]) * 1.0)


        dt_train = dt_interp_init.copy()
        
        # dt_train[0,:] = -1 * np.random.exponential(0.707, (1,train_size))
        #dt_train[0,:] = np.random.uniform(0,1,(1, train_size)) * 2

        # ivim stuff
        dt_train = np.vstack((
            np.random.uniform(0,1,(1, train_size)) * 3,
            dt_train))
        
        f_pri = np.random.uniform(0,1,(1, train_size)) * 0.6
        D_pri = np.random.uniform(0,1,(1, train_size)) * 200 + 3
        
        bval = self.grad[:,-1]
        s_training = dt_train[:,0] * ((1-f_pri) * np.exp(b @ dt_train).T +
                                      f_pri * np.exp(-bval * D_pri))
        s_training += np.abs(sigma *
                        (np.random.standard_normal(s_training.shape) + 
                        1j*np.random.standard_normal(s_training.shape)))

        # # original prior stuff
        # dt_train = np.vstack((
        #         -1 * np.random.exponential(0.707, (1,train_size)),
        #         dt_train))
        # s_training = np.exp(b @ dt_train).T
        # s_training += np.abs(sigma *
        #                 (np.random.standard_normal(s_training.shape) + 
        #                 1j*np.random.standard_normal(s_training.shape)))
        
        # import matplotlib.pyplot as plt
        # fig, axs = plt.subplots(3,7)
        # axs = axs.ravel()
        # for i in range(21):
        #     axs[i].hist(dt_train[i,:],100)
        # plt.show()
        # import pdb; pdb.set_trace()
        
        print('PR')
        # compress the training data for speed
        npars = 15
        x = np.log(s_training)
        #u,s,v = np.linalg.svd((x - x.mean()) / x.std(), full_matrices=False)
        #reduced_dim = v[:,npars:]
        data_training = x #@ reduced_dim

        # polynomial regression
        order_poly_fit = 1
        x_moments = self.compute_extended_moments(data_training, order_poly_fit)
        pinv_x = np.linalg.pinv(x_moments)
        
        coeffs_dt = np.zeros((x_moments.shape[1], dt_train.shape[0]))
        for i in range(dt_train.shape[0]):
            coeffs_dt[:,i] = pinv_x @ dt_train[i,:].T

        # apply the trained coefficients to the original data
        dwi = gaussian_filter(dwi, sigma=(0.425,.425,.425,0), truncate=2)
        data_testing = np.log(vectorize(dwi, mask) / s0).T #@ reduced_dim
        x_moments_test = self.compute_extended_moments(data_testing, order_poly_fit)
        dt_poly = np.zeros((dt_train.shape[0], n_brain_voxels))
        #dt_poly_train = np.zeros_like(dt_train)
        for i in range(dt_train.shape[0]):
            kn = x_moments_test @ coeffs_dt[:,i]
            kn[kn < np.min(dt_train[i,:])] = np.min(dt_train[i,:])
            kn[kn > np.max(dt_train[i,:])] = np.max(dt_train[i,:])
            dt_poly[i,:] = kn
            #dt_poly_train[i,:] = x_moments @ coeffs_dt[:,i]

        # import matplotlib.pyplot as plt
        # plt.plot(dt_train[1,:],dt_poly_train[1,:],'o'); 
        # plt.xlim([-3,3])
        # plt.ylim([-3,3])
        # # fig, ax = plt.subplots(7,1)
        # # for i in range(dt.shape[0]):
        # #     ax[i].plot(dt_train[i,:], dt_poly_train[i,:],'o')
            
        # plt.show()
        # import pdb; pdb.set_trace()

        dt_poly[~np.isfinite(dt_poly)] = 0
        s0_poly = dt_poly[0,:]
        dt_poly = dt_poly[1:,:]
        D_apprSq = 1/(np.sum(dt_poly[(0,3,5),:], axis=0)/3)**2
        if not flag_dti == 'True':
            dt_poly[6:,:] *= np.tile(D_apprSq, (15,1))

        return dt_poly
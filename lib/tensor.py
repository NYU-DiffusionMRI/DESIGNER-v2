from joblib import Parallel, delayed
import numpy as np
from tqdm import tqdm
from lib.mpunits import vectorize
import os

import scipy.io as sio
import cvxpy as cp

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
        Compute the approximate diffusion coefficient
        """
        # compute ADC
        dcnt, dind = self.create_tensor_order(2)
        ndir = dir.shape[0]
        bD = np.tile(dcnt,(ndir, 1)) * dir[:,dind[:, 0]] * dir[:,dind[:, 1]]
        adc = bD @ dt
        return adc


    def kurtosis_coeff(self, dt, dir):
        """
        Function to Compute the aproximate kurtosis coefficient
        """
        # compute AKC
        wcnt, wind = self.create_tensor_order(4)
        ndir = dir.shape[0]
        adc = self.diffusion_coeff(dt[:6], dir)
        md = np.sum(dt[np.array([0,3,5])], 0)/3
        bW = np.tile(wcnt,(ndir, 1)) * dir[:,wind[:, 0]] * dir[:,wind[:, 1]] * dir[:,wind[:, 2]] * dir[:,wind[:, 3]]
        akc = bW @ dt[6:]
        akc = (akc * np.tile(md**2, (adc.shape[0], 1))) / (adc**2)
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
        theta = np.arange(0,2*np.pi-dt,dt)
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
                dt = np.linalg.pinv(w @ b) @ (w @ np.log(dwi))
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
            dt = np.zeros((22))

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
            dirs = np.array(self.fibonacci_sphere(256, True))
            
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
        Diffusion parameter estimation using weightel linear least squares fitting
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
        nblocks = 10

        try:
            N = 10000
            akc_mask = Parallel(n_jobs=self.n_cores, prefer='processes')\
                (delayed(self.compute_outliers)
                (dt, dir[int(N/nblocks*(i-1)):int(N/nblocks*i),:]) for i in range(1, nblocks + 1)
                )
        except:
            N = 5000
            akc_mask = Parallel(n_jobs=8, prefer='processes')\
                (delayed(self.compute_outliers)
                (dt, dir[int(N/nblocks*(i-1)):int(N/nblocks*i),:]) for i in range(1, nblocks + 1)
                )

        return np.sum(akc_mask, axis=0)



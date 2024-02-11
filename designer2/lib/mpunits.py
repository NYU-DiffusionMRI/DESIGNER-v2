from re import L
import numpy as np
import os
from scipy.linalg import svd
from numpy_groupies import aggregate
from skimage.util.shape import view_as_windows


def vectorize(image, mask):
    """
    - If the input is 1D or 2D: unpatch it to 3D or 4D using a mask.
    - If the input is 3D or 4D: vectorize to 2D using a mask.
    - Used by all classes so adding as global function
    """
    if mask is None:
        mask = np.ones((image.shape[0], image.shape[1], image.shape[2]), order='F')

    if image.ndim == 1:
        image = np.expand_dims(image, axis=0)

    if image.ndim == 2:
        n = image.shape[0]
        s = np.zeros((mask.shape[0], mask.shape[1], mask.shape[2], n), order='F')
        for i in range(0, n):
            dummy = np.zeros((mask.shape))
            dummy[mask] = image[i,:]
            s[:,:,:,i] = dummy

    if image.ndim == 3:
        image = np.expand_dims(image, axis=-1)

    if image.ndim == 4:
        s = np.zeros((image.shape[-1], np.sum(mask).astype(int)), order='F')
        for i in range(0, image.shape[-1]):
            tmp = image[:,:,:,i]
            s[i,:] = tmp[mask]

    return np.squeeze(s)


def unpad(x, pad_width):
    """
    If there is padding present, remove it.
    """
    slices = []
    for c in pad_width:
        e = None if c[1] == 0 else -c[1]
        slices.append(slice(c[0], e))
    return x[tuple(slices)]


def get_image_dimensions(X):
    """
    M and generally the number of pixels in a patch, N is the number of image volumes
    M and N will vary from case to case, as we generally create a patch with minimum size:
    product(patch_shape) > n measurements
    """
    M = X.shape[0]
    N = X.shape[1]
    
    return (M, N)


def maybe_transpose(M, N, X):
    """
    Compute conjugate matrix and transpose if needed.
    """
    if M < N:
        X = np.conj(X).T
        
    return X


def compute_svd(X):
    """
    Compute singular value decomposition on matrix,
    square result and convert to float32.
    """
    u, vals, v = svd(X, full_matrices=False)
    vals = (vals**2).astype('float32')
    
    return u, v, vals


def select_msv(vals, u, v):
    """
    Select most significant SVD value:
    We are sorting the singular values from largest to smallest, the largest singular 
    value generally corresponds to the largest "thing" in an image that corresponds to 
    signal (often the mean of the image). The smallest SVs will correspond to Gaussian Noise,
    the goal of this function is to automatically detect the cuttoff between noise carrying 
    and signal carryings SVs.
    """
    order = np.argsort(vals)[::-1]
    
    u = u[:,order]
    v = v[:,order]
    vals = vals[order]
    
    return u, v, vals


def init_crop_factor(Mp):
    """
    tn is a crop factor we are not using (cropping out the very smallest SVs) 
    In the case of heteroscedastic noise, There can be a small number of non-noise
    singular values smaller than noise carrying ones. It is sometimes nice to
    therefore include a "crop factor".
    """
    tn = 0
    return (0, np.arange(0,Mp-tn).T)


def create_signal_carrying_vector(Mp):
    """
    Create a vector of potential "signal carrying SVs".
    """
    return np.arange(0,Mp).T


def get_csum(Mp, vals):
    """
    Vector of cumulative sums of singular values.
    """
    return np.cumsum(vals[::-1])[::-1]


def find_signal_carrying_svs(Mp, tn, sigmasq_1, sigmasq_2):
    """
    Find the intersection between cumsum noise variance and MP noise variance.
    If the two curves intersect, let the first intersection correspond to the noisy SV cutoff.
    Otherwise, keep all singular values assuming that none correspond to noise.
    """
    t = np.where(sigmasq_2 < sigmasq_1[:Mp-tn])

    if t[0].any():
        t = t[0][0]
    else:
        t = Mp - 1

    return t


def get_sigma(sigmasq_1, t):
    """
    Sigma is the square root of noise variance at index t.
    """
    return np.sqrt(sigmasq_1[t])


def get_default_kernel(nvols):
    """
    Return a default kernel if not specified in the module configuration. Max size is 9x9x9
    """
    p = np.arange(3, nvols, 2)
    pf = np.where(p**3 >= nvols)[0][0]
    defaultKernel = p[pf]
    return np.array([defaultKernel, defaultKernel, defaultKernel])


def is_diffusion_shelled(grad):
    """
    Determine the number of unique shells in a diffusion series and ensure
    that there is at least one b=0 and one b>0 shell
    """
    b0inds = grad[:,-1] == 0
    if sum(b0inds) > 0:
        at_least_one_bzero = True
    else:
        at_least_one_bzero = False

    bgt0inds = grad[:,-1] > 0
    if sum(bgt0inds) > 0:
        at_least_one_nonbzero = True
    else:
        at_least_one_nonbzero = False

    if at_least_one_bzero and at_least_one_nonbzero:
        return True
    else:
        return False


def compute_noise_variance_vector(vals, Mp, tn, rangeMP):
    """
    Compute vector of noise variance according to Marchneco Pastur theory
    """
    rangeData = vals[:Mp-tn] - vals[Mp-1]
    sigmasq_2 = rangeData / rangeMP
        
    return sigmasq_2


def veraart(Mp, Np, tn, ptn, p, csum):
    """
    Compute veraart
    """
    sigmasq_1 = csum / ((Mp - p) * Np)
    rangeMP = 4 * np.sqrt((Mp - ptn) * (Np - tn))
    
    return (sigmasq_1, rangeMP)


def cordero_grande(Mp, Np, ptn, p, csum):
    """
    Compute Cordero-Grande
    """
    sigmasq_1 = csum / ((Mp - p) * (Np - p))
    rangeMP = 4 * np.sqrt((Mp - ptn) * (Np - ptn))
    
    return (sigmasq_1, rangeMP)


def jespersen(Mp, Np, tn, p, csum):
    """
    Compute Jespersen
    """
    sigmasq_1 = csum / ((Mp - p) * (Np - p))
    rangeMP = 4 * np.sqrt((Np - tn) * (Mp))
    
    return (sigmasq_1, rangeMP)


def eig_shrink(vals, gamma):
    """
    eig_shrink
    """
    t = 1 + np.sqrt(gamma)
    s = np.zeros((vals.shape))
    x = vals[vals > t]
    s[vals > t] = np.sqrt((x**2 - gamma - 1)**2 - 4*gamma) / x
    return np.diag(s)


def crop_image(vals, Mp, Np, u, v, sigma):
    """
    Crop SVs corresponding to noise directly.
    """
    vals_norm = np.sqrt(vals)/(np.sqrt(Mp) * sigma)
    vals_frobnorm = (np.sqrt(Mp)*sigma) * eig_shrink(vals_norm, Np/Mp)
    return np.matrix(u) * vals_frobnorm * np.matrix(v)


def shrink_image(vals, u, v, t):
    """
    Use shrinkage to minimize the influence of noise carrying SVs proportionally.
    """
    vals[t:] = 0
    return np.matrix(u) * np.diag(np.sqrt(vals)) * np.matrix(v)


def compute_weighted_patchav(wp, signal, temp, nvoxels, dims):
    """
    Weighted patch averaging is performed to account for the contribution of
    voxels that occur in multiple patches.
    """
    return aggregate(
        temp, wp*signal, func='sum', size=nvoxels, dtype=signal.dtype
    ).reshape(dims)


def get_weighted_patch_average(temp, dims, nvoxels, kernel, step):
    """
    Get distance based wights for patch averaging.
    """

    (i, j, k) = np.unravel_index(temp, dims)

    distance = (i - np.mean(i, axis=1, keepdims=True)) ** 2 + \
        (j - np.mean(j, axis=1, keepdims=True)) ** 2 + \
        (k - np.mean(k, axis=1, keepdims=True)) ** 2

    count = np.histogram(temp, np.arange(nvoxels + 1))[0]

    denominator = view_as_windows(
        compute_weighted_patchav(1, distance.flatten(), temp.flatten(), nvoxels, dims),
        kernel, step
    ).reshape(-1, np.prod(kernel))
    
    cbl = view_as_windows(
        count.reshape(dims), kernel, step
    ).reshape(-1, np.prod(kernel))

    num = denominator - distance
    wp = num / (denominator * (cbl - 1) + np.finfo(float).eps)
    wp[cbl == 1] = 1
    wp = wp.flatten()

    return wp

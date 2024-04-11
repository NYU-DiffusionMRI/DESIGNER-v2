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

import pytest
import numpy as np
from designer_fit_wrappers import parallel_outlier_smooth, refit_or_smooth

def test_parallel_outlier_smooth():
    # Replace with actual test data
    inds = (1, 1, 1)
    kernel = 3
    outlier_locations = np.zeros((5, 5, 5))
    dwi_norm = np.ones((5, 5, 5, 3))
    dwi = np.ones((5, 5, 5, 3))
    smoothlevel = None
    
    result = parallel_outlier_smooth(inds, kernel, outlier_locations, dwi_norm, dwi, smoothlevel)
    # Add assertions to verify the expected behavior
    assert result is not None

def test_refit_or_smooth():
    # Replace with actual test data
    outlier_locations = np.zeros((5, 5, 5))
    dwi = np.ones((5, 5, 5, 3))
    mask = None
    smoothlevel = None
    n_cores = 1
    
    result = refit_or_smooth(outlier_locations, dwi, mask, smoothlevel, n_cores)
    # Add assertions to verify the expected behavior
    assert result.shape == dwi.shape
import numpy as np
from scipy import ndimage
from skimage.transform import resize
from joblib import Parallel, delayed
from numba import njit, prange
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans


def normalize_data(data):
    """
    Vectorized normalization of 4D MRI data based on Median Absolute Deviation (MAD) and global maximum.

    Parameters:
    - data: 4D numpy array [Nx, Ny, Nz, Nt]

    Returns:
    - data_normalized: Normalized data
    - norm_factor: Array of MAD factors per time frame
    - max_norm: Global maximum for normalization
    """
    # Ensure data is in float format to handle divisions properly
    data = data.astype(np.float32, copy=False)
    
    masked_data = np.where(data != 0, data, np.nan)

    medians = np.nanmedian(masked_data, axis=(0, 1, 2))

    # Compute per-frame MAD: median of absolute deviations from the median
    abs_devs = np.abs(masked_data - medians[np.newaxis, np.newaxis, np.newaxis, :])
    mad = np.nanmedian(abs_devs, axis=(0, 1, 2))
    mad = np.where(np.isnan(mad) | (mad == 0), 1.0, mad)
    data_normalized = data / mad[np.newaxis, np.newaxis, np.newaxis, :]

    max_norm = np.nanmax(data_normalized)
    if max_norm == 0:
        max_norm = 1.0

    data_normalized /= max_norm

    return data_normalized, mad, max_norm

def denormalize_data(data, norm_factor, max_norm):
    """
    De-normalize the data based on MAD factors and global maximum.

    Parameters:
    - data: Normalized 4D data
    - norm_factor: Array of MAD factors per time frame
    - max_norm: Global maximum used in normalization

    Returns:
    - data_denormalized: De-normalized data
    """
    # Reshape norm_factor for broadcasting
    norm_factor_reshaped = norm_factor.reshape(1, 1, 1, -1)    
    data *= max_norm * norm_factor_reshaped

    return data

def downsample_data_3d(data, factor, n_jobs=-1):
    """
    Downsample 4D data spatially by the given factor.

    Parameters:
    - data: 4D numpy array [Nx, Ny, Nz, Nt]
    - factor: Downsampling factor (integer)

    Returns:
    - data_ds: Downsampled data
    """
    if factor == 1:
        return data.copy()
    
    Nx, Ny, Nz, Nt = data.shape
    new_size = (
        max(1, round(Nx / factor)),
        max(1, round(Ny / factor)),
        max(1, round(Nz / factor))
    )
    
    def resize_frame(t):
        return resize(
            data[..., t],
            new_size,
            order=3,  # cubic interpolation
            mode='edge',
            anti_aliasing=True,
            preserve_range=True  # Preserve original data range
        )
    
    # Perform parallel resizing
    resized_frames = Parallel(n_jobs=n_jobs)(
        delayed(resize_frame)(t) for t in range(Nt)
    )
    
    # Stack resized frames along the last axis
    data_ds = np.stack(resized_frames, axis=-1).astype(data.dtype)
    
    return data_ds

def upsample_data_3d(data, target_size, n_jobs=-1):
    """
    Upsample 4D data spatially to the target size using parallel processing.

    Parameters:
    - data: 4D numpy array [Nx, Ny, Nz, Nt]
    - target_size: Tuple of target spatial dimensions (Nx_new, Ny_new, Nz_new)
    - n_jobs: Number of parallel jobs (-1 uses all available cores)

    Returns:
    - data_us: Upsampled data
    """
    Nx_new, Ny_new, Nz_new = target_size
    Nx, Ny, Nz, Nt = data.shape

    if Nt == 0:
        raise ValueError("The input data has zero time frames (Nt=0).")

    # Define a function to resize a single time frame
    def resize_frame(t):
        return resize(
            data[..., t],
            target_size,
            order=3,  # cubic interpolation
            mode='edge',
            anti_aliasing=True,
            preserve_range=True  # Preserve original data range
        )

    # Perform parallel resizing
    resized_frames = Parallel(n_jobs=n_jobs)(
        delayed(resize_frame)(t) for t in range(Nt)
    )

    # Stack resized frames along the last axis
    data_us = np.stack(resized_frames, axis=-1).astype(data.dtype)

    return data_us

def warp_data_3d(I, u, v, w, n_jobs=-1):
    """
    Warps the 4D image I using the displacement fields u, v, w using parallel processing.

    Parameters:
    - I: 4D numpy array [Nx, Ny, Nz, Nt]
    - u, v, w: 4D displacement fields [Nx, Ny, Nz, Nt]
    - n_jobs: Number of parallel jobs (-1 uses all available cores)

    Returns:
    - warped_image: Warped 4D image
    """
    Nx, Ny, Nz, Nt = I.shape
    warped_image = np.empty_like(I)
    
    # Create meshgrid for coordinates
    X, Y, Z = np.meshgrid(np.arange(Ny), np.arange(Nx), np.arange(Nz), indexing='ij')
    
    def warp_frame(t):
        X_new = X + u[..., t]
        Y_new = Y + v[..., t]
        Z_new = Z + w[..., t]
        
        # scipy's map_coordinates expects the coordinates in the order (z, y, x)
        coords = [Z_new.ravel(), Y_new.ravel(), X_new.ravel()]
        warped = ndimage.map_coordinates(
            I[..., t],
            coords,
            order=3,
            mode='nearest'
        )
        return warped.reshape(Nx, Ny, Nz)
    
    # Perform parallel warping
    resized_frames = Parallel(n_jobs=n_jobs)(
        delayed(warp_frame)(t) for t in range(Nt)
    )
    
    # Stack resized frames along the last axis
    warped_image = np.stack(resized_frames, axis=-1)
    
    return warped_image

def unpad(data, kernel):
    """
    Removes padding from the data based on kernel size.

    Parameters:
    - data: 4D numpy array [Nx, Ny, Nz, Nt]
    - kernel: List or array of kernel sizes [kx, ky, kz]

    Returns:
    - data_unpadded: Unpadded data
    """
    k = (np.array(kernel) - 1) // 2
    return data[
        int(k[0]):-int(k[0]) if k[0] !=0 else None,
        int(k[1]):-int(k[1]) if k[1] !=0 else None,
        int(k[2]):-int(k[2]) if k[2] !=0 else None,
        :
    ]

def process_patch_indices(patch_pos, blocksize, stepsize, I_shape, I):
    """
    Extracts a single patch from the 4D image and computes its linear indices.
    
    Parameters:
    - patch_pos: Tuple (x_start, y_start, z_start, t_start)
    - blocksize: Tuple (nrows, ncols, nslcs, ntimes)
    - stepsize: Tuple (d_row, d_col, d_slc, d_time)
    - I_shape: Shape of the 4D image [m, n, r, t]
    - I: 4D numpy array [x, y, z, t]
    
    Returns:
    - patch_flat: Flattened patch data [nrows * ncols * nslcs * (ntimes + 1), ]
    - combined_indices: Flattened linear indices corresponding to patch elements [nrows * ncols * nslcs * (ntimes + 1), ]
    """
    nrows, ncols, nslcs, ntimes = blocksize
    m, n_dim, r, t_dim = I_shape
    
    x_start, y_start, z_start, t_start = patch_pos
    # Compute end positions, ensuring they do not exceed dimensions
    x_end = min(x_start + nrows, m)
    y_end = min(y_start + ncols, n_dim)
    z_end = min(z_start + nslcs, r)
    t_end = min(t_start + ntimes, t_dim)
    
    actual_nrows = x_end - x_start
    actual_ncols = y_end - y_start
    actual_nslcs = z_end - z_start
    actual_ntimes = t_end - t_start
    
    # Extract reference patch (frame 0)
    ref_patch = I[x_start:x_start+actual_nrows, y_start:y_start+actual_ncols, z_start:z_start+actual_nslcs, 0]  # Shape: [nrows, ncols, nslcs]
    
    # Extract current temporal patches
    curr_patch = I[x_start:x_start+actual_nrows, y_start:y_start+actual_ncols, z_start:z_start+actual_nslcs, t_start:t_start+actual_ntimes]  # Shape: [nrows, ncols, nslcs, ntimes]
    
    # Concatenate along temporal axis
    combined_patch = np.concatenate((ref_patch[..., np.newaxis], curr_patch), axis=3)  # Shape: [nrows, ncols, nslcs, ntimes + 1]
    
    # Flatten the combined patch
    patch_flat = combined_patch.flatten()
    
    # Compute linear indices for reference patch
    # In NumPy, linear indices are row-major (C order)
    # Compute indices for [x_start:x_end, y_start:y_end, z_start:z_end, t=0]
    ref_indices = np.ravel_multi_index(
        (
            x_start + np.arange(actual_nrows)[:, None, None],
            y_start + np.arange(actual_ncols)[None, :, None],
            z_start + np.arange(actual_nslcs)[None, None, :]
        ),
        dims=I_shape[:3]
    ).flatten()
    
    # Compute linear indices for current patches
    # Formula: linear_index = x * n * r * t + y * r * t + z * t + t_index
    # where x, y, z are spatial indices and t_index is the temporal index
    x_indices = x_start + np.arange(actual_nrows).reshape(actual_nrows, 1, 1)
    y_indices = y_start + np.arange(actual_ncols).reshape(1, actual_ncols, 1)
    z_indices = z_start + np.arange(actual_nslcs).reshape(1, 1, actual_nslcs)
    
    # Compute spatial linear indices
    spatial_lin = (
        x_indices * n_dim * r * t_dim +
        y_indices * r * t_dim +
        z_indices * t_dim
    )
    
    # Expand spatial_lin to include temporal frames
    spatial_lin_expanded = spatial_lin[:, :, :, np.newaxis]  # (nrows, ncols, nslcs, 1)
    
    # Compute temporal indices
    t_indices = t_start + np.arange(actual_ntimes).reshape(1, 1, 1, actual_ntimes)
    
    # Compute linear indices for current patches
    linear_indices_curr = spatial_lin_expanded + t_indices  # Shape: (nrows, ncols, nslcs, ntimes)
    
    # Flatten the current linear indices
    linear_indices_curr_flat = linear_indices_curr.flatten()
    
    # Combine reference and current indices
    combined_indices = np.concatenate((ref_indices, linear_indices_curr_flat))
    
    return patch_flat, combined_indices

def im2col4(I, blocksize, stepsize, n_jobs=-1):
    """
    Optimized version of im2col4, extracting 4D patches from a 4D image into columns.
    Includes the reference frame and ensures unique frames within each patch.
    
    Parameters:
    - I: 4D numpy array [x, y, z, t].
    - blocksize: Tuple (nrows, ncols, nslcs, ntimes), size of the patches.
    - stepsize: Tuple (d_row, d_col, d_slc, d_time), steps between patches.
    - n_jobs: Number of parallel jobs (-1 uses all available cores).
    
    Returns:
    - out: 2D numpy array where each column is a vectorized patch from the image.
    - index_patches: 2D numpy array of linear indices for each patch.
    """
    nrows, ncols, nslcs, ntimes = blocksize
    d_row, d_col, d_slc, d_time = stepsize
    m, n_dim, r, t_dim = I.shape

    # Adjust ntimes to include the reference frame
    total_frames_in_patch = ntimes + 1  # +1 for the reference frame

    # Calculate the number of patches along each dimension
    num_patches_x = (m - nrows) // d_row + 1
    num_patches_y = (n_dim - ncols) // d_col + 1
    num_patches_z = (r - nslcs) // d_slc + 1
    num_patches_t = (t_dim - ntimes) // d_time + 1

    # Generate all patch starting positions
    patch_positions_x = np.arange(0, num_patches_x * d_row, d_row)
    patch_positions_y = np.arange(0, num_patches_y * d_col, d_col)
    patch_positions_z = np.arange(0, num_patches_z * d_slc, d_slc)
    patch_positions_t = np.arange(1, num_patches_t * d_time + 1, d_time)  # Start from frame 1 (zero-based)

    # Ensure the last patch includes the edge by adjusting if necessary
    if patch_positions_t[-1] + ntimes > t_dim:
        patch_positions_t[-1] = t_dim - ntimes

    # Create a meshgrid of all patch starting positions
    grid_x, grid_y, grid_z, grid_t = np.meshgrid(
        patch_positions_x,
        patch_positions_y,
        patch_positions_z,
        patch_positions_t,
        indexing='ij'
    )

    # Flatten the grid to obtain a list of starting positions
    patch_positions = np.stack([
        grid_x.flatten(),
        grid_y.flatten(),
        grid_z.flatten(),
        grid_t.flatten()
    ], axis=1)  # Shape: (num_patches, 4)

    num_patches = patch_positions.shape[0]
    psize = nrows * ncols * nslcs * total_frames_in_patch

    # Preallocate output arrays
    out = np.empty((psize, num_patches), dtype=I.dtype)
    index_patches = np.empty((psize, num_patches), dtype=np.int64)

    # Define the worker function
    def worker(patch_pos):
        patch_flat, combined_indices = process_patch_indices(patch_pos, blocksize, stepsize, I.shape, I)
        return patch_flat, combined_indices

    # Parallel processing of patches
    results = Parallel(n_jobs=n_jobs)(
        delayed(worker)(pos) for pos in patch_positions
    )

    # Assemble the results
    for idx, (patch_flat, combined_indices) in enumerate(results):
        out[:, idx] = patch_flat
        index_patches[:, idx] = combined_indices

    return out, index_patches

def compute_gradients(patch):
    """
    Computes spatial and temporal gradients of the input data.

    Parameters:
    - patch: 4D numpy array [psize_x, psize_y, psize_z, N]

    Returns:
    - Gx, Gy, Gz, It: Spatial and temporal gradients
    """

    # Compute spatial gradients using central differences with vectorization
    Gx = (np.roll(patch, -1, axis=1) - np.roll(patch, 1, axis=1)) / 2.0
    Gy = (np.roll(patch, -1, axis=0) - np.roll(patch, 1, axis=0)) / 2.0
    Gz = (np.roll(patch, -1, axis=2) - np.roll(patch, 1, axis=2)) / 2.0

    # Handle boundary conditions by setting edge gradients to zero
    Gx[:, 0, :, :] = 0
    Gx[:, -1, :, :] = 0
    Gy[0, :, :, :] = 0
    Gy[-1, :, :, :] = 0
    Gz[:, :, 0, :] = 0
    Gz[:, :, -1, :] = 0

    # Compute temporal gradients using vectorized slicing
    It = np.zeros_like(patch)
    It[..., 1:] = patch[..., 1:] - patch[..., :-1]
    It[..., 0] = 0  # Optionally set to zero or replicate the first temporal difference

    return Gx, Gy, Gz, It

def design_matrix_for_patch(gx_patch, gy_patch, gz_patch, gt_patch, brightness_order):
    """
    Constructs the design matrix A and vector b for motion estimation.

    Parameters:
    - gx_patch, gy_patch, gz_patch, gt_patch: Gradient arrays
    - brightness_order: Integer for brightness modeling

    Returns:
    - A: Design matrix
    - b: Vector
    """
    N = gx_patch.shape[3]
    psize = gx_patch.shape[0] * gx_patch.shape[1] * gx_patch.shape[2]
    num_brightness_terms = brightness_order + 1
    total_terms_per_frame = 3 + num_brightness_terms

    # Preallocate A and b
    A = np.zeros((psize * (N - 1), total_terms_per_frame * (N - 1)), dtype=np.float32)
    b = np.zeros((psize * (N - 1), 1), dtype=np.float32)

    # Extract all necessary frames at once
    gx = gx_patch[..., 1:N].reshape(psize, N - 1)  
    gy = gy_patch[..., 1:N].reshape(psize, N - 1) 
    gz = gz_patch[..., 1:N].reshape(psize, N - 1)  
    gt = gt_patch[..., 1:N].reshape(psize, N - 1)  

    # Compute brightness terms for all frames
    E = gt[:, :, np.newaxis] ** np.arange(num_brightness_terms)  # Shape: (psize, N-1, num_brightness_terms)
    E = E.reshape(psize, (N - 1) * num_brightness_terms)       # Shape: (psize, (N-1)*num_brightness_terms)

    # Assign b
    b[:, 0] = -gt.flatten()

    # Assign gradients and brightness terms to A
    for k in range(N - 1):
        row_start = k * psize
        row_end = (k + 1) * psize
        col_start = k * total_terms_per_frame
        col_end = (k + 1) * total_terms_per_frame

        # Assign gx, gy, gz
        A[row_start:row_end, col_start] = gx[:, k]
        A[row_start:row_end, col_start + 1] = gy[:, k]
        A[row_start:row_end, col_start + 2] = gz[:, k]

        # Assign brightness terms
        A[row_start:row_end, col_start + 3:col_start + 3 + num_brightness_terms] = E[:, k*num_brightness_terms:(k+1)*num_brightness_terms].reshape(psize, num_brightness_terms)

    return A, b

def compute_motion_fields_for_patch(coeffs, patch_size, brightness_order):
    """
    Computes the motion fields U, V, W for a patch from coefficients using vectorized operations.

    Parameters:
    - coeffs: 1D numpy array, length = (3 + brightness_order + 1)*(N-1)
    - patch_size: Tuple [psize_x, psize_y, psize_z, N]
    - brightness_order: Integer for brightness modeling

    Returns:
    - U, V, W: 4D numpy arrays [psize_x, psize_y, psize_z, N]
    """
    psize_x, psize_y, psize_z, N = patch_size
    num_brightness_terms = brightness_order + 1
    total_terms_per_transition = 3 + num_brightness_terms

    # Validate the length of coeffs
    expected_length = total_terms_per_transition * (N - 1)
    if coeffs.size != expected_length:
        raise ValueError(f"Coefficient vector length {coeffs.size} does not match expected {(total_terms_per_transition)*(N-1)}.")

    # Initialize U, V, W with zeros
    U = np.zeros((psize_x, psize_y, psize_z, N), dtype=np.float32)
    V = np.zeros_like(U)
    W = np.zeros_like(U)

    # Reshape coeffs to (N-1, total_terms_per_transition)
    coeffs_reshaped = coeffs.reshape(N-1, total_terms_per_transition)

    # Extract u, v, w coefficients
    u_coeffs = coeffs_reshaped[:, 0].reshape(1, 1, 1, N-1)  # Shape: (1, 1, 1, N-1)
    v_coeffs = coeffs_reshaped[:, 1].reshape(1, 1, 1, N-1)  # Shape: (1, 1, 1, N-1)
    w_coeffs = coeffs_reshaped[:, 2].reshape(1, 1, 1, N-1)  # Shape: (1, 1, 1, N-1)

    # Assign the coefficients to U, V, W for frames 1 to N-1
    U[..., 1:N] = u_coeffs
    V[..., 1:N] = v_coeffs
    W[..., 1:N] = w_coeffs

    return U, V, W

def get_mask(data, k=5):
    """
    Generates a mask excluding the reference cluster using k-means clustering.

    Parameters:
    - data: 4D numpy array [Nx, Ny, Nz, Nd]
    - k: Number of clusters

    Returns:
    - mask: 3D boolean numpy array [Nx, Ny, Nz]
    """
    Nx, Ny, Nz, Nd = data.shape
    reshaped_data = data.reshape(-1, Nd)

    # Perform k-means clustering
    kmeans = KMeans(n_clusters=k, max_iter=1000, n_init=5, random_state=42)
    km = kmeans.fit_predict(reshaped_data)

    # Select the largest cluster as reference
    unique, counts = np.unique(km, return_counts=True)
    reference_cluster = unique[np.argmax(counts)]

    # Create mask
    mask = km != reference_cluster
    mask = mask.reshape(Nx, Ny, Nz)

    return mask

def process_patch(patch_flat, brightness_order):
    """
    Wrapper function to process a flattened patch.

    Parameters:
    - patch_flat: 1D numpy array
    - brightness_order: Integer

    Returns:
    - du_patch, dv_patch, dw_patch
    """
    # Reshape patch
    # Placeholder: define actual patch shape
    patch = patch_flat.reshape((5, 5, 5, -1))  # Example shape

    # Compute gradients
    Gx, Gy, Gz, It = compute_gradients(patch)

    # Design matrix
    A, b = design_matrix_for_patch(Gx, Gy, Gz, It, brightness_order)

    # Regularization parameter
    lambda_reg = 1.0

    # Solve (A'A + lambda I) c = A'b using numpy's solver
    AtA = A.T @ A + lambda_reg * np.eye(A.shape[1], dtype=A.dtype)
    Atb = A.T @ b
    coeffs = np.linalg.solve(AtA, Atb)

    # Compute motion fields
    U_patch, V_patch, W_patch = compute_motion_fields_for_patch(coeffs, patch.shape, brightness_order)

    # Flatten the motion fields
    return U_patch.flatten(), V_patch.flatten(), W_patch.flatten()

def offast_lkv_3d_4d_patchavg_pyramid_v3(
    data,
    kernel=[5, 5, 5],
    step=[3, 3, 3],
    brightness_order=1,
    window_size=5,
    step_time=2,
    pyramid_levels=[4, 2],
    iterations_per_level=[3, 3]
):
    """
    Optical flow estimation for 4D MRI data using a multi-scale approach.

    Parameters:
    - data: 4D numpy array [Nx, Ny, Nz, Nt]
    - kernel: List of kernel sizes [kx, ky, kz]
    - step: List of step sizes [sx, sy, sz]
    - brightness_order: Integer for brightness modeling
    - window_size: Integer for window size
    - step_time: Integer for time step
    - pyramid_levels: List of downsampling factors
    - iterations_per_level: List of iterations per pyramid level

    Returns:
    - WarpedData: Corrected 4D data
    - U, V, W: Motion fields in x, y, z directions
    - COST: Cost metric per iteration and level
    """

    if isinstance(iterations_per_level, int):
        iterations_per_level = [iterations_per_level] * len(pyramid_levels)
    if len(iterations_per_level) != len(pyramid_levels):
        raise ValueError("iterations_per_level must be a scalar or have the same length as pyramid_levels.")

    adjusted_window_size = window_size + 1
    data_normalized, norm_factor, max_norm = normalize_data(data)

    Nx_full, Ny_full, Nz_full, Nt = data.shape

    # Build the Resolution Pyramid
    if pyramid_levels == [1]:
        pyramid_levels = [1]
    num_levels = len(pyramid_levels)
    data_pyramid = []
    scale_factors = pyramid_levels

    for factor in scale_factors:
        data_pyramid.append(downsample_data_3d(data, factor))

    # Define kernel and step sizes
    max_kernel_size = np.array([9, 9, 9])
    max_step_size = np.array([3, 3, 3])
    min_kernel_size = np.array([3, 3, 3])
    min_step_size = np.array([1, 1, 1])    

    kernel_sizes = np.zeros((num_levels, 3), dtype=int)
    step_sizes = np.zeros((num_levels, 3), dtype=int)

    for level in range(num_levels):
        scale_ratio = level / (num_levels - 1) if num_levels > 1 else 0
        kernel_sizes[level] = np.round(min_kernel_size + scale_ratio * (max_kernel_size - min_kernel_size)).astype(int)
        step_sizes[level] = np.round(min_step_size + scale_ratio * (max_step_size - min_step_size)).astype(int)
        # Ensure kernel sizes are odd
        kernel_sizes[level] += 1 - (kernel_sizes[level] % 2)

    if num_levels == 1:
        kernel_sizes = np.array(kernel)
        step_sizes = np.array(step)

    # Initialize motion fields
    U_total = np.zeros((Nx_full, Ny_full, Nz_full, Nt), dtype=data.dtype)
    V_total = np.zeros((Nx_full, Ny_full, Nz_full, Nt), dtype=data.dtype)
    W_total = np.zeros((Nx_full, Ny_full, Nz_full, Nt), dtype=data.dtype)

    COST = np.zeros((max(iterations_per_level), num_levels))

    # Process each level in the pyramid
    for level in range(num_levels):
        current_data = data_pyramid[level]
        current_factor = scale_factors[level]

        scaled_kernel = kernel_sizes[level]
        scaled_step = step_sizes[level]

        Nx_curr, Ny_curr, Nz_curr, _ = current_data.shape
        scaled_kernel = np.minimum(scaled_kernel, [Nx_curr, Ny_curr, Nz_curr])
        k = (scaled_kernel - 1) / 2

        blocksize = list(scaled_kernel) + [adjusted_window_size]
        stepsize = list(scaled_step) + [step_time]
        #import pdb; pdb.set_trace()

        # Pad current data
        pad_width = [(int(k[0]), int(k[0])), (int(k[1]), int(k[1])), (int(k[2]), int(k[2])), (0,0)]
        data_padded = np.pad(current_data, pad_width, mode='edge')

        if level == 0:
            # Initialize motion fields at current level
            U_level = np.zeros_like(data_padded)
            V_level = np.zeros_like(data_padded)
            W_level = np.zeros_like(data_padded)
        else:
            # Use motion fields from previous level
            U_level = np.pad(U_level_upsampled, [(int(k[0]), int(k[0])), (int(k[1]), int(k[1])), 
                                                  (int(k[2]), int(k[2])), (0,0)], mode='edge')
            V_level = np.pad(V_level_upsampled, [(int(k[0]), int(k[0])), (int(k[1]), int(k[1])), 
                                                  (int(k[2]), int(k[2])), (0,0)], mode='edge')
            W_level = np.pad(W_level_upsampled, [(int(k[0]), int(k[0])), (int(k[1]), int(k[1])), 
                                                  (int(k[2]), int(k[2])), (0,0)], mode='edge')

        num_iterations = iterations_per_level[level]

        # Warp data using accumulated motion fields
        warped_data = warp_data_3d(data_padded, U_level, V_level, W_level)

        for iter_num in range(num_iterations):
            print(f'Level {level+1}/{num_levels}, Iteration {iter_num+1}/{num_iterations}')

            # Extract patches and corresponding indices
            data_patches, index_patches = im2col4(warped_data, blocksize, stepsize)
            import pdb; pdb.set_trace()

            Npatch = data_patches.shape[1]

            # Initialize arrays to store motion fields
            u_patches = Parallel(n_jobs=-1)(delayed(process_patch)(
                data_patches[:, nn], brightness_order) for nn in range(Npatch))
            
            # Unpack patch results
            du_patches, dv_patches, dw_patches = zip(*u_patches)

            # Combine incremental motion fields from patches
            dU = np.zeros_like(warped_data)
            dV = np.zeros_like(warped_data)
            dW = np.zeros_like(warped_data)
            counts_total = np.zeros_like(warped_data)

            for nn in range(Npatch):
                du_patch = du_patches[nn]
                dv_patch = dv_patches[nn]
                dw_patch = dw_patches[nn]
                indices = index_patches[:, nn]

                dU.flat[indices] += du_patch
                dV.flat[indices] += dv_patch
                dW.flat[indices] += dw_patch
                counts_total.flat[indices] += 1

            # Average motion fields where overlaps occur
            dU /= counts_total
            dV /= counts_total
            dW /= counts_total

            # Handle NaNs resulting from division by zero
            dU = np.nan_to_num(dU)
            dV = np.nan_to_num(dV)
            dW = np.nan_to_num(dW)

            # Update the motion fields
            U_level += dU
            V_level += dV
            W_level += dW

            # Warp data using accumulated motion fields
            warped_data = warp_data_3d(data_padded, U_level, V_level, W_level)

            # Compute COST (example metric; adjust as needed)
            COST[iter_num, level] = np.linalg.norm(warped_data - data_padded[..., 0])
        
        # Plot COST (optional)
        colors = ['r', 'b', 'g', 'm']
        plt.figure(100)
        plt.subplot(1, num_levels, level+1)
        plt.plot(COST[:num_iterations, level] / COST[0, level], color=colors[level % len(colors)], linewidth=1)
        plt.xlabel('Iteration')
        plt.ylabel('Normalized COST')
        plt.title(f'Level {level+1}')
        plt.grid(True)
        plt.show()

        # Remove padding
        U_level = unpad(U_level, scaled_kernel)
        V_level = unpad(V_level, scaled_kernel)
        W_level = unpad(W_level, scaled_kernel)

        # Upsample motion fields to the next finer level or to full resolution
        if level < num_levels - 1:
            next_data = data_pyramid[level + 1]
            target_size = next_data.shape[:3]
            U_level_upsampled = upsample_data_3d(U_level, target_size)
            V_level_upsampled = upsample_data_3d(V_level, target_size)
            W_level_upsampled = upsample_data_3d(W_level, target_size)

            # Calculate scaling ratios
            scale_factors_ratio = np.array(target_size) / np.array(data_pyramid[level].shape[:3])
            U_level_upsampled *= scale_factors_ratio[0]
            V_level_upsampled *= scale_factors_ratio[1]
            W_level_upsampled *= scale_factors_ratio[2]
        else:
            # At the finest level, adjust to full resolution if necessary
            if current_factor != 1 or num_levels == 1:
                target_size = data.shape[:3]
                U_level_upsampled = upsample_data_3d(U_level, target_size)
                V_level_upsampled = upsample_data_3d(V_level, target_size)
                W_level_upsampled = upsample_data_3d(W_level, target_size)

                # Calculate scaling ratios
                scale_factors_ratio = np.array(target_size) / np.array(data_pyramid[level].shape[:3])
                U_level_upsampled *= scale_factors_ratio[0]
                V_level_upsampled *= scale_factors_ratio[1]
                W_level_upsampled *= scale_factors_ratio[2]

    # Return the final motion fields
    U = U_level_upsampled
    V = V_level_upsampled
    W = W_level_upsampled
    WarpedData = warp_data_3d(data, U, V, W)

    # De-normalize Data
    WarpedData = denormalize_data(WarpedData, norm_factor, max_norm)

    return WarpedData, U, V, W, COST



def main():
    from ants import image_read
    import os
    import time

    root = "/Users/benaron/Documents/datasets"
    raw = os.path.join(root, "nii_raw_complex")
    dwi1 = os.path.join(raw, "sub-030_tp_1_MS_FWF_RMR_2iso_TE60_LTE_dwi.nii")
    nii = image_read(dwi1)
    img = nii.numpy()

    start_time = time.time()
    data_normalized,_,_,_,cost = offast_lkv_3d_4d_patchavg_pyramid_v3(img)
    end_time = time.time()
    elapsed_time = end_time - start_time
    hours, rem = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(rem, 60)
    print(f"test1 completed in {int(hours):02} hours, {int(minutes):02} minutes, {float(seconds):05} seconds.")


main()



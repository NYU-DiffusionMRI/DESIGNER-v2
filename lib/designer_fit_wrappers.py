import os
import warnings

import nibabel as nib
import numpy as np

import lib.io as mio


# =============================================================================
# Shared utilities
# =============================================================================

def _get_outlier_indices(outlier_locations, mask=None):
    """Return coordinates of voxels to process as a (3, N) integer array."""
    selection = outlier_locations if mask is None else mask
    return np.asarray(np.where(selection == 1))


def _make_csf_mask(md, threshold=1.2):
    """Create a boolean mask for high-MD/CSF-like voxels."""
    return np.asarray(md >= threshold, dtype=bool)


def _safe_log(values, eps=1e-12):
    """Log-transform signal values while preventing log(0) or log(negative)."""
    return np.log(np.clip(values, eps, None))


def _expand_patch_bounds(x, y, z, kernel, outlier_locations):
    """
    Grow a local patch until it contains at least one non-outlier voxel,
    or until the complete image extent has been reached.

    Returns
    -------
    xmin, xmax, ymin, ymax, zmin, zmax : int
        Slice bounds for the selected patch.
    """
    k = max(int(kernel) // 2, 1)
    shape = outlier_locations.shape

    while True:
        xmin = max(x - k - 1, 0)
        xmax = min(x + k, shape[0])
        ymin = max(y - k - 1, 0)
        ymax = min(y + k, shape[1])
        zmin = max(z - k - 1, 0)
        zmax = min(z + k, shape[2])

        outlier_patch = outlier_locations[xmin:xmax, ymin:ymax, zmin:zmax]

        # Patch contains at least one valid/non-outlier voxel.
        if not np.all(outlier_patch):
            return xmin, xmax, ymin, ymax, zmin, zmax

        # Entire volume has been reached; do not grow indefinitely.
        reached_image_bounds = (
            xmin == 0
            and xmax == shape[0]
            and ymin == 0
            and ymax == shape[1]
            and zmin == 0
            and zmax == shape[2]
        )

        if reached_image_bounds:
            return xmin, xmax, ymin, ymax, zmin, zmax

        k += 2

def _fixed_patch_bounds(x, y, z, kernel, outlier_locations):
    """
    Return a fixed-size patch centered on (x, y, z).

    The patch does not expand when no valid neighbors exist.
    """
    k = int(kernel) // 2
    shape = outlier_locations.shape

    xmin = max(x - k, 0)
    xmax = min(x + k + 1, shape[0])

    ymin = max(y - k, 0)
    ymax = min(y + k + 1, shape[1])

    zmin = max(z - k, 0)
    zmax = min(z + k + 1, shape[2])

    return xmin, xmax, ymin, ymax, zmin, zmax

def _get_low_b_shells(bval):
    """
    Identify the b=0 shell and the highest nominal shell below b=1.6.

    Nominal shells are defined at increments of 0.1
    (e.g., 0.9, 1.0, 1.1, 1.2).

    Values within ±0.1 of a nominal shell are grouped together.
    Assumes b=1 corresponds to 1000 s/mm².
    """
    bval = np.asarray(bval)

    # Identify b=0 shell.
    b_low_shell = bval <= 0.1

    # Find highest nominal nonzero shell below 1.6.
    low_bvals = bval[(bval > 0.1) & (bval < 1.6)]

    if low_bvals.size == 0:
        raise ValueError("No nonzero shells below b=1.6.")

    # Assign values to nominal shells (increments of 0.1).
    nominal_shells = np.round(low_bvals / 0.1) * 0.1

    # Highest nominal shell.
    b_high = np.max(nominal_shells)

    # Select measurements within ±0.1 of that nominal shell.
    b_high_shell = np.abs(bval - b_high) <= 0.1

    return (
        bval,
        b_high,
        0.0,
        b_high_shell,
        b_low_shell,
    )

def _neighbor_selection(
    x,
    y,
    z,
    kernel,
    outlier_locations,
    dwi_norm,
    bval,
    fa,
    md,
    md_mask,
    include_md_distance=True,
    smoothlevel=None,
):
    """
    Select similar, valid neighboring voxels around one outlier voxel.

    Returns
    -------
    bounds : tuple
        xmin, xmax, ymin, ymax, zmin, zmax
    goodidx : ndarray of bool
        Boolean vector selecting usable neighboring voxels.
    """
    (
        xmin,
        xmax,
        ymin,
        ymax,
        zmin,
        zmax,
    ) = _fixed_patch_bounds(x, y, z, kernel, outlier_locations)

    outlier_patch = outlier_locations[xmin:xmax, ymin:ymax, zmin:zmax]
    psize = outlier_patch.size

    # Exclude b=0 volumes from neighbor signal similarity.
    non_b0 = np.asarray(bval) >= 0.01

    if not np.any(non_b0):
        raise ValueError("No diffusion-weighted volumes available for neighbor selection.")

    # Extract only b > 0 volumes.
    patch = dwi_norm[
        xmin:xmax, ymin:ymax, zmin:zmax, :
    ][..., non_b0].reshape(psize, -1)

    ref = dwi_norm[x, y, z, non_b0][None, :]

    # Calculate signal similarity using only diffusion-weighted volumes.
    dwi_diff = patch - ref
    dwi_distance = np.sum(dwi_diff * dwi_diff, axis=1) / patch.shape[1]

    fapatch = fa[xmin:xmax, ymin:ymax, zmin:zmax].ravel()
    mdpatch = md[xmin:xmax, ymin:ymax, zmin:zmax].ravel()
    csfpatch = md_mask[xmin:xmax, ymin:ymax, zmin:zmax].ravel()

    if include_md_distance:
        md_std = mdpatch.std()
        md_distance = (mdpatch - md[x, y, z]) / (md_std + 1e-12)
        intensities = np.sqrt(dwi_distance + md_distance * md_distance)
    else:
        intensities = np.sqrt(dwi_distance)

    fa_std = fapatch.std()
    md_std = mdpatch.std()

    # Omit voxels with fa, md more than one stdev away from target outlier voxel (similar tissue type)
    omit = (
        (fapatch > fa[x, y, z] + fa_std)
        | (fapatch < fa[x, y, z] - fa_std)
        | (mdpatch > md[x, y, z] + md_std)
        | (mdpatch < md[x, y, z] - md_std)
    )

    invalid = outlier_patch.ravel().astype(bool) | csfpatch | omit

    valid_intensities = intensities.copy()
    max_distance = valid_intensities.max() if valid_intensities.size else 0.0
    valid_intensities[invalid] = max_distance

    if smoothlevel is None or smoothlevel is False:
        threshold = valid_intensities.mean()
        goodidx = (valid_intensities <= threshold) & ~invalid
    else:
        threshold = np.percentile(valid_intensities, smoothlevel)
        goodidx = (valid_intensities <= threshold) & ~invalid

    # If thresholding excludes every candidate, retain all valid candidates.
    if not np.any(goodidx):
        goodidx = ~invalid

    return (xmin, xmax, ymin, ymax, zmin, zmax), goodidx


# =============================================================================
# DWI smoothing / refitting
# =============================================================================

def parallel_outlier_smooth(
    inds,
    kernel,
    outlier_locations,
    dwi_norm,
    dwi,
    smoothlevel=None,
):
    """
    Estimate a replacement DWI signal using weighted similar neighboring voxels.
    """
    x, y, z = map(int, inds)

    (
        xmin,
        xmax,
        ymin,
        ymax,
        zmin,
        zmax,
    ) = _expand_patch_bounds(x, y, z, kernel, outlier_locations)

    outlier_patch = outlier_locations[xmin:xmax, ymin:ymax, zmin:zmax].ravel()
    psize = outlier_patch.size

    patch_norm = dwi_norm[
        xmin:xmax,
        ymin:ymax,
        zmin:zmax,
        :,
    ].reshape(psize, -1)

    patch_orig = dwi[
        xmin:xmax,
        ymin:ymax,
        zmin:zmax,
        :,
    ].reshape(psize, -1)

    ref = dwi_norm[x, y, z, :][None, :]
    distances = np.sqrt(np.sum((patch_norm - ref) ** 2, axis=1))
    distances /= max(patch_norm.shape[1], 1)

    valid = ~outlier_patch.astype(bool)

    if not np.any(valid):
        return dwi[x, y, z, :].copy()

    valid_distances = distances[valid]

    if smoothlevel is None or smoothlevel is False:
        threshold = np.median(valid_distances)
    else:
        threshold = np.percentile(valid_distances, smoothlevel)

    selected = valid & (distances <= threshold)

    if not np.any(selected):
        selected = valid

    selected_distances = distances[selected]
    max_distance = selected_distances.max()

    # More similar voxels receive larger weights.
    weights = max_distance - selected_distances

    # Equal distances would otherwise result in an all-zero weight vector.
    if np.allclose(weights.sum(), 0):
        weights = np.ones_like(weights)

    weights /= weights.sum()

    return np.sum(patch_orig[selected, :] * weights[:, None], axis=0)


def refit_or_smooth(
    outlier_locations,
    dwi,
    mask=None,
    smoothlevel=None,
    n_cores=-3,
):
    """Replace outlier DWI voxels with weighted neighboring estimates."""
    from joblib import Parallel, delayed

    outinds = _get_outlier_indices(outlier_locations, mask)

    if outinds.shape[1] == 0:
        return dwi.copy()

    dwi_norm = np.abs(dwi) / (np.amax(np.abs(dwi), axis=(0, 1, 2)) + 1e-12)
    dwi_new = dwi.copy()
    kernel = 7

    values = Parallel(n_jobs=n_cores, prefer="processes")(
        delayed(parallel_outlier_smooth)(
            outinds[:, i],
            kernel,
            outlier_locations,
            dwi_norm,
            dwi,
            smoothlevel,
        )
        for i in range(outinds.shape[1])
    )

    dwi_new[outinds[0], outinds[1], outinds[2], :] = np.asarray(values)
    return dwi_new

def replace_outlier_tensors(dki, dwi, dt, mask, outliers, L=None):

    mask = np.asarray(mask, dtype=bool)
    outliers = np.asarray(outliers, dtype=bool)

    print("\n--- replace_outlier_tensors ---")
    print("DWI shape:", dwi.shape)
    print("Mask shape:", mask.shape)
    print("Outliers shape:", outliers.shape)
    print("DT shape:", dt.shape)
    print("Number of outliers:", np.sum(outliers))

    if not np.any(outliers):
        print("No outliers detected!")
        return dwi, dt

    # Correct outlier DWI signals
    dwi_new = refit_or_smooth(
        outliers, dwi, n_cores=dki.n_cores
    )

    print("DWI changed:", np.any(dwi_new != dwi))

    # Fit corrected DWI at outlier locations
    replacement, _, _ = dki.dki_fit(
        dwi_new, outliers, L=L
    )

    outlier_idx = outliers[mask]

    print("Replacement shape:", replacement.shape)
    print("Outlier index shape:", outlier_idx.shape)
    print("Number of replacement voxels:", np.sum(outlier_idx))

    # Compare coefficients before replacement
    dt_before = dt[:, outlier_idx].copy()

    print("Replacement differs from original:",
          np.any(~np.isclose(replacement, dt_before, equal_nan=True)))

    # Replace coefficients
    dt[:, outlier_idx] = replacement

    print("DT changed:",
          np.any(~np.isclose(dt[:, outlier_idx], dt_before, equal_nan=True)))

    return dwi_new, dt
  
# =============================================================================
# AKC outlier mask processing
# =============================================================================

def parallel_outlier_akc(inds, akc_mask_tmp, akc_dirs, rk, md, fa):
    """
    Update an AKC mask value based on radial kurtosis and mean diffusivity.

    `akc_dirs` and `fa` are retained in the signature for compatibility with
    the original calling interface.
    """
    x, y, z = map(int, inds)

    output = akc_mask_tmp[x, y, z]

    if rk[x, y, z] < 0.5 and md[x, y, z] < 2:
        output = 1

    return int(output)


def akc_out(outlier_inds, akc_mask_tmp, akc_dirs, rk, md, fa, n_cores=-3):
    """
    Update the AKC outlier mask using RK and MD thresholds.

    Vectorized implementation:
    - No multiprocessing or pickling
    - Preserves the original mask
    - Updates only voxels specified by outlier_inds
    """

    new_akc_mask = akc_mask_tmp.copy()

    spatial_shape = new_akc_mask.shape

    if rk.shape != spatial_shape:
        raise ValueError(
            f"RK shape {rk.shape} does not match "
            f"AKC mask shape {spatial_shape}")

    if md.shape != spatial_shape:
        raise ValueError(
            f"MD shape {md.shape} does not match "
            f"AKC mask shape {spatial_shape}")

    if fa.shape != spatial_shape:
        raise ValueError(
            f"FA shape {fa.shape} does not match "
            f"AKC mask shape {spatial_shape}")

    if outlier_inds.shape[1] == 0:
        return new_akc_mask

    if outlier_inds.shape[0] != len(spatial_shape):
        raise ValueError(
            f"Expected {len(spatial_shape)} coordinate dimensions, "
            f"got {outlier_inds.shape[0]}"
        )

    coords = tuple(outlier_inds)

    additional_outliers = (
        (rk[coords] < 0.5) &
        (md[coords] < 2)
    )

    selected_coords = tuple(
        coord[additional_outliers]
        for coord in coords
    )

    new_akc_mask[selected_coords] = 1

    return new_akc_mask

# =============================================================================
# Slope-based b0 restoration
# =============================================================================

def parallel_outlier_slope(
    inds,
    kernel,
    outlier_locations,
    bval,
    dwi_norm,
    dwi,
    fa,
    md,
    md_mask,
    smoothlevel,
):
    """
    Estimate b=0 signals using diffusion-signal slopes from similar neighbors.

    Returns
    -------
    restored_b0 : ndarray
        Restored signal for each detected b=0 volume.
    original_slope : float
        Mean original outlier-voxel slope.
    restored_slope : float
        Mean selected neighboring slope.
    """
    warnings.filterwarnings(
        "ignore",
        message="Mean of empty slice",
        category=RuntimeWarning,
    )

    x, y, z = map(int, inds)

    (
        _,
        b_high,
        b_low,
        b_high_shell,
        b_low_shell,
    ) = _get_low_b_shells(bval)

    (
        xmin,
        xmax,
        ymin,
        ymax,
        zmin,
        zmax,
    ), goodidx = _neighbor_selection(
        x,
        y,
        z,
        kernel,
        outlier_locations,
        dwi_norm,
        bval,
        fa,
        md,
        md_mask,
        include_md_distance=True,
        smoothlevel=smoothlevel,
    )

    b_high_signal = np.mean(_safe_log(dwi[x, y, z, b_high_shell]))
    b_low_signal = _safe_log(dwi[x, y, z, b_low_shell])

    original_slope = np.abs(
        (b_high_signal - b_low_signal) / (b_high - b_low)
    )

    patch_high = _safe_log(
        dwi[xmin:xmax, ymin:ymax, zmin:zmax, :][..., b_high_shell]
    ).reshape(-1, np.sum(b_high_shell))

    patch_low = _safe_log(
        dwi[xmin:xmax, ymin:ymax, zmin:zmax, :][..., b_low_shell]
    ).reshape(-1, np.sum(b_low_shell))

    neighbor_high = patch_high[goodidx].mean(axis=1)
    neighbor_low = patch_low[goodidx]

    if neighbor_high.size == 0:
        selected_slope = original_slope
    else:
        neighbor_slopes = np.abs(
            (neighbor_high[:, None] - neighbor_low) / (b_high - b_low)
        )

        slope_std = np.nanstd(neighbor_slopes)

        # Preserve the newer, less restrictive slope criterion:
        # include neighbors > or = the outlier slope.
        eligible = neighbor_slopes >= original_slope

        masked_slopes = np.where(eligible, neighbor_slopes, np.nan)
        selected_slope = np.nanmean(masked_slopes, axis=0)

        selected_slope = np.where(
            np.isnan(selected_slope),
            original_slope,
            selected_slope,
        )

    restored_log_b0 = np.abs(selected_slope) * b_high + b_high_signal
    restored_b0 = np.exp(restored_log_b0).squeeze()

    return (
        restored_b0,
        float(np.nanmean(original_slope)),
        float(np.nanmean(selected_slope)),
    )


def db0restore_slope(
    outlier_locations,
    dwi,
    bval,
    k,
    perc,
    fa,
    md,
    mask=None,
    n_cores=-3,
    return_slope_maps=True,
):
    """
    Restore b=0 volumes using neighboring diffusion-signal slopes.

    Parameters
    ----------
    return_slope_maps : bool, default=True
        If True, return the original and replacement slope maps in addition to
        the corrected DWI. 
    Returns
    -------
        dwi_new, original_slope_map, restored_slope_map

    """
    from joblib import Parallel, delayed
    from tqdm import tqdm

    outinds = _get_outlier_indices(outlier_locations, mask)
    dwi_new = dwi.copy()

    if outinds.shape[1] == 0:
        if return_slope_maps:
            return dwi_new, np.zeros_like(md, dtype=float), np.zeros_like(md, dtype=float)
        return dwi_new

    dwi_norm = np.abs(dwi) / (np.amax(np.abs(dwi), axis=(0, 1, 2)) + 1e-12)
    md_mask = _make_csf_mask(md)

    results = Parallel(n_jobs=n_cores, prefer="processes")(
        delayed(parallel_outlier_slope)(
            outinds[:, i],
            k,
            outlier_locations,
            bval,
            dwi_norm,
            dwi,
            fa,
            md,
            md_mask,
            perc,
        )
        for i in tqdm(range(outinds.shape[1]))
    )

    restored_values, original_slopes, restored_slopes = zip(*results)

    restored_values = np.asarray(restored_values)
    original_slopes = np.asarray(original_slopes)
    restored_slopes = np.asarray(restored_slopes)

    b0_indices = np.where(np.asarray(bval) < 0.01)[0]

    for b0_position, b0_index in enumerate(b0_indices):
        dwi_new[
            outinds[0],
            outinds[1],
            outinds[2],
            b0_index,
        ] = restored_values[:, b0_position]

    if not return_slope_maps:
        return dwi_new

    original_slope_map = np.zeros_like(md, dtype=float)
    restored_slope_map = np.zeros_like(md, dtype=float)

    original_slope_map[outinds[0], outinds[1], outinds[2]] = original_slopes
    restored_slope_map[outinds[0], outinds[1], outinds[2]] = restored_slopes

    return dwi_new, original_slope_map, restored_slope_map
  
  
# ============================================================================
# Output utilities
# =============================================================================

def save_params(paramDict, niiex, model, outdir, format="nifti"):
    """
    Save fitted parameter volumes as MIF or NIfTI files.

    Existing output files are not overwritten.
    """
    os.makedirs(outdir, exist_ok=True)

    for key, volume in paramDict.items():
        if "L" in key:
            outpath = os.path.join(outdir, f"{key}.mif")
        else:
            outpath = os.path.join(outdir, f"{key}_{model}.mif")

        if os.path.exists(outpath):
            continue

        mif = mio.Image(
            data=volume,
            vox=niiex.vox,
            transform=niiex.transform,
            grad=niiex.grad,
        )

        if format == "mif":
            mif.save(outpath)

        elif format == "nifti":
            # Get spatial voxel sizes.
            voxel_sizes = np.asarray(mif.vox, dtype=float).ravel()

            # NIfTI affine matrices are 4x4.
            # Ensure that we have three spatial voxel dimensions.
            spatial_voxels = np.ones(3, dtype=float)

            n_vox = min(len(voxel_sizes), 3)
            spatial_voxels[:n_vox] = voxel_sizes[:n_vox]

            # The fourth dimension is not a spatial dimension.
            voxel_sizes_4d = np.append(spatial_voxels, 1.0)

            transform = np.asarray(mif.transform, dtype=float)

            if transform.shape != (4, 4):
                raise ValueError(
                    f"Expected a 4x4 affine transform, got {transform.shape}"
                )

            affine = transform @ np.diag(voxel_sizes_4d)

            nii = nib.Nifti1Image(
                mif.data,
                affine,
            )

            nib.save(nii, outpath.replace(".mif", ".nii"))

        else:
            raise ValueError(
                f"Unsupported output format: {format!r}. "
                "Use 'mif' or 'nifti'."
            )
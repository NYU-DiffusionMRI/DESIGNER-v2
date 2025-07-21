from pathlib import Path
from typing import Union, Optional, Union, Dict

import numpy as np
import nibabel as nib
from nibabel.nifti1 import Nifti1Image

from tests.types import StatsDict, DWIStage, DiffusionModelType, is_valid_dwi_stage, is_valid_diffusion_model_type


def assert_stats(
    stats: StatsDict,
    expected_stats: StatsDict,
    *,
    tolerances: Dict,
    context: Union[DWIStage, DiffusionModelType]
):
    """
    Assert that two StatsDict objects are approximately equal for each ROI,
    using tolerances for mean and std from the provided dictionary.

    Args:
        stats: Actual statistics dict mapping ROI names to [mean, (std)].
        expected_stats: Expected statistics dict, same format as stats.
        tolerances: Dict of tolerance values for mean and std.
        context: DWIStage or DiffusionModelType, used for tolerance selection and logging.
    """
    print(f"context: {context}")

    # Determine context type and set up mean/std comparison functions
    if is_valid_dwi_stage(context):
        def _assert_mean(roi, actual, expected):
            default_key = f"b0_mean_rtol"
            key = f"b0_{roi}_mean_rtol"
            rtol = tolerances.get(key, tolerances[default_key])
            assert np.isclose(actual, expected, rtol=rtol)
        def _assert_std(roi, actual, expected):
            default_key = f"b0_std_rtol"
            key = f"b0_{roi}_std_rtol"
            rtol = tolerances.get(key, tolerances[default_key])
            assert np.isclose(actual, expected, rtol=rtol)
    elif is_valid_diffusion_model_type(context):
        def _assert_mean(roi, actual, expected):
            default_key = f"fa_mean_atol"
            key = f"fa_{roi}_mean_atol"
            atol = tolerances.get(key, tolerances[default_key])
            assert np.isclose(actual, expected, atol=atol)
        def _assert_std(roi, actual, expected):
            default_key = f"fa_std_atol"
            key = f"fa_{roi}_std_atol"
            atol = tolerances.get(key, tolerances[default_key])
            assert np.isclose(actual, expected, atol=atol)
    else:
        raise ValueError(f"Invalid context: {context}")

    for roi, stats_values in stats.items():
        print(f"roi: {roi}")
        print(f"actual mean: {stats_values[0]}, expected mean: {expected_stats[roi][0]}")
        _assert_mean(roi, stats_values[0], expected_stats[roi][0])

        if len(stats_values) > 1 and len(expected_stats[roi]) > 1:   # standard deviation exists
            print(f"actual std: {stats_values[1]}, expected std: {expected_stats[roi][1]}")
            _assert_std(roi, stats_values[1], expected_stats[roi][1])



def create_binary_mask_from_fa(
    fa: Union[Path, Nifti1Image],
    brain_mask: Optional[Union[Path, Nifti1Image]] = None,
    threshold: float = 0.3
) -> Nifti1Image:
    """
    Generate a binary mask from a fractional anisotropy (FA) image, optionally restricted by a brain mask.

    Args:
        fa: Path or Nifti1Image of the FA image.
        threshold: Voxels with FA >= threshold are set to 1 (default: 0.3).
        brain_mask: Optional path or Nifti1Image of a brain mask. If provided, restricts the mask to this region.

    Returns:
        Nifti1Image: Binary mask (uint8) in the same space as the canonicalized FA image.
    """
    img = Nifti1Image.from_filename(fa) if isinstance(fa, Path) else fa
    fa_can = nib.as_closest_canonical(img)
    data = fa_can.get_fdata()

    # optionally load & apply skull-strip mask
    if brain_mask is not None:
        mask_img = Nifti1Image.from_filename(brain_mask) if isinstance(brain_mask, Path) else brain_mask
        mask_can = nib.as_closest_canonical(mask_img)
        data *= mask_can.get_fdata()

    data_no_nan = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)
    binary = (data_no_nan >= threshold).astype(np.uint8)

    return Nifti1Image(binary, fa_can.affine, fa_can.header)


def extract_mean_b0(img_path: Path, bval_path: Path) -> Nifti1Image:
    """Extract mean b0 volumes ensuring correct diffusion axis handling.
    
    Args:
        img_path: Path to the 4D diffusion image
        bval_path: Path to the corresponding bvals file
        
    Returns:
        Nifti1Image containing the mean b0 image with preserved header and affine
        
    Raises:
        ValueError: If no b=0 volumes found or if image dimensions don't match bvals
    """
    img = Nifti1Image.from_filename(img_path)
    data = img.get_fdata()
    
    img_shape = data.shape
    if len(img_shape) != 4:
        raise ValueError(f"Expected 4D image with diffusion volumes, but got shape {img_shape}")
        
    bvals = np.loadtxt(bval_path)
    n_volumes = img_shape[-1]  # Last dimension should be diffusion
    if len(bvals) != n_volumes:
        raise ValueError(f"Number of bvals ({len(bvals)}) doesn't match number of volumes ({n_volumes})")

    b0_indices = np.where(bvals == 0)[0]
    if b0_indices.size == 0:
        raise ValueError("No b=0 volumes found in the data")

    b0_data = data[..., b0_indices]
    b0_data_no_nan = np.nan_to_num(b0_data, nan=0.0, posinf=0.0, neginf=0.0)
    
    mean_b0 = np.mean(b0_data_no_nan, axis=-1)

    mean_b0_img = Nifti1Image(mean_b0, img.affine, header=img.header)
    mean_b0_img.header.set_data_shape(mean_b0.shape)
    
    return mean_b0_img


def compute_roi_mean_and_std(image: Nifti1Image, roi: Nifti1Image) -> tuple[float, float]:
    """Compute the mean and standard deviation of data within a region of interest (ROI).

    Args:
        image: Input image (Nifti1Image)
        roi: ROI image (Nifti1Image)

    Returns:
        Tuple containing (mean, std) of the data values within the ROI

    Note:
        Only non-zero values in the masked data are used for computation.
        The ROI should be a binary mask where 1 indicates ROI voxels and 0 indicates background.
    """
    image_can = nib.as_closest_canonical(image)
    roi_can = nib.as_closest_canonical(roi)
    data_masked = image_can.get_fdata() * roi_can.get_fdata()

    m = np.mean(data_masked[data_masked != 0])
    s = np.std(data_masked[data_masked != 0])

    return float(m), float(s)

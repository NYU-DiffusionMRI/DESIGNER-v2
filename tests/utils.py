from pathlib import Path
from typing import Union, List

import numpy as np
import nibabel as nib
from nibabel.nifti1 import Nifti1Image

from tests.types import StatsDict


def assert_roi_mean_and_std(image: Nifti1Image, roi: Union[Path, Nifti1Image], expected_values: List[float]):
    """Assert the mean and standard deviation of the data within the ROI are close to the expected values.

    Args:
        image: Input image
        roi: ROI specification - can be file path (Path) or Nifti1Image
        expected_values: Expected [mean, std] or [mean] values
    """
    if isinstance(roi, Path):
        roi_image = nib.load(roi)
    elif isinstance(roi, Nifti1Image):
        roi_image = roi
    else:
        raise ValueError(f"Invalid ROI type: {type(roi)}")
    
    mean, std = compute_roi_mean_and_std(image, roi_image)
    # print(f"mean: {mean}, std: {std}")

    if len(expected_values) == 2:
        assert np.isclose(mean, expected_values[0])
        assert np.isclose(std, expected_values[1])
    elif len(expected_values) == 1:
        assert np.isclose(mean, expected_values[0])
    else:
        raise ValueError(f"Expected values must be a list of length 1 or 2, got {len(expected_values)}")


def assert_stats(stats: StatsDict, expected_stats: StatsDict):
    for roi, stats_values in stats.items():
        for val, expected_val in zip(stats_values, expected_stats[roi]):
            assert np.isclose(val, expected_val)


def create_binary_mask_from_fa(fa_file: Path, threshold: float = 0.3) -> Nifti1Image:
    """Create a binary mask from a fractional anisotropy (FA) image.

    Args:
        fa_file: Path to the FA NIfTI file
        threshold: FA threshold value to create binary mask (voxels >= threshold become 1, default: 0.3)

    Returns:
        Binary mask as Nifti1Image
    """
    img = nib.load(str(fa_file))
    data = img.get_fdata()

    data_no_nan = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)
    binary_data = (data_no_nan >= threshold).astype(np.uint8)

    return Nifti1Image(binary_data, img.affine, img.header)


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
    img = nib.load(img_path)
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

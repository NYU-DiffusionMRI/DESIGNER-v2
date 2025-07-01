from pathlib import Path
from typing import Optional, Union, List

import numpy as np
import nibabel as nib


def assert_roi_mean_and_std(data: np.ndarray, roi: Union[Path, np.ndarray], expected_values: List[float]):
    """Assert the mean and standard deviation of the data within the ROI are close to the expected values.

    Args:
        data: Input data array (e.g., image data) with shape (X, Y, Z)
        roi: ROI specification - can be file path (Path) or numpy array
        expected_values: Expected [mean, std] or [mean] values
    """
    if isinstance(roi, np.ndarray):
        roi_data = roi
    elif isinstance(roi, Path):
        roi_data = nib.load(roi).get_fdata()
    else:
        raise ValueError(f"Invalid ROI type: {type(roi)}")

    mean, std = compute_roi_mean_and_std(data, roi_data)

    if len(expected_values) == 2:
        assert np.isclose(mean, expected_values[0])
        assert np.isclose(std, expected_values[1])
    elif len(expected_values) == 1:
        assert np.isclose(mean, expected_values[0])
    else:
        raise ValueError(f"Expected values must be a list of length 1 or 2, got {len(expected_values)}")


def create_binary_mask_from_fa(fa_file: Path, output_file: Optional[Path] = None, *, threshold: float) -> np.ndarray:
    """Create a binary mask from a fractional anisotropy (FA) image.

    Args:
        fa_file: Path to the FA NIfTI file
        output_file: Optional path to save the binary mask as NIfTI file
        threshold: FA threshold value to create binary mask (voxels >= threshold become 1)

    Returns:
        Binary mask as numpy array (uint8, 0s and 1s)

    Note:
        NaN values in the FA image are converted to 0 before thresholding
    """
    img = nib.load(str(fa_file))
    data = img.get_fdata()

    data_no_nan = np.nan_to_num(data, nan=0.0)
    binary_data = (data_no_nan >= threshold).astype(np.uint8)

    if output_file is not None:
        new_img = nib.Nifti1Image(binary_data, img.affine, img.header)
        nib.save(new_img, str(output_file))

    return binary_data


def extract_mean_b0(img_path: Path, bval_path: Path) -> np.ndarray:
    """Extract the mean b=0 volume from a diffusion-weighted image.

    Args:
        img_path: Path to the DWI NIfTI file
        bval_path: Path to the b-value text file

    Returns:
        Mean b=0 volume as numpy array with shape (X, Y, Z)

    Raises:
        ValueError: If no b=0 volumes are found in the data
    """
    img = nib.load(img_path)
    data = img.get_fdata()  # shape: (X, Y, Z, N)

    bvals = np.loadtxt(bval_path)

    b0_indices = np.where(bvals == 0)[0]
    if b0_indices.size == 0:
        raise ValueError("No b=0 volumes found in the data")

    b0_data = data[..., b0_indices]
    b0_data_no_nan = np.nan_to_num(b0_data, nan=0.0)
    mean_b0 = np.mean(b0_data_no_nan, axis=3)

    return mean_b0


def compute_roi_mean_and_std(data: np.ndarray, roi: np.ndarray) -> tuple[float, float]:
    """Compute the mean and standard deviation of data within a region of interest (ROI).

    Args:
        data: Input data array (e.g., image data) with shape (X, Y, Z)
        roi: Binary ROI mask array with same spatial dimensions as data (0s and 1s)

    Returns:
        Tuple containing (mean, std) of the data values within the ROI

    Note:
        Only non-zero values in the masked data are used for computation.
        The ROI should be a binary mask where 1 indicates ROI voxels and 0 indicates background.
    """
    data_masked = data * roi

    m = np.mean(data_masked[data_masked != 0])
    s = np.std(data_masked[data_masked != 0])

    return float(m), float(s)

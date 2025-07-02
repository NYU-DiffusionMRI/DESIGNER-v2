import subprocess
import shutil
import json
from pathlib import Path
from typing import List, Tuple, Dict

import numpy as np
import nibabel as nib
from nibabel.nifti1 import Nifti1Image
import pytest

from tests.utils import create_binary_mask_from_fa, extract_mean_b0, assert_roi_mean_and_std


ground_truth = json.load(open("tests/ground_truth_statistics/D2.json"))


@pytest.fixture(scope="module")
def paths():
    data_dir = Path("tests/data/D2/")
    tmp_dir = Path("tests/tmp_D2/")
    processing_dir = tmp_dir / "processing"
    params_dir = tmp_dir / "params"

    return {
        # processing and parameters directories for the pipeline
        "tmp_dir": tmp_dir,
        "processing_dir": processing_dir,
        "params_dir": params_dir,
        # images to be processed by pipeline
        "dwi_images": [data_dir / "meso_ds.nii.gz", data_dir / "research_ds.nii.gz"],
        # Phase Encoding image
        "phase_encoding_image": data_dir / "pa_ds.nii.gz",
        # Designer output image
        "dwi_designer": tmp_dir / "dwi_designer.nii",
        "dwi_designer_bvec": tmp_dir / "dwi_designer.bvec",
        "dwi_designer_bval": tmp_dir / "dwi_designer.bval",
        # FA images
        "fa_dki": params_dir / "fa_dki.nii",
        "fa_dti": params_dir / "fa_dti.nii",
        "fa_wdki": params_dir / "fa_wdki.nii",
        # ROI images
        "roi1": data_dir / "roi1.nii.gz",
        "roi2": data_dir / "roi2.nii.gz",
        "voxel": data_dir / "voxel.nii.gz",
        # Intermediate image (after topup)
        "dwi_topup": processing_dir / "eddy_processing" / "dwi_pe_0_applytopup.nii.gz",
        "dwi_topup_bval": tmp_dir / "dwi_designer.bval",    # same as designer bval
    }


@pytest.fixture(scope="module")
def ground_truth_data():
    return ground_truth


@pytest.fixture(scope="module", autouse=True)
def run_pipeline(paths):
    dwi_args=",".join([str(path) for path in paths["dwi_images"]])

    processing_dir=paths["processing_dir"]
    params_dir=paths["params_dir"]
    designer_image_path = paths["dwi_designer"]

    if processing_dir.exists():
        shutil.rmtree(processing_dir)

    try:
        # TODO: passing relative path for phase_encoding_image doesn't work. need to fix it in DESIGNER app code.
        ret = subprocess.run(["designer", "-set_seed", "-eddy", "-rpe_pair", str(paths["phase_encoding_image"].resolve()), "-mask", "-scratch", str(processing_dir), "-nocleanup", dwi_args, str(designer_image_path)])
        assert ret.returncode == 0
        assert designer_image_path.exists()

        if params_dir.exists():
            shutil.rmtree(params_dir)

        ret = subprocess.run(["tmi", "-SMI", "-DKI", "-WDKI", "-DTI", "-mask", str(processing_dir / "brain_mask.nii"), str(designer_image_path), str(params_dir)])
        assert ret.returncode == 0
        assert params_dir.exists()

        yield

    finally:
        if paths["tmp_dir"].exists():
            shutil.rmtree(paths["tmp_dir"])


@pytest.fixture(scope="module")
def white_matter_roi(paths) -> Nifti1Image:
    return create_binary_mask_from_fa(paths["fa_dki"])


def test_white_matter_voxel_count(white_matter_roi, ground_truth_data):
    wm_voxel_cnt = np.count_nonzero(white_matter_roi.get_fdata())
    expected_count = ground_truth_data["white_matter_voxel_count"]
    assert wm_voxel_cnt == expected_count


def get_b0_test_params(ground_truth) -> List[Tuple[str, str, Dict]]:
    """Generate test parameters for B0 stats each step from ground truth data."""
    b0_stats = ground_truth["b0_stats"]
    params = []
    for dwi_key, values in b0_stats.items():
        bval_key = f"{dwi_key}_bval"
        params.append((dwi_key, bval_key, values))
    return params

@pytest.mark.parametrize("dwi_key, bval_key, expected_values", get_b0_test_params(ground_truth))
def test_b0_stats(paths, white_matter_roi, dwi_key, bval_key, expected_values):
    b0_image = extract_mean_b0(paths[dwi_key], paths[bval_key])

    assert_roi_mean_and_std(b0_image, white_matter_roi, expected_values["wm"])
    assert_roi_mean_and_std(b0_image, paths["roi1"], expected_values["roi1"])
    assert_roi_mean_and_std(b0_image, paths["roi2"], expected_values["roi2"])
    assert_roi_mean_and_std(b0_image, paths["voxel"], [expected_values["voxel"]])


def get_fa_test_params(ground_truth) -> List[Tuple[str, Dict]]:
    """Generate test parameters for FA stats from ground truth data."""
    fa_stats = ground_truth["fa_stats"]
    return [(fa_type, values) for fa_type, values in fa_stats.items()]


@pytest.mark.parametrize("fa_type, expected_values", get_fa_test_params(ground_truth))
def test_fa_stats(paths, white_matter_roi, fa_type, expected_values):
    fa_image = nib.load(paths[fa_type])
    fa_data = fa_image.get_fdata()
    fa_data_no_nan = np.nan_to_num(fa_data, nan=0.0, posinf=0.0, neginf=0.0)
    fa_image_no_nan = Nifti1Image(fa_data_no_nan, fa_image.affine, fa_image.header)

    assert_roi_mean_and_std(fa_image_no_nan, white_matter_roi, expected_values["wm"])
    assert_roi_mean_and_std(fa_image_no_nan, paths["roi1"], expected_values["roi1"])
    assert_roi_mean_and_std(fa_image_no_nan, paths["roi2"], expected_values["roi2"])
    assert_roi_mean_and_std(fa_image_no_nan, paths["voxel"], [expected_values["voxel"]])

import subprocess
import shutil
import json
from pathlib import Path

import numpy as np
import nibabel as nib
import pytest

from tests.utils import create_binary_mask_from_fa, extract_mean_b0, assert_roi_mean_and_std


ground_truth = json.load(open("tests/ground_truth_statistics/D4.json"))


@pytest.fixture(scope="module")
def paths():
    data_dir = Path("tests/data/D4/")
    tmp_dir = Path("tests/tmp_D4/")
    processing_dir = tmp_dir / "processing"
    params_dir = tmp_dir / "params"

    return {
        # processing and parameters directories for the pipeline
        "tmp_dir": tmp_dir,
        "processing_dir": processing_dir,
        "params_dir": params_dir,
        # images to be processed by pipeline
        "dwi_images": [data_dir / "dif1_slice.nii.gz", data_dir / "dif2_slice.nii.gz"],
        "phase_images": [data_dir / "dif1phase_slice.nii.gz", data_dir / "dif2phase_slice.nii.gz"],
        # Designer output image
        "dwi_designer": tmp_dir / "dwi_designer.nii",
        "bvec": tmp_dir / "dwi_designer.bvec",
        "bval": tmp_dir / "dwi_designer.bval",
        # FA images
        "fa_dki": params_dir / "fa_dki.nii",
        "fa_dti": params_dir / "fa_dti.nii",
        "fa_wdki": params_dir / "fa_wdki.nii",
        # ROI images
        "roi1": data_dir / "roi1.nii.gz",
        "roi2": data_dir / "roi2.nii.gz",
        "voxel": data_dir / "voxel.nii.gz",
        # Intermediate image (after denoising)
        "dwidn": processing_dir / "dwidn.nii",
        "dwidn_bvec": processing_dir / "dwidn.bvec",
        "dwidn_bval": processing_dir / "dwidn.bval",
        # Intermediate image (after degibbs)
        "dwidg": processing_dir / "working_rpg.nii",
        "dwidg_bvec": processing_dir / "working_rpg.bvec",
        "dwidg_bval": processing_dir / "working_rpg.bval",
    }


@pytest.fixture(scope="module")
def ground_truth_data():
    return ground_truth


@pytest.fixture(scope="module", autouse=True)
def run_pipeline(paths):
    dwi_args=",".join([str(path) for path in paths["dwi_images"]])
    phase_args=",".join([str(path) for path in paths["phase_images"]])

    processing_dir=paths["processing_dir"]
    params_dir=paths["params_dir"]
    designer_image_path = paths["dwi_designer"]

    if processing_dir.exists():
        shutil.rmtree(processing_dir)

    try:
        ret = subprocess.run(["designer", "-denoise", "-phase", phase_args, "-degibbs", "-normalize", "-mask", "-scratch", str(processing_dir), "-nocleanup", dwi_args, str(designer_image_path)])
        assert ret.returncode == 0
        assert designer_image_path.exists()

        if params_dir.exists():
            shutil.rmtree(params_dir)

        ret = subprocess.run(["tmi", "-SMI", "-DKI", "-WDKI", "-DTI", "-sigma", str(processing_dir / "sigma.nii"), "-mask", str(processing_dir / "brain_mask.nii"), str(designer_image_path), str(params_dir)])
        assert ret.returncode == 0
        assert params_dir.exists()

        yield

    finally:
        if paths["tmp_dir"].exists():
            shutil.rmtree(paths["tmp_dir"])


@pytest.fixture(scope="module")
def white_matter_roi(paths):
    return create_binary_mask_from_fa(paths["fa_dki"], threshold=0.3)


def test_white_matter_voxel_count(white_matter_roi, ground_truth_data):
    wm_voxel_cnt = np.count_nonzero(white_matter_roi)
    expected_count = ground_truth_data["white_matter_voxel_count"]
    assert wm_voxel_cnt == expected_count


def test_b0_stats(paths, white_matter_roi, ground_truth_data):
    b0_data = extract_mean_b0(paths["dwi_designer"], paths["bval"])
    expected_values = ground_truth_data["b0_stats"]

    assert_roi_mean_and_std(b0_data, white_matter_roi, expected_values["wm"])
    assert_roi_mean_and_std(b0_data, paths["roi1"], expected_values["roi1"])
    assert_roi_mean_and_std(b0_data, paths["roi2"], expected_values["roi2"])
    assert_roi_mean_and_std(b0_data, paths["voxel"], [expected_values["voxel"]])


def get_fa_test_params(ground_truth):
    """Generate test parameters for FA stats from ground truth data."""
    fa_stats = ground_truth["fa_stats"]
    return [(fa_type, values) for fa_type, values in fa_stats.items()]


@pytest.mark.parametrize("fa_type, expected_values", get_fa_test_params(ground_truth))
def test_fa_stats(paths, white_matter_roi, fa_type, expected_values):
    fa_data = nib.load(paths[fa_type]).get_fdata()
    fa_data = np.nan_to_num(fa_data, nan=0.0)

    assert_roi_mean_and_std(fa_data, white_matter_roi, expected_values["wm"])
    assert_roi_mean_and_std(fa_data, paths["roi1"], expected_values["roi1"])
    assert_roi_mean_and_std(fa_data, paths["roi2"], expected_values["roi2"])
    assert_roi_mean_and_std(fa_data, paths["voxel"], [expected_values["voxel"]])


def get_b0_each_step_test_params(ground_truth):
    """Generate test parameters for B0 stats each step from ground truth data."""
    b0_stats = ground_truth["b0_stats_each_step"]
    params = []
    for dwi_key, values in b0_stats.items():
        bval_key = f"{dwi_key}_bval"
        params.append((dwi_key, bval_key, values))
    return params


@pytest.mark.parametrize("dwi_key, bval_key, expected_values", get_b0_each_step_test_params(ground_truth))
def test_b0_stats_each_step(paths, white_matter_roi, dwi_key, bval_key, expected_values):
    b0_data = extract_mean_b0(paths[dwi_key], paths[bval_key])

    assert_roi_mean_and_std(b0_data, white_matter_roi, expected_values["wm"])
    assert_roi_mean_and_std(b0_data, paths["roi1"], expected_values["roi1"])
    assert_roi_mean_and_std(b0_data, paths["roi2"], expected_values["roi2"])
    assert_roi_mean_and_std(b0_data, paths["voxel"], [expected_values["voxel"]])

import subprocess
import shutil
import json
from typing import List
from pathlib import Path
from functools import reduce

import numpy as np
import nibabel as nib
import pytest

from tests.utils import create_binary_mask_from_fa, extract_mean_b0, assert_roi_mean_and_std


ground_truth = json.load(open("tests/ground_truth_statistics/D5.json"))
te_list = [60.0, 92.0]

@pytest.fixture(scope="module")
def paths():
    data_dir = Path("tests/data/D5/")
    tmp_dir = Path("tests/tmp_D5/")
    processing_dir = tmp_dir / "processing"
    params_dir = tmp_dir / "params"

    return {
        # processing and parameters directories for the pipeline
        "tmp_dir": tmp_dir,
        "processing_dir": processing_dir,
        "params_dir": params_dir,
        # images to be processed by pipeline
        "dwi_images": [data_dir / "dwi1_ds.nii.gz", data_dir / "dwi2_ds.nii.gz", data_dir / "dwi3_ds.nii.gz"],
        "phase_encoding_image": data_dir / "rpe_ds.nii.gz",
        "phase_images": [data_dir / "phase1_ds.nii.gz", data_dir / "phase2_ds.nii.gz", data_dir / "phase3_ds.nii.gz"],
        # Designer output image
        "dwi_designer": tmp_dir / "dwi_designer.nii",
        "bvec": tmp_dir / "dwi_designer.bvec",
        "bval": tmp_dir / "dwi_designer.bval",
        # FA images
        "fa_dki": [params_dir / f"fa_dki_te{te}.nii" for te in te_list],
        "fa_dti": [params_dir / f"fa_dti_te{te}.nii" for te in te_list],
        "fa_wdki": [params_dir / f"fa_wdki_te{te}.nii" for te in te_list],
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
        # Intermediate image (after eddy/topup)
        "dwiec": processing_dir / "dwiec.nii",
        "dwiec_bvec": processing_dir / "dwiec.bvec",
        "dwiec_bval": processing_dir / "dwiec.bval",
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
        # TODO: passing relative path for phase_encoding_image doesn't work. need to fix it in DESIGNER app code.
        ret = subprocess.run(["designer", "-set_seed", "-denoise", "-shrinkage", "frob", "-degibbs", "-eddy", "-rpe_pair", str(paths["phase_encoding_image"].resolve()), "-normalize", "-mask", "-scratch", str(processing_dir), "-nocleanup", "-eddy_fakeb", "0.85,1.2,1", "-rpe_te", "60", "-bshape", "1,0.6,1", "-echo_time", "0.060,0.078,0.092", "-phase", phase_args, dwi_args, str(designer_image_path)])
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
def white_matter_roi_list(paths) -> List[np.ndarray]:
    return [create_binary_mask_from_fa(fa_file, threshold=0.3) for fa_file in paths["fa_dki"]]

@pytest.fixture(scope="module")
def white_matter_roi(white_matter_roi_list) -> np.ndarray:
    return reduce(np.bitwise_or, white_matter_roi_list)


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
    ret_params = []
    
    for fa_type, echo_time_stats in fa_stats.items():
        sorted_echo_times = sorted(echo_time_stats.keys())
        fa_stats_by_echo_time = [echo_time_stats[echo_time] for echo_time in sorted_echo_times]
        ret_params.append((fa_type, fa_stats_by_echo_time))

    return ret_params


@pytest.mark.parametrize("fa_type, expected_values_list", get_fa_test_params(ground_truth))
def test_fa_stats(paths, white_matter_roi_list, fa_type, expected_values_list):
    fa_file_list = paths[fa_type]
    fa_data_list = []
    for fa_file in fa_file_list:
        fa_data = nib.load(fa_file).get_fdata()
        no_nan_fa_data = np.nan_to_num(fa_data, nan=0.0)

        fa_data_list.append(no_nan_fa_data)

    for fa_data, white_matter_roi, expected_values in zip(fa_data_list, white_matter_roi_list, expected_values_list):
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

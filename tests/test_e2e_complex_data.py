import subprocess
import shutil
import json
from typing import List
from pathlib import Path
from typing import List, Tuple, Dict
from functools import reduce

import numpy as np
import nibabel as nib
from nibabel.nifti1 import Nifti1Image
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
        "dwi_designer_bvec": tmp_dir / "dwi_designer.bvec",
        "dwi_designer_bval": tmp_dir / "dwi_designer.bval",
        # FA images
        "fa_dki": [params_dir / f"fa_dki_te{te}.nii" for te in te_list],
        "fa_dti": [params_dir / f"fa_dti_te{te}.nii" for te in te_list],
        "fa_wdki": [params_dir / f"fa_wdki_te{te}.nii" for te in te_list],
        # ROI images
        "roi1": data_dir / "roi1.nii.gz",
        "roi2": data_dir / "roi2.nii.gz",
        "voxel": data_dir / "voxel.nii.gz",
        # Intermediate image (after denoising)
        "dwi_denoising": processing_dir / "dwidn.nii",
        "dwi_denoising_bvec": processing_dir / "dwidn.bvec",
        "dwi_denoising_bval": processing_dir / "dwidn.bval",
        # Intermediate image (after degibbs)
        "dwi_degibbs": processing_dir / "working_rpg.nii",
        "dwi_degibbs_bvec": processing_dir / "working_rpg.bvec",
        "dwi_degibbs_bval": processing_dir / "working_rpg.bval",
        # Intermediate image (after eddy/topup)
        "dwi_eddy": processing_dir / "dwiec.nii",
        "dwi_eddy_bvec": processing_dir / "dwiec.bvec",
        "dwi_eddy_bval": processing_dir / "dwiec.bval",
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
def white_matter_roi_list(paths) -> List[Nifti1Image]:
    return [create_binary_mask_from_fa(fa_file) for fa_file in paths["fa_dki"]]

@pytest.fixture(scope="module")
def white_matter_roi(white_matter_roi_list) -> Nifti1Image:
    wm_roi = reduce(np.bitwise_or, [roi.get_fdata().astype(bool) for roi in white_matter_roi_list])

    return Nifti1Image(wm_roi, white_matter_roi_list[0].affine, white_matter_roi_list[0].header)


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


def get_fa_test_params(ground_truth) -> List[Tuple[str, Dict, int]]:
    """Generate test parameters for FA stats from ground truth data."""
    fa_stats = ground_truth["fa_stats"]
    params = []
    
    for fa_type, values in fa_stats.items():
        if '_te' in fa_type:    # multi-TE case
            base_type, te_part = fa_type.split('_te')
            te_idx = te_list.index(float(te_part))
        else:   # TODO: not needed if there will be separate multi-TE test base class
            base_type = fa_type
            te_idx = 0
            
        params.append((base_type, values, te_idx))
            
    return params


@pytest.mark.parametrize("fa_type, expected_values, te_idx", get_fa_test_params(ground_truth))
def test_fa_stats(paths, white_matter_roi_list, fa_type, expected_values, te_idx):
    fa_path = paths[fa_type][te_idx]
    fa_image = nib.load(fa_path)
    fa_data = fa_image.get_fdata()
    fa_data_no_nan = np.nan_to_num(fa_data, nan=0.0, posinf=0.0, neginf=0.0)
    fa_image_no_nan = Nifti1Image(fa_data_no_nan, fa_image.affine, fa_image.header)
    
    white_matter_roi = white_matter_roi_list[te_idx]
    
    assert_roi_mean_and_std(fa_image_no_nan, white_matter_roi, expected_values["wm"])
    assert_roi_mean_and_std(fa_image_no_nan, paths["roi1"], expected_values["roi1"])
    assert_roi_mean_and_std(fa_image_no_nan, paths["roi2"], expected_values["roi2"])
    assert_roi_mean_and_std(fa_image_no_nan, paths["voxel"], [expected_values["voxel"]])
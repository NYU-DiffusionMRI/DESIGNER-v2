import subprocess
import shutil
from pathlib import Path

import numpy as np
import nibabel as nib
import pytest

from tests.utils import create_binary_mask_from_fa, extract_mean_b0, compute_roi_mean_and_std

@pytest.fixture(scope="module")
def paths():
    data_dir = Path("tests/data/D1/")

    tmp_dir = Path("tests/tmp/D1/")
    scratch_dir = tmp_dir / "processing"
    params_dir = tmp_dir / "params"

    return {
        # scratch and parameters directories for the pipeline
        "scratch": scratch_dir,
        "params": params_dir,
        # images to be processed by pipeline
        "dwi_images": [data_dir / "meso_slice_crop.nii", data_dir / "research_slice_crop.nii"],
        # Designer output image
        "dwi_designer": tmp_dir / "dwi_designer.nii",
        "bvec": tmp_dir / "dwi_designer.bvec",
        "bval": tmp_dir / "dwi_designer.bval",
        "b0": params_dir / "b0.nii",
        # FA images
        "fa_dki": params_dir / "fa_dki.nii",
        "fa_dti": params_dir / "fa_dti.nii",
        "fa_wdki": params_dir / "fa_wdki.nii",
        # ROI images
        "roi1": data_dir / "roi1.nii",
        "roi2": data_dir / "roi2.nii",
        "voxel": data_dir / "voxel.nii",
        # Intermediate image (after denoising)
        "dwidn": scratch_dir / "dwidn.nii",
        "dwidn_bvec": scratch_dir / "dwidn.bvec",
        "dwidn_bval": scratch_dir / "dwidn.bval",
        # Intermediate image (after degibbs)
        "dwidg": scratch_dir / "working_rpg.nii",
        "dwidg_bvec": scratch_dir / "working_rpg.bvec",
        "dwidg_bval": scratch_dir / "working_rpg.bval",
        # Intermediate image (after b1correct)
        "dwibc": scratch_dir / "dwibc.nii",
        "dwibc_bvec": scratch_dir / "dwibc.bvec",
        "dwibc_bval": scratch_dir / "dwibc.bval",
        # Intermediate image (after rician correction)
        "dwirc": scratch_dir / "dwirc.nii",
        "dwirc_bvec": scratch_dir / "dwirc.bvec",
        "dwirc_bval": scratch_dir / "dwirc.bval",
    }


@pytest.fixture(scope="module", autouse=True)
def run_pipeline(paths):
    dwi_args=",".join([str(path) for path in paths["dwi_images"]])

    processing_dir=paths["scratch"]
    params_dir=paths["params"]
    designer_image_path = paths["dwi_designer"]

    if processing_dir.exists():
        shutil.rmtree(processing_dir)

    subprocess.run(["designer", "-denoise", "-shrinkage", "frob", "-adaptive_patch", "-rician", "-degibbs", "-b1correct", "-normalize", "-mask", "-scratch", str(processing_dir), "-nocleanup", dwi_args, str(designer_image_path)])
    assert designer_image_path.exists()

    if params_dir.exists():
        shutil.rmtree(params_dir)

    subprocess.run(["tmi", "-SMI", "-DKI", "-WDKI", "-DTI", "-mask", str(processing_dir / "brain_mask.nii"), str(designer_image_path), str(params_dir)])
    assert params_dir.exists()

    yield

    # cleanup?


@pytest.fixture(scope="module")
def white_matter_roi(paths):
    return create_binary_mask_from_fa(paths["fa_dki"], threshold=0.3)


def test_white_matter_voxel_count(white_matter_roi):
    wm_voxel_cnt = np.count_nonzero(white_matter_roi)

    assert wm_voxel_cnt == 10789


# full pipeline
def test_b0_stats(paths, white_matter_roi):
    b0_data = extract_mean_b0(paths["dwi_designer"], paths["bval"])

    wm_mean, wm_std = compute_roi_mean_and_std(b0_data, white_matter_roi)
    assert np.isclose(wm_mean, 231.368)
    assert np.isclose(wm_std, 105.784, rtol=2e-3)

    roi1 = nib.load(paths["roi1"]).get_fdata()
    roi1_mean, roi1_std = compute_roi_mean_and_std(b0_data, roi1)
    assert np.isclose(roi1_mean, 222.263)
    assert np.isclose(roi1_std, 29.9254, rtol=2e-3)
    
    roi2 = nib.load(paths["roi2"]).get_fdata()
    roi2_mean, roi2_std = compute_roi_mean_and_std(b0_data, roi2)
    assert np.isclose(roi2_mean, 237.977)
    assert np.isclose(roi2_std, 31.7704, rtol=2e-3)

    single_voxel = nib.load(paths["voxel"]).get_fdata()
    single_voxel_mean, _ = compute_roi_mean_and_std(b0_data, single_voxel)
    assert np.isclose(single_voxel_mean, 192.126)


@pytest.mark.parametrize("fa_type, expected_values", [
    ("fa_dti", {
        "wm": (0.455159, 0.17331),
        "roi1": (0.701169, 0.123107),
        "roi2": (0.615107, 0.0898764),
        "voxel": 0.867286
    }),
    ("fa_dki", {
        "wm": (0.520315, 0.1688541),
        "roi1": (0.75445, 0.126123),
        "roi2": (0.681549, 0.108101),
        "voxel": 0.94263
    }),
    ("fa_wdki", {
        "wm": (0.520315, 0.1688541),
        "roi1": (0.75445, 0.126123),
        "roi2": (0.681549, 0.108101),
        "voxel": 0.94263
    })
])
def test_fa_stats(paths, white_matter_roi, fa_type, expected_values):
    fa_data = nib.load(paths[fa_type]).get_fdata()
    fa_data = np.nan_to_num(fa_data, nan=0.0)

    wm_mean, wm_std = compute_roi_mean_and_std(fa_data, white_matter_roi)
    assert np.isclose(wm_mean, expected_values["wm"][0])
    assert np.isclose(wm_std, expected_values["wm"][1], rtol=2e-3)

    roi1 = nib.load(paths["roi1"]).get_fdata()
    roi1_mean, roi1_std = compute_roi_mean_and_std(fa_data, roi1)
    assert np.isclose(roi1_mean, expected_values["roi1"][0])
    assert np.isclose(roi1_std, expected_values["roi1"][1], rtol=2e-3)
    
    roi2 = nib.load(paths["roi2"]).get_fdata()
    roi2_mean, roi2_std = compute_roi_mean_and_std(fa_data, roi2)
    assert np.isclose(roi2_mean, expected_values["roi2"][0])
    assert np.isclose(roi2_std, expected_values["roi2"][1], rtol=2e-3)

    single_voxel = nib.load(paths["voxel"]).get_fdata()
    single_voxel_mean, _ = compute_roi_mean_and_std(fa_data, single_voxel)
    assert np.isclose(single_voxel_mean, expected_values["voxel"])


@pytest.mark.parametrize("dwi_key, bval_key, expected_values", [
    ("dwidn", "dwidn_bval", {
        "wm": (227.365, 111.864),
        "roi1": (210.935, 23.1792),
        "roi2": (237.642, 33.3991),
        "voxel": 212.14
    }),
    ("dwidg", "dwidg_bval", {
        "wm": (229.873, 105.359),
        "roi1": (219.175, 29.4502),
        "roi2": (237.716, 32.1328),
        "voxel": 190.403
    }),
    ("dwibc", "dwibc_bval", {
        "wm": (233.302, 107.624),
        "roi1": (223.835, 30.2657),
        "roi2": (242.031, 32.2505),
        "voxel": 197.618
    }),
    ("dwirc", "dwirc_bval", {
        "wm": (232.783, 107.745),
        "roi1": (223.161, 30.3397),
        "roi2": (241.496, 32.2885),
        "voxel": 196.919
    })
])
def test_b0_stats_each_step(paths, white_matter_roi, dwi_key, bval_key, expected_values):
    b0_data = extract_mean_b0(paths[dwi_key], paths[bval_key])

    wm_mean, wm_std = compute_roi_mean_and_std(b0_data, white_matter_roi)
    assert np.isclose(wm_mean, expected_values["wm"][0])
    assert np.isclose(wm_std, expected_values["wm"][1], rtol=2e-3)

    roi1 = nib.load(paths["roi1"]).get_fdata()
    roi1_mean, roi1_std = compute_roi_mean_and_std(b0_data, roi1)
    assert np.isclose(roi1_mean, expected_values["roi1"][0])
    assert np.isclose(roi1_std, expected_values["roi1"][1], rtol=2e-3)

    roi2 = nib.load(paths["roi2"]).get_fdata()
    roi2_mean, roi2_std = compute_roi_mean_and_std(b0_data, roi2)
    assert np.isclose(roi2_mean, expected_values["roi2"][0])
    assert np.isclose(roi2_std, expected_values["roi2"][1], rtol=2e-3)

    single_voxel = nib.load(paths["voxel"]).get_fdata()
    single_voxel_mean, _ = compute_roi_mean_and_std(b0_data, single_voxel)
    assert np.isclose(single_voxel_mean, expected_values["voxel"])

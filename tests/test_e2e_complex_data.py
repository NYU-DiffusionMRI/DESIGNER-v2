import json
from pathlib import Path

import pytest

from tests.e2e_runner import prepare_complex_data_e2e_runner, E2ERunner
from tests.types import DWIStage, FAType, StatsDict
from tests.utils import assert_stats


ground_truth = json.load(open("tests/ground_truth_statistics/D5_complex_data.json"))


@pytest.fixture(scope="module", autouse=True)
def pipeline():
    test_dir = Path("tests/")
    data_dir = test_dir / "data/D5/"
    scratch_dir = test_dir / "tmp_D5/"

    test_runner = prepare_complex_data_e2e_runner(scratch_dir, data_dir)

    try:
        input_files = ["dwi1_ds.nii.gz", "dwi2_ds.nii.gz", "dwi3_ds.nii.gz"]
        input_paths = [data_dir / file for file in input_files]

        test_runner.run(input_paths)

        yield test_runner

    finally:
        pass
        # test_runner.cleanup()


def test_white_matter_voxel_count(pipeline: E2ERunner):
    wm_voxel_cnt = pipeline.count_white_matter_voxels()
    expected_cnt = ground_truth["white_matter_voxel_count"]
    assert wm_voxel_cnt == expected_cnt


@pytest.mark.parametrize("preproc_stage, expected_stats", ground_truth["b0_stats"].items())
def test_b0_stats(pipeline: E2ERunner, preproc_stage: DWIStage, expected_stats: StatsDict):
    b0_stats = pipeline.compute_b0_roi_stats(preproc_stage)
    assert_stats(b0_stats, expected_stats)


def parse_fa_with_te(fa_with_te: str) -> tuple[str, float]:
    fa_type, echo_time = fa_with_te.split("_te")
    return fa_type, float(echo_time)


fa_test_cases = [
    (*parse_fa_with_te(fa_with_te), gt_stats) 
    for fa_with_te, gt_stats in ground_truth["fa_stats"].items()
]

@pytest.mark.parametrize("fa_type, echo_time, expected_stats", fa_test_cases)
def test_fa_stats(pipeline: E2ERunner, fa_type: FAType, echo_time: float, expected_stats: StatsDict):
    fa_stats = pipeline.compute_fa_roi_stats(fa_type, echo_time)
    assert_stats(fa_stats, expected_stats)

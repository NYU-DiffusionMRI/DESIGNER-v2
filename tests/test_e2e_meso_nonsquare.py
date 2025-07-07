import json
from pathlib import Path

import pytest

from tests.e2e_runner import E2ERunner
from tests.e2e_runner_factory import prepare_meso_nonsquare_e2e_runner
from tests.types import DWIStage, FAType, StatsDict
from tests.utils import assert_stats


ground_truth = json.load(open("tests/benchmark/D1_meso_nonsquare.json"))


@pytest.fixture(scope="module", autouse=True)
def pipeline():
    test_dir = Path("tests/")
    data_dir = test_dir / "data/D1/"
    scratch_dir = test_dir / "tmp_D1/"

    test_runner: E2ERunner = prepare_meso_nonsquare_e2e_runner(scratch_dir, data_dir)

    try:
        input_files = ["meso_slice_crop.nii.gz", "research_slice_crop.nii.gz"]
        input_paths = [data_dir / file for file in input_files]

        test_runner.run(input_paths)

        yield test_runner

    finally:
        test_runner.cleanup()


def test_white_matter_voxel_count(pipeline: E2ERunner):
    wm_voxel_cnt = pipeline.count_white_matter_voxels()
    expected_cnt = ground_truth["white_matter_voxel_count"]
    assert wm_voxel_cnt == expected_cnt


@pytest.mark.parametrize("preproc_stage, expected_stats", ground_truth["b0_stats"].items())
def test_b0_stats(pipeline: E2ERunner, preproc_stage: DWIStage, expected_stats: StatsDict):
    b0_stats = pipeline.compute_b0_roi_stats(preproc_stage)
    assert_stats(b0_stats, expected_stats)


@pytest.mark.parametrize("fa_type, expected_stats", ground_truth["fa_stats"].items())
def test_fa_stats(pipeline: E2ERunner, fa_type: FAType, expected_stats: StatsDict):
    fa_stats = pipeline.compute_fa_roi_stats(fa_type)
    assert_stats(fa_stats, expected_stats)

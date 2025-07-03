import json
from pathlib import Path

import pytest

from tests.utils import assert_stats
from tests.e2e_runner import MESONonsquareRunner
from tests.types import DWIStage, FAType, StatsDict
from tests.e2e_runner import DesignerRunner, TMIRunner


ground_truth = json.load(open("tests/ground_truth_statistics/D1_MESO_nonsquare.json"))


@pytest.fixture(scope="module", autouse=True)
def runner():

    tmp_dir = Path("tests") / "tmp_D1"

    designer = DesignerRunner.init_meso_nonsquare(tmp_dir)
    tmi = TMIRunner.init_all_fa(tmp_dir)
    meso_runner = MESONonsquareRunner(designer, tmi, tmp_dir=tmp_dir)

    try:
        meso_runner.run_pipeline()
        yield meso_runner

    finally:
        meso_runner.cleanup()


def test_white_matter_voxel_count(runner: MESONonsquareRunner):
    wm_voxel_cnt = runner.count_white_matter_voxels()
    expected_cnt = ground_truth["white_matter_voxel_count"]
    assert wm_voxel_cnt == expected_cnt


@pytest.mark.parametrize("preproc_stage, expected_stats", ground_truth["b0_stats"].items())
def test_b0_stats(runner: MESONonsquareRunner, preproc_stage: DWIStage, expected_stats: StatsDict):
    b0_stats = runner.compute_dwi_b0_stats(preproc_stage)
    assert_stats(b0_stats, expected_stats)


@pytest.mark.parametrize("fa_type, expected_stats", ground_truth["fa_stats"].items())
def test_fa_stats(runner: MESONonsquareRunner, fa_type: FAType, expected_stats: StatsDict):
    fa_stats = runner.compute_fa_stats(fa_type)
    assert_stats(fa_stats, expected_stats)

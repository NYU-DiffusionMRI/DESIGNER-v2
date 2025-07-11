import json
from pathlib import Path

import pytest

from tests.e2e_runner import E2ERunner, StatsComputer
from tests.e2e_runner_factory import prepare_meso_eddy_e2e_runner
from tests.types import DWIStage, DiffusionModelType, StatsDict
from tests.utils import assert_stats


ground_truth = json.load(open("tests/benchmark/D2_meso_eddy_wo_bids.json"))

test_dir = Path("tests/")
data_dir = test_dir / "data/D2/"
scratch_dir = test_dir / "tmp_D2_wo_bids/"


@pytest.fixture(scope="module", autouse=True)
def runner(request: pytest.FixtureRequest):
    test_runner: E2ERunner = prepare_meso_eddy_e2e_runner(scratch_dir, data_dir, without_bids=True)

    try:
        test_runner.run()
        yield test_runner

    finally:
        no_cleanup = request.config.getoption("--no-cleanup")

        if not no_cleanup:
            test_runner.cleanup()


@pytest.fixture(scope="module")
def stats_computer(runner: E2ERunner):
    stats_computer = StatsComputer(runner, roi_dir=data_dir)
    return stats_computer


def test_white_matter_voxel_count(stats_computer: StatsComputer):
    wm_voxel_cnt = stats_computer.count_white_matter_voxels()
    expected_cnt = ground_truth["white_matter_voxel_count"]
    assert wm_voxel_cnt == expected_cnt


@pytest.mark.parametrize("preproc_stage, expected_stats", ground_truth["b0_stats"].items())
def test_b0_stats(stats_computer: StatsComputer, preproc_stage: DWIStage, expected_stats: StatsDict):
    b0_stats = stats_computer.compute_b0_roi_stats(preproc_stage)
    assert_stats(b0_stats, expected_stats)


@pytest.mark.parametrize("fa_model, expected_stats", ground_truth["fa_stats"].items())
def test_fa_stats(stats_computer: StatsComputer, fa_model: DiffusionModelType, expected_stats: StatsDict):
    fa_stats = stats_computer.compute_fa_roi_stats(fa_model)
    assert_stats(fa_stats, expected_stats)

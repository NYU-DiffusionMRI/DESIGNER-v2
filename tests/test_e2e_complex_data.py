import json
from pathlib import Path

import pytest

from tests.e2e_runner import E2ERunner, StatsComputer
from tests.e2e_runner_factory import prepare_complex_data_e2e_runner
from tests.types import DWIStage, DiffusionModelType, StatsDict
from tests.utils import assert_stats


ground_truth = json.load(open("tests/benchmark/D5_complex_data.json"))

test_dir = Path("tests/")
data_dir = test_dir / "data/D5/"
scratch_dir = test_dir / "tmp_D5/"


@pytest.fixture(scope="module", autouse=True)
def runner(request: pytest.FixtureRequest):
    test_runner: E2ERunner = prepare_complex_data_e2e_runner(scratch_dir, data_dir)

    try:
        test_runner.run()
        yield test_runner

    finally:
        no_cleanup = request.config.getoption("--no-cleanup")

        if not no_cleanup:
            test_runner.cleanup()


@pytest.fixture(scope="module")
def stats_computer(runner: E2ERunner):
    valid_echo_times = [60.0]
    stats_computer = StatsComputer(runner, roi_dir=data_dir, valid_echo_times=valid_echo_times)
    return stats_computer


def test_white_matter_voxel_count(stats_computer: StatsComputer):
    wm_voxel_cnt = stats_computer.count_white_matter_voxels()
    expected_cnt = ground_truth["white_matter_voxel_count"]
    assert wm_voxel_cnt == expected_cnt


@pytest.mark.parametrize("preproc_stage, expected_stats", ground_truth["b0_stats"].items())
def test_b0_stats(stats_computer: StatsComputer, preproc_stage: DWIStage, expected_stats: StatsDict):
    b0_stats = stats_computer.compute_b0_roi_stats(preproc_stage)
    assert_stats(b0_stats, expected_stats)


def parse_fa_with_te(fa_with_te: str) -> tuple[str, float]:
    fa_model, echo_time = fa_with_te.split("_te")
    return fa_model, float(echo_time)


fa_test_cases = [
    (*parse_fa_with_te(fa_with_te), gt_stats) 
    for fa_with_te, gt_stats in ground_truth["fa_stats"].items()
]

@pytest.mark.parametrize("fa_model, echo_time, expected_stats", fa_test_cases)
def test_fa_stats(stats_computer: StatsComputer, fa_model: DiffusionModelType, echo_time: float, expected_stats: StatsDict):
    fa_stats = stats_computer.compute_fa_roi_stats(fa_model, echo_time)
    assert_stats(fa_stats, expected_stats)

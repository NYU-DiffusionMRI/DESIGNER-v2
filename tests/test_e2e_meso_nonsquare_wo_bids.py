import json
from pathlib import Path

import numpy as np
import pytest

from tests.e2e_runner import E2ERunner, StatsComputer
from tests.e2e_runner_factory import prepare_meso_nonsquare_e2e_runner
from tests.types import DWIStage, DiffusionModelType, StatsDict
from tests.utils import assert_stats


# D1 MESO non-square without BIDS benchmark is the same as the one with BIDS.
ground_truth = json.load(open("tests/benchmark/D1_meso_nonsquare.json"))
# ground_truth = json.load(open("benchmark_windows_x86_64/D1_meso_nonsquare.json"))

test_dir = Path("tests/")
data_dir = test_dir / "data/D1/"
scratch_dir = test_dir / "tmp_D1_wo_bids/"


@pytest.fixture(scope="module", autouse=True)
def runner(request: pytest.FixtureRequest):
    test_runner: E2ERunner = prepare_meso_nonsquare_e2e_runner(scratch_dir, data_dir, without_bids=True)

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


@pytest.fixture(scope="module")
def tolerances(tolerance_config: dict):
    profile = tolerance_config["dataset_profiles"]["D1"]
    return tolerance_config["tolerances"][profile]


def test_wm_ratio(stats_computer: StatsComputer, tolerances: dict[str, float]):
    print("E2E TEST: D1 MESO NON-SQUARE WO BIDS")
    wm_ratio = stats_computer.compute_wm_ratio()
    expected_ratio = ground_truth["wm_ratio"]
    print(f"wm_ratio: {wm_ratio}, expected_ratio: {expected_ratio}")
    assert np.isclose(wm_ratio, expected_ratio, atol=tolerances["wm_ratio_atol"])


@pytest.mark.parametrize("preproc_stage, expected_stats", ground_truth["b0_stats"].items())
def test_b0_stats(stats_computer: StatsComputer, preproc_stage: DWIStage, expected_stats: StatsDict, tolerances: dict[str, float]):
    print("E2E TEST: D1 MESO NON-SQUARE WO BIDS")
    b0_stats = stats_computer.compute_b0_roi_stats(preproc_stage)
    assert_stats(b0_stats, expected_stats, mean_tol=tolerances["b0_mean_rtol"], std_tol=tolerances["b0_std_rtol"], is_relative=True, context=preproc_stage)
    print()


@pytest.mark.parametrize("fa_model, expected_stats", ground_truth["fa_stats"].items())
def test_fa_stats(stats_computer: StatsComputer, fa_model: DiffusionModelType, expected_stats: StatsDict, tolerances: dict[str, float]):
    print("E2E TEST: D1 MESO NON-SQUARE WO BIDS")
    fa_stats = stats_computer.compute_fa_roi_stats(fa_model)
    assert_stats(fa_stats, expected_stats, mean_tol=tolerances["fa_mean_atol"], std_tol=tolerances["fa_std_atol"], is_relative=False, context=fa_model)
    print()

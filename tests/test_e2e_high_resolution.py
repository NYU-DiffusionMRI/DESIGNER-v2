import json
from pathlib import Path

import numpy as np
import pytest

from tests.e2e_runners import E2ERunner, StatsComputer
from tests.e2e_runner_factory import prepare_high_resolution_e2e_runner
from tests.types import DWIStage, DiffusionModelType, StatsDict, ToleranceProfile
from tests.utils import assert_stats


ground_truth = json.load(open("tests/benchmark/D4_high_resolution.json"))

test_dir = Path("tests/")
data_dir = test_dir / "data/D4/"
scratch_dir = test_dir / "tmp_D4/"


@pytest.fixture(scope="module", autouse=True)
def runner(request: pytest.FixtureRequest):
    test_runner: E2ERunner = prepare_high_resolution_e2e_runner(scratch_dir, data_dir)

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
def tol_profile(tolerance_config: dict) -> ToleranceProfile:
    profile = tolerance_config["dataset_profiles"]["D4"]
    return tolerance_config["tolerances"][profile]


def test_wm_ratio(stats_computer: StatsComputer, tol_profile: ToleranceProfile):
    print("E2E TEST: D4 HIGH RESOLUTION")
    wm_ratio = stats_computer.compute_wm_ratio()
    expected_ratio = ground_truth["wm_ratio"]
    print(f"wm_ratio: {wm_ratio}, expected_ratio: {expected_ratio}")
    assert np.isclose(wm_ratio, expected_ratio, atol=tol_profile["wm_ratio_atol"])


@pytest.mark.parametrize("preproc_stage, expected_stats", ground_truth["b0_stats"].items())
def test_b0_stats(stats_computer: StatsComputer, preproc_stage: DWIStage, expected_stats: StatsDict, tol_profile: ToleranceProfile):
    print("E2E TEST: D4 HIGH RESOLUTION")
    b0_stats = stats_computer.compute_b0_roi_stats(preproc_stage)
    assert_stats(b0_stats, expected_stats, tol_profile=tol_profile, context=preproc_stage)
    print()


@pytest.mark.parametrize("fa_model, expected_stats", ground_truth["fa_stats"].items())
def test_fa_stats(stats_computer: StatsComputer, fa_model: DiffusionModelType, expected_stats: StatsDict, tol_profile: ToleranceProfile):
    print("E2E TEST: D4 HIGH RESOLUTION")
    fa_stats = stats_computer.compute_fa_roi_stats(fa_model)
    assert_stats(fa_stats, expected_stats, tol_profile=tol_profile, context=fa_model)
    print()

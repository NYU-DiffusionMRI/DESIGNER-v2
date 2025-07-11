import argparse
import multiprocessing
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Callable, Dict, List, Tuple

from tests.e2e_runner import E2ERunner, StatsComputer
from tests.e2e_runner_factory import (
    prepare_meso_nonsquare_e2e_runner,
    prepare_meso_eddy_e2e_runner,
    prepare_complex_data_e2e_runner,
    prepare_heal_coronal_e2e_runner,
    prepare_meso_degibbs_e2e_runner,
    prepare_high_resolution_e2e_runner,
)
from tests.types import DWIStage, DiffusionModelType


# No need of MESO non-square without BIDS benchmark since it is the same as the one with BIDS.
def save_meso_nonsquare_benchmark(save_path: Path):
    test_dir = Path("tests/")
    data_dir = test_dir / "data/D1/"
    scratch_dir = test_dir / "tmp_D1_benchmark/"
    benchmark_runner: E2ERunner | None = None

    try:
        benchmark_runner = prepare_meso_nonsquare_e2e_runner(scratch_dir, data_dir)
        benchmark_runner.run()

        stages: List[DWIStage] = ["denoising", "degibbs", "b1correct", "rician", "designer"]
        fa_models: List[DiffusionModelType] = ["dti", "dki", "wdki"]

        stats_computer = StatsComputer(benchmark_runner, roi_dir=data_dir)
        stats_computer.save_benchmark(save_path, stages, fa_models)
        
        return benchmark_runner

    except Exception as e:
        raise RuntimeError(f"Failed to run meso nonsquare benchmark: {str(e)}") from e
    finally:
        if benchmark_runner is not None:
            benchmark_runner.cleanup()


def save_meso_eddy_benchmark(save_path: Path, *, without_bids: bool = False):
    test_dir = Path("tests/")
    data_dir = test_dir / "data/D2/"
    benchmark_runner: E2ERunner | None = None

    if without_bids:
        scratch_dir = test_dir / "tmp_D2_wo_bids_benchmark/"
    else:
        scratch_dir = test_dir / "tmp_D2_benchmark/"
   
    try:
        benchmark_runner = prepare_meso_eddy_e2e_runner(scratch_dir, data_dir, without_bids=without_bids)
        benchmark_runner.run()

        stages: List[DWIStage] = ["topup", "designer"]
        fa_models: List[DiffusionModelType] = ["dti", "dki", "wdki"]
        
        stats_computer = StatsComputer(benchmark_runner, roi_dir=data_dir)
        stats_computer.save_benchmark(save_path, stages, fa_models)
        
        return benchmark_runner

    except Exception as e:
        raise RuntimeError(f"Failed to run meso eddy benchmark: {str(e)}") from e
    finally:
        if benchmark_runner is not None:
            benchmark_runner.cleanup()


def save_meso_degibbs_benchmark(save_path: Path):
    test_dir = Path("tests/")
    data_dir = test_dir / "data/D3/"
    scratch_dir = test_dir / "tmp_D3_benchmark/"
    benchmark_runner: E2ERunner | None = None

    try:
        benchmark_runner = prepare_meso_degibbs_e2e_runner(scratch_dir, data_dir)
        benchmark_runner.run()

        stages: List[DWIStage] = ["designer"]
        fa_models: List[DiffusionModelType] = ["dti", "dki", "wdki"]

        stats_computer = StatsComputer(benchmark_runner, roi_dir=data_dir)
        stats_computer.save_benchmark(save_path, stages, fa_models)
        
        return benchmark_runner

    except Exception as e:
        raise RuntimeError(f"Failed to run meso degibbs benchmark: {str(e)}") from e
    finally:
        if benchmark_runner is not None:
            benchmark_runner.cleanup()


def save_high_resolution_benchmark(save_path: Path):
    test_dir = Path("tests/")
    data_dir = test_dir / "data/D4/"
    scratch_dir = test_dir / "tmp_D4_benchmark/"
    benchmark_runner: E2ERunner | None = None

    try:
        benchmark_runner = prepare_high_resolution_e2e_runner(scratch_dir, data_dir)
        benchmark_runner.run()

        stages: List[DWIStage] = ["denoising", "degibbs", "designer"]
        fa_models: List[DiffusionModelType] = ["dti", "dki", "wdki"]

        stats_computer = StatsComputer(benchmark_runner, roi_dir=data_dir)
        stats_computer.save_benchmark(save_path, stages, fa_models)
        
        return benchmark_runner

    except Exception as e:
        raise RuntimeError(f"Failed to run high resolution benchmark: {str(e)}") from e
    finally:
        if benchmark_runner is not None:
            benchmark_runner.cleanup()


def save_complex_data_benchmark(save_path: Path):
    test_dir = Path("tests/")
    data_dir = test_dir / "data/D5/"
    scratch_dir = test_dir / "tmp_D5_benchmark/"
    benchmark_runner: E2ERunner | None = None

    try:
        benchmark_runner = prepare_complex_data_e2e_runner(scratch_dir, data_dir)
        benchmark_runner.run()

        stages: List[DWIStage] = ["denoising", "degibbs", "eddy", "designer"]
        fa_models: List[DiffusionModelType] = ["dti", "dki", "wdki"]

        valid_echo_times = [60.0]
        stats_computer = StatsComputer(benchmark_runner, roi_dir=data_dir, valid_echo_times=valid_echo_times)
        stats_computer.save_benchmark(save_path, stages, fa_models, valid_echo_times)
        
        return benchmark_runner

    except Exception as e:
        raise RuntimeError(f"Failed to run complex data benchmark: {str(e)}") from e
    finally:
        if benchmark_runner is not None:
            benchmark_runner.cleanup()


def save_heal_coronal_benchmark(save_path: Path):
    test_dir = Path("tests/")
    data_dir = test_dir / "data/D6/"
    scratch_dir = test_dir / "tmp_D6_benchmark/"
    benchmark_runner: E2ERunner | None = None

    try:
        benchmark_runner = prepare_heal_coronal_e2e_runner(scratch_dir, data_dir)
        benchmark_runner.run()

        stages: List[DWIStage] = ["denoising", "degibbs", "designer"]
        fa_models: List[DiffusionModelType] = ["dti"]

        stats_computer = StatsComputer(benchmark_runner, roi_dir=data_dir, with_skull_stripping=True)
        stats_computer.save_benchmark(save_path, stages, fa_models)
        
        return benchmark_runner

    except Exception as e:
        raise RuntimeError(f"Failed to run heal coronal benchmark: {str(e)}") from e
    finally:
        if benchmark_runner is not None:
            benchmark_runner.cleanup()


# This wrapper is needed since the labmda function cannot be pickled (not serializable) during multiprocessing.
def save_meso_eddy_wo_bids_benchmark(save_path: Path):
    """Wrapper function for save_meso_eddy_benchmark with without_bids=True"""
    return save_meso_eddy_benchmark(save_path, without_bids=True)


def get_available_benchmarks() -> Dict[str, Tuple[Callable, str, str]]:
    """Returns a mapping of benchmark IDs to their (function, filename, description)."""
    return {
        "D1": (
            save_meso_nonsquare_benchmark,
            "D1_meso_nonsquare.json",
            "MESO non-square matrix"
        ),
        "D2": (
            save_meso_eddy_benchmark,
            "D2_meso_eddy.json",
            "MESO with eddy"
        ),
        "D2_wo_bids": (
            save_meso_eddy_wo_bids_benchmark,
            "D2_meso_eddy_wo_bids.json",
            "MESO with eddy (without BIDS)"
        ),
        "D3": (
            save_meso_degibbs_benchmark,
            "D3_meso_degibbs.json",
            "MESO with degibbs"
        ),
        "D4": (
            save_high_resolution_benchmark,
            "D4_high_resolution.json",
            "High resolution"
        ),
        "D5": (
            save_complex_data_benchmark,
            "D5_complex_data.json",
            "Complex data"
        ),
        "D6": (
            save_heal_coronal_benchmark,
            "D6_heal_coronal.json",
            "HEAL coronal data"
        )
    }

def parse_args() -> argparse.Namespace:
    available_benchmarks = get_available_benchmarks()
    
    # Generate help text dynamically from available benchmarks
    benchmark_descriptions = [
        f"  {benchmark_id}: {description}"
        for benchmark_id, (_, _, description) in sorted(available_benchmarks.items())
    ]
    benchmark_help = (
        'List of benchmarks to run or "all" for all benchmarks. Available benchmarks:\n' +
        '\n'.join(benchmark_descriptions)
    )
    
    parser = argparse.ArgumentParser(
        description="Generate benchmarks for DESIGNER tests",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Directory where benchmark files will be saved"
    )
    parser.add_argument(
        "--benchmarks",
        type=str,
        nargs="+",
        default=["all"],
        help=benchmark_help
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=multiprocessing.cpu_count() // 2,
        help=f"Number of CPU cores to use (default: half of available cores)"
    )
    return parser.parse_args()


def run_single_benchmark(func: Callable, save_path: Path) -> Tuple[str, bool, str]:
    """Run a single benchmark with proper error handling and cleanup.
    Returns: (benchmark_name, success, error_message)"""
    try:
        func(save_path)  # Don't store the runner - it's cleaned up in the function
        return (save_path.name, True, "")
    except Exception as e:
        return (save_path.name, False, str(e))


def main():
    start_time = time.time()
    
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    available_benchmarks = get_available_benchmarks()
    selected_benchmarks = list(available_benchmarks.keys()) if "all" in args.benchmarks else args.benchmarks
    
    # Validate benchmark selection
    invalid_benchmarks = [b for b in selected_benchmarks if b not in available_benchmarks]
    if invalid_benchmarks:
        raise ValueError(f"Invalid benchmark(s): {', '.join(invalid_benchmarks)}. "
                        f"Available benchmarks are: {', '.join(available_benchmarks.keys())}")
    
    # Prepare benchmark list
    benchmarks = [
        (available_benchmarks[benchmark_id][0], output_dir / available_benchmarks[benchmark_id][1])
        for benchmark_id in selected_benchmarks
    ]
    
    print(f"🚀 Running {len(benchmarks)} benchmark(s)...")
    print(f"📁 Output directory: {output_dir}")
    print(f"📊 Selected benchmarks:")
    for benchmark_id in selected_benchmarks:
        _, _, description = available_benchmarks[benchmark_id]
        print(f"   • {benchmark_id}: {description}")

    max_workers = args.cores
    print(f"👥 Using {max_workers} worker processes")

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for func, save_path in benchmarks:
            futures.append(executor.submit(run_single_benchmark, func, save_path))

        completed = 0
        total = len(futures)
        failed_benchmarks = []

        for future in as_completed(futures):
            benchmark_name, success, error_msg = future.result()
            completed += 1
            if success:
                print(f"✅ Benchmark {completed}/{total} finished successfully: {benchmark_name}")
            else:
                print(f"❌ Benchmark {completed}/{total} failed: {benchmark_name}")
                print(f"   Error: {error_msg}")
                failed_benchmarks.append(benchmark_name)

    elapsed_time = time.time() - start_time
    minutes = int(elapsed_time // 60)
    seconds = int(elapsed_time % 60)
    
    print(f"\n✨ All selected benchmarks completed in {minutes}m {seconds}s")
    print(f"📁 Results saved in: {output_dir}")
    
    if failed_benchmarks:
        print(f"\n⚠️  {len(failed_benchmarks)} benchmark(s) failed:")
        for name in failed_benchmarks:
            print(f"   - {name}")
        sys.exit(1)


if __name__ == "__main__":
    main()

# End-to-End Testing Setup

This directory provides the end-to-end testing setup for DESIGNER-v2. It checks the full preprocessing pipeline on a variety of datasets and scenarios. These tests run automatically as part of the CI pipeline to help keep the pipeline stable and catch any bugs that might be introduced by future changes.

## Overview

The testing suite is designed to:
- Validate preprocessing steps. (denoising, degibbs, eddy correction, etc.)
- Test various input data configurations and scenarios.
- Ensure consistent results within the same platform. (e.g., Mac ARM64, and Windows x86_64)
- Ensure stable results across different platforms.
- Compare outputs against established benchmarks.
- Ease the process of updating benchmarks when app logic changes.

## Directory Structure

```
tests/
├── benchmark/               # Ground truth benchmark data
│   ├── D1_meso_nonsquare.json
│   ├── D1_meso_nonsquare_wo_bids.json
│   ├── D2_meso_eddy.json
│   └── ...
├── data/                   # Test datasets
│   ├── D1/                 # MESO non-square data
│   ├── D2/                 # MESO eddy data
│   └── ...
├── models/                 # ML models for testing
│   └── synthstrip_v7.3.2.pt
├── test_e2e_*.py           # Individual test modules for each benchmark
├── e2e_runners.py          # Base classes for test execution
├── e2e_runner_factory.py   # Factory functions for DESIGNER and TMI configurations
├── types.py                # Type definitions and validation
├── utils.py                # Common utility functions
├── tolerance_config.yaml   # Tolerance profiles for different datasets
└── requirements_test.txt   # Test dependencies
```

## Test Datasets

The suite includes tests for multiple datasets with different characteristics:

1. **D1: MESO Non-Square**
   - Tests non-square matrix acquisition
   - Validates with and without BIDS metadata

2. **D2: MESO Eddy**
   - Focuses on eddy current correction
   - Tests both BIDS and non-BIDS workflows
   - Uses downsampled data (0.38 factor) for CI pipeline efficiency
   - Has specific tolerance profile due to downsampling effects

3. **D3: MESO Degibbs**
   - Validates Gibbs ringing correction
   - Tests partial Fourier acquisition

4. **D4: High Resolution**
   - Tests processing of high-resolution data
   - Validates memory management

5. **D5: Complex Data**
   - Tests multi-shell acquisitions
   - Validates multiple echo times
   - Uses downsampled data (0.38 factor) for CI pipeline efficiency
   - Has specific tolerance profile due to downsampling effects

6. **D6: HEAL Coronal**
   - Tests coronal acquisition
   - Uses downsampled data (0.5 factor)
   - Validates skull stripping using Freesurfer's synthstrip model

## Running Tests

### Prerequisites

Test dependencies are automatically installed in the Docker dev container. If not using the dev container, run the following command in the project root to install the required dependencies:
```bash
pip install -r tests/requirements_test.txt
```

### Running Tests Locally

From the project root:
```bash
# Run all tests
pytest -s -v [--no-cleanup] tests/

# Faster execution: run in parallel (recommended)
pytest -s -v --dist loadscope -n $NUM_CORES [--no-cleanup] tests/

# Run a specific test module (e.g., D1 MESO non-square)
pytest -s -v tests/test_e2e_D1_meso_nonsquare.py
```
- `--no-cleanup` option preserves temporary files (DESIGNER's `processing/` directory and TMI's `params/` directory) for debugging.
- `$NUM_CORES` is the number of cores to use for parallel execution. Can be set to `auto` to use all available cores.


## Benchmarks and Tolerances

### Benchmark Ground Truth

Benchmark data is stored as JSON files in `tests/benchmark/` and includes:
- **White matter ratio** measurements
- **B0 statistics** for different preprocessing steps (denoising, eddy correction, etc.)
- **FA statistics** for different models (DTI, DKI, WDKI)

cf.) B0 and FA statistics are computed from the same white matter ROI, two randomly drawn ROIs and one random voxel. Please refer to `tests/benchmark/D1_meso_nonsquare.json` for an example.

### Tolerance Configuration for Different Datasets

The `tolerance_config.yaml` file manages test tolerances counting benchmark variations across different platforms:
- Tolerances are constructed based on variations between Mac ARM64, Windows x86_64, and Linux x86_64.
- Higer tolerances are used for D2 and D5 datasets to account for platform variations due to downsampling effects.
- **Relative tolerances** are used for B0 statistics and **absolute tolerances** are used for FA statistics and white matter ratio.

### Updating Benchmarks

When DESIGNER or TMI app logic changes, benchmarks may need to be updated:

1. Verify the code changes are correct and intended.
2. Generate new benchmarks using:
   ```bash
   # Generate all benchmarks
   python tests/scripts/generate_test_benchmark.py --cores $NUM_CORES --output-dir tests/benchmark
   
   # Generate benchmark for specific benchmarks (e.g., D1_wo_bids and D2)
   python tests/scripts/generate_test_benchmark.py --benchmarks D1_wo_bids D2 --cores $NUM_CORES --output-dir tests/benchmark

   # For more details
   python tests/scripts/generate_test_benchmark.py --help
   ```

## Adding a New E2E Test Configuration

1. Create test data in `tests/data/` (DWI files, ROI masks, etc. in **.nii.gz** format).
2. Create factory function in `e2e_runner_factory.py` for the new test configuration.
3. Add benchmark saving logic in `tests/scripts/generate_test_benchmark.py`.
4. Generate benchmark data in `tests/benchmark/` following the [Updating Benchmarks](#updating-benchmarks) section.
5. Add test module following existing patterns (prefix with `test_e2e_`).
6. Update tolerance configuration as needed experimenting on different platforms.

## Scripts

Helper scripts in `scripts/` (in the project root):
- `generate_test_benchmark.py`: Generate/update benchmark data when app logic changes.
- `download_freesurfer_synthstrip.sh`: Download Freesurfer's deep learning based synthstrip model for D6 skull stripping.

## Notes

- Sythnstrip model (`tests/models/synthstrip_v7.3.2.pt`), which is used for D6 skull stripping, is downloaded by running `bash scripts/download_freesurfer_synthstrip.sh`

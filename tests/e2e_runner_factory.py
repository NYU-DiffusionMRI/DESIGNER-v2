from pathlib import Path
from typing import List

from tests.e2e_runner import E2ERunner
from tests.types import DiffusionModelType


# BELOW ARE FACTORY FUNCTIONS FOR E2E RUNNERS

def prepare_meso_nonsquare_e2e_runner(scratch_dir: Path, data_dir: Path, *, without_bids: bool = False) -> E2ERunner:
    designer_config = [
        "-set_seed",
        "-denoise", 
        "-shrinkage", "frob", 
        "-adaptive_patch", 
        "-rician", 
        "-degibbs", 
        "-b1correct", 
        "-normalize"
    ]

    if without_bids:
        designer_config.extend([
            "-pf", "6/8",
            "-pe_dir", "AP",
        ])

    if without_bids:
        # input files are suffixed with 2 since DESIGNER automatically checks BIDS with the same input file names.
        input_files = ["meso_slice_crop2.nii.gz", "research_slice_crop2.nii.gz"]
    else:
        input_files = ["meso_slice_crop.nii.gz", "research_slice_crop.nii.gz"]
    input_paths = [data_dir / file for file in input_files]

    diffusion_models: List[DiffusionModelType] = ["smi", "dti", "dki", "wdki"]
    
    runner = E2ERunner(
        scratch_dir=scratch_dir,
        input_dwi_paths=input_paths,
        designer_cli_options=designer_config,
        diffusion_models=diffusion_models
    )

    return runner


def prepare_meso_eddy_e2e_runner(scratch_dir: Path, data_dir: Path, *, without_bids: bool = False) -> E2ERunner:
    rpe_image_path = data_dir / "pa_ds.nii.gz"
    
    designer_config = [
        "-set_seed", 
        "-eddy", 
        # TODO: passing relative path doesn't work. need to fix it in DESIGNER code.
        "-rpe_pair", str(rpe_image_path.resolve()),
    ]

    if without_bids:
        designer_config.extend([
            "-pf", "6/8",
            "-pe_dir", "AP",
        ])
        # input files are suffixed with 2 since DESIGNER automatically checks BIDS with the same input file names.
        input_files = ["meso_ds2.nii.gz", "research_ds2.nii.gz"]
    else:
        input_files = ["meso_ds.nii.gz", "research_ds.nii.gz"]
    input_paths = [data_dir / file for file in input_files]

    diffusion_models: List[DiffusionModelType] = ["smi", "dti", "dki", "wdki"]

    runner = E2ERunner(
        scratch_dir=scratch_dir,
        input_dwi_paths=input_paths,
        designer_cli_options=designer_config,
        diffusion_models=diffusion_models,
        with_sigma=False,
    )

    return runner


def prepare_meso_degibbs_e2e_runner(scratch_dir: Path, data_dir: Path) -> E2ERunner:
    designer_config = ["-set_seed", "-degibbs"]

    input_files = ["meso_slice.nii.gz", "research_slice.nii.gz"]
    input_paths = [data_dir / file for file in input_files]

    diffusion_models: List[DiffusionModelType] = ["smi", "dti", "dki", "wdki"]

    runner = E2ERunner(
        scratch_dir=scratch_dir,
        input_dwi_paths=input_paths,
        designer_cli_options=designer_config,
        diffusion_models=diffusion_models,
        with_sigma=False,
    )

    return runner


def prepare_high_resolution_e2e_runner(scratch_dir: Path, data_dir: Path) -> E2ERunner:
    phase_image_paths = [
        data_dir / "dif1phase_slice.nii.gz", 
        data_dir / "dif2phase_slice.nii.gz", 
    ]

    phase_args = ",".join([str(path) for path in phase_image_paths])
    
    designer_config = [
        "-set_seed",
        "-denoise", 
        "-phase", phase_args,
        "-degibbs", 
        "-normalize"
    ]

    input_files = ["dif1_slice.nii.gz", "dif2_slice.nii.gz"]
    input_paths = [data_dir / file for file in input_files]

    diffusion_models: List[DiffusionModelType] = ["smi", "dti", "dki", "wdki"]

    runner = E2ERunner(
        scratch_dir=scratch_dir,
        input_dwi_paths=input_paths,
        designer_cli_options=designer_config,
        diffusion_models=diffusion_models,
    )

    return runner


def prepare_complex_data_e2e_runner(scratch_dir: Path, data_dir: Path) -> E2ERunner:
    phase_image_paths = [
        data_dir / "phase1_ds.nii.gz", 
        data_dir / "phase2_ds.nii.gz", 
    ]

    rpe_image_path = data_dir / "rpe_ds.nii.gz"

    echo_times = [60.0, 78.0]

    phase_args = ",".join([str(path) for path in phase_image_paths])
    te_args = ",".join([str(time / 1000) for time in echo_times])   # convert ms to s
    
    designer_config = [
        "-set_seed", 
        "-denoise", 
        "-shrinkage", "frob", 
        "-phase", phase_args,
        "-degibbs", 
        "-eddy", 
        "-rpe_pair", str(rpe_image_path.resolve()), 
        "-eddy_fakeb", "0.85,1.2", 
        "-rpe_te", "60", 
        "-bshape", "1,0.6", 
        "-echo_time", te_args,
        "-normalize"
    ]

    input_files = ["dwi1_ds.nii.gz", "dwi2_ds.nii.gz"]
    input_paths = [data_dir / file for file in input_files]

    diffusion_models: List[DiffusionModelType] = ["smi", "dti", "dki", "wdki"]

    runner = E2ERunner(
        scratch_dir=scratch_dir,
        input_dwi_paths=input_paths,
        designer_cli_options=designer_config,
        diffusion_models=diffusion_models,
    )

    return runner


def prepare_heal_coronal_e2e_runner(scratch_dir: Path, data_dir: Path) -> E2ERunner:
    designer_config = [
        "-set_seed",
        "-denoise", 
        "-shrinkage", "frob", 
        "-adaptive_patch", 
        "-degibbs", 
        "-pf", "1", 
        "-normalize",
    ]

    input_files = ["dwi1_ds.nii.gz", "dwi2_ds.nii.gz"]
    input_paths = [data_dir / file for file in input_files]

    diffusion_models: List[DiffusionModelType] = ["dti"]

    runner = E2ERunner(
        scratch_dir=scratch_dir,
        input_dwi_paths=input_paths,
        designer_cli_options=designer_config,
        diffusion_models=diffusion_models,
        with_sigma=False,
        with_mask=True
    )

    return runner

from pathlib import Path

from tests.e2e_runner import E2ERunner, DesignerRunner, TMIRunner, StatsComputer


# BELOW ARE FACTORY FUNCTIONS FOR E2E RUNNERS

def prepare_meso_nonsquare_e2e_runner(scratch_dir: Path, data_dir: Path, *, without_bids: bool = False) -> E2ERunner:
    cmd_config = [
        "-denoise", 
        "-shrinkage", "frob", 
        "-adaptive_patch", 
        "-rician", 
        "-degibbs", 
        "-b1correct", 
        "-normalize"
    ]

    if without_bids:
        cmd_config.extend([
            "-pf", "6/8",
            "-pe_dir", "AP",
        ])
    
    designer = DesignerRunner(scratch_dir, cmd_config)
    tmi = TMIRunner.init_all_fa(scratch_dir)
    stats_computer = StatsComputer(designer, tmi, roi_dir=data_dir)
    runner = E2ERunner(designer, tmi, stats_computer, scratch_dir=scratch_dir)

    return runner


def prepare_meso_eddy_e2e_runner(scratch_dir: Path, data_dir: Path, *, without_bids: bool = False) -> E2ERunner:
    rpe_image_path = data_dir / "pa_ds.nii.gz"
    
    cmd_config = [
        "-set_seed", 
        "-eddy", 
        # TODO: passing relative path doesn't work. need to fix it in DESIGNER code.
        "-rpe_pair", str(rpe_image_path.resolve()),
    ]

    if without_bids:
        cmd_config.extend([
            "-pf", "6/8",
            "-pe_dir", "AP",
        ])

    designer = DesignerRunner(scratch_dir, cmd_config)
    tmi = TMIRunner.init_all_fa(scratch_dir, with_sigma=False)
    stats_computer = StatsComputer(designer, tmi, roi_dir=data_dir)
    runner = E2ERunner(designer, tmi, stats_computer, scratch_dir=scratch_dir)

    return runner


def prepare_meso_degibbs_e2e_runner(scratch_dir: Path, data_dir: Path) -> E2ERunner:
    cmd_config = ["-degibbs"]
    
    designer = DesignerRunner(scratch_dir, cmd_config)
    tmi = TMIRunner.init_all_fa(scratch_dir, with_sigma=False)
    stats_computer = StatsComputer(designer, tmi, roi_dir=data_dir)
    runner = E2ERunner(designer, tmi, stats_computer, scratch_dir=scratch_dir)

    return runner


def prepare_high_resolution_e2e_runner(scratch_dir: Path, data_dir: Path) -> E2ERunner:
    phase_image_paths = [
        data_dir / "dif1phase_slice.nii.gz", 
        data_dir / "dif2phase_slice.nii.gz", 
    ]

    phase_args = ",".join([str(path) for path in phase_image_paths])
    
    cmd_config = [
        "-denoise", 
        "-phase", phase_args,
        "-degibbs", 
        "-normalize"
    ]
    
    designer = DesignerRunner(scratch_dir, cmd_config)
    tmi = TMIRunner.init_all_fa(scratch_dir)
    stats_computer = StatsComputer(designer, tmi, roi_dir=data_dir)
    runner = E2ERunner(designer, tmi, stats_computer, scratch_dir=scratch_dir)

    return runner


def prepare_complex_data_e2e_runner(scratch_dir: Path, data_dir: Path) -> E2ERunner:
    phase_image_paths = [
        data_dir / "phase1_ds.nii.gz", 
        data_dir / "phase2_ds.nii.gz", 
        data_dir / "phase3_ds.nii.gz"
    ]

    rpe_image_path = data_dir / "rpe_ds.nii.gz"

    echo_times = [60.0, 78.0, 92.0]
    valid_echo_times = [60.0, 92.0]

    phase_args = ",".join([str(path) for path in phase_image_paths])
    te_args = ",".join([str(time / 1000) for time in echo_times])   # convert ms to s
    
    cmd_config = [
    "-set_seed", 
        "-denoise", 
        "-shrinkage", "frob", 
        "-phase", phase_args,
        "-degibbs", 
        "-eddy", 
        "-rpe_pair", str(rpe_image_path.resolve()), 
        "-eddy_fakeb", "0.85,1.2,1", 
        "-rpe_te", "60", 
        "-bshape", "1,0.6,1", 
        "-echo_time", te_args,
        "-normalize"
    ]

    designer = DesignerRunner(scratch_dir, cmd_config)
    tmi = TMIRunner.init_all_fa(scratch_dir)
    stats_computer = StatsComputer(designer, tmi, roi_dir=data_dir, valid_echo_times=valid_echo_times)
    runner = E2ERunner(designer, tmi, stats_computer, scratch_dir=scratch_dir)

    return runner


def prepare_heal_coronal_e2e_runner(scratch_dir: Path, data_dir: Path) -> E2ERunner:
    cmd_config = [
        "-denoise", 
        "-shrinkage", "frob", 
        "-adaptive_patch", 
        "-degibbs", 
        "-pf", "1", 
        "-normalize",
    ]

    designer = DesignerRunner(scratch_dir, cmd_config)
    tmi = TMIRunner.init_dti(scratch_dir, with_sigma=False, with_mask=True)
    stats_computer = StatsComputer(designer, tmi, roi_dir=data_dir, with_skull_stripping=True)
    runner = E2ERunner(designer, tmi, stats_computer, scratch_dir=scratch_dir)

    return runner

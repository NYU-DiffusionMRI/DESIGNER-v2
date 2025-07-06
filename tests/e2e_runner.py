from pathlib import Path
import shutil
import subprocess
from typing import List, Tuple, Optional, Dict
import json
from functools import reduce

from nibabel.nifti1 import Nifti1Image
import numpy as np

from tests.types import DWIStage, FAType, DWIImagePath, StatsDict
from tests.utils import compute_roi_mean_and_std, create_binary_mask_from_fa, extract_mean_b0


class DesignerRunner:

    def __init__(self, scratch_dir: Path, cli_options: List[str]):
        self.scratch_dir = scratch_dir
        self.processing_dir = self.scratch_dir / "processing"
        self.cli_options = cli_options

        self.out_path = DWIImagePath(
            nifti=self.scratch_dir / "dwi_designer.nii",
            bval=self.scratch_dir / "dwi_designer.bval",
            bvec=self.scratch_dir / "dwi_designer.bvec"
        )


    def run(self,  input_image_paths: List[Path]) -> Tuple[DWIImagePath, Path, Path]:
        dwi_args=",".join([str(path) for path in input_image_paths])

        if self.processing_dir.exists():
            shutil.rmtree(self.processing_dir)

        ret = subprocess.run([
            "designer", 
            *self.cli_options,
            "-mask", 
            "-scratch", str(self.processing_dir), 
            "-nocleanup", 
            dwi_args, 
            str(self.out_path.nifti)
        ])
        assert ret.returncode == 0
        assert self.out_path.nifti.exists()

        sigma_path = self.processing_dir / "sigma.nii"
        brain_mask_path = self.processing_dir / "brain_mask.nii"

        return self.out_path, sigma_path, brain_mask_path


    def get_dwi_path(self, stage: DWIStage) -> DWIImagePath:
        if stage == "designer":
            return self.out_path
        elif stage == "topup":
            return DWIImagePath(
                nifti=self.processing_dir / "eddy_processing" / "dwi_pe_0_applytopup.nii.gz",
                bval=self.scratch_dir / "dwi_designer.bval",    # same as designer bval
                bvec=self.scratch_dir / "dwi_designer.bvec"
            )
        else:
            # For intermediate stages, files are in processing directory
            stage_to_base = {
                "denoising": "dwidn",
                "degibbs": "working_rpg",
                "eddy": "dwiec",
                "b1correct": "dwibc",
                "rician": "dwirc"
            }
            base_name = stage_to_base[stage]
            
            return DWIImagePath(
                nifti=self.processing_dir / f"{base_name}.nii",
                bval=self.processing_dir / f"{base_name}.bval",
                bvec=self.processing_dir / f"{base_name}.bvec"
            )


class TMIRunner:

    # TODO: improve removing boolean flags
    def __init__(self, scratch_dir: Path, cli_options: List[str], *, with_sigma: bool = True, with_mask: bool = True):
        self.params_dir = scratch_dir / "params"
        self.cli_options = cli_options

        self.with_sigma = with_sigma
        self.with_mask = with_mask
    

    def run(self, input_image_path: Path, sigma_path: Path, brain_mask_path: Path):
        if self.params_dir.exists():
            shutil.rmtree(self.params_dir)

        cmd = ["tmi", *self.cli_options]
        # if sigma_path is not None:
        if self.with_sigma:
            cmd.extend(["-sigma", str(sigma_path)])
        # if brain_mask_path is not None:
        if self.with_mask:
            cmd.extend(["-mask", str(brain_mask_path)])
        cmd.extend([str(input_image_path), str(self.params_dir)])

        ret = subprocess.run(cmd)
        assert ret.returncode == 0
        assert self.params_dir.exists()


    def get_fa_path(self, fa_type: FAType, echo_time: Optional[float] = None) -> Path:
        if echo_time is None:
            return self.params_dir / f"fa_{fa_type}.nii"
        else:
            return self.params_dir / f"fa_{fa_type}_te{echo_time}.nii"


    # TODO: improve
    @classmethod
    def init_all_fa(cls, scratch_dir: Path, *, with_sigma: bool = True, with_mask: bool = True):
        return cls(scratch_dir, ["-SMI", "-DKI", "-WDKI", "-DTI"], with_sigma=with_sigma, with_mask=with_mask)


    # TODO: improve
    @classmethod
    def init_dti(cls, scratch_dir: Path, *, with_sigma: bool = True, with_mask: bool = True):
        return cls(scratch_dir, ["-DTI"], with_sigma=with_sigma, with_mask=with_mask)    


class StatsComputer:

    def __init__(self, designer: DesignerRunner, tmi: TMIRunner, roi_dir: Path, *, valid_echo_times: Optional[List[float]] = None, with_skull_stripping: bool = False):
        self.designer = designer
        self.tmi = tmi
        self.valid_echo_times = valid_echo_times
        self.with_skull_stripping = with_skull_stripping
        
        self.roi1 = Nifti1Image.from_filename(roi_dir / "roi1.nii.gz")
        self.roi2 = Nifti1Image.from_filename(roi_dir / "roi2.nii.gz")
        self.single_voxel = Nifti1Image.from_filename(roi_dir / "voxel.nii.gz")
        
        self.wm_roi: Optional[Nifti1Image] = None
        self.echo_to_wm_roi: Optional[Dict[float, Nifti1Image]] = None


    def init_wm_roi(self):
        """Initialize white matter ROI(s) after FA images have been generated."""
        if self.valid_echo_times is None and not self.with_skull_stripping:
            self.wm_roi = create_binary_mask_from_fa(self.tmi.get_fa_path("dki"))
        elif self.valid_echo_times is not None and not self.with_skull_stripping:
            self.echo_to_wm_roi = {te: create_binary_mask_from_fa(self.tmi.get_fa_path("dki", te)) for te in self.valid_echo_times}
            merged_wm_roi = reduce(np.bitwise_or, [roi.get_fdata().astype(bool) for roi in self.echo_to_wm_roi.values()])
            first_echo_time = self.valid_echo_times[0]
            self.wm_roi = Nifti1Image(merged_wm_roi, self.echo_to_wm_roi[first_echo_time].affine, self.echo_to_wm_roi[first_echo_time].header)
        elif self.valid_echo_times is None and self.with_skull_stripping:
            pass
        else:
            raise ValueError("Invalid combination of echo times and skull stripping. Currently not supported.")


    def count_white_matter_voxels(self) -> int:
        if self.wm_roi is None:
            raise RuntimeError("white matter ROI has not been initialized. Call init_wm_roi() first.")
        return np.count_nonzero(self.wm_roi.get_fdata())


    def compute_b0_roi_stats(self, stage: DWIStage) -> StatsDict:
        if self.wm_roi is None:
            raise RuntimeError("white matter ROI has not been initialized. Call init_wm_roi() first.")
            
        image_path = self.designer.get_dwi_path(stage)
        b0_image = extract_mean_b0(image_path.nifti, image_path.bval)

        wm_mean, wm_std = compute_roi_mean_and_std(b0_image, self.wm_roi)
        roi1_mean, roi1_std = compute_roi_mean_and_std(b0_image, self.roi1)
        roi2_mean, roi2_std = compute_roi_mean_and_std(b0_image, self.roi2)
        voxel_mean, _ = compute_roi_mean_and_std(b0_image, self.single_voxel)

        return {
            "wm": [wm_mean, wm_std],
            "roi1": [roi1_mean, roi1_std],
            "roi2": [roi2_mean, roi2_std],
            "voxel": [voxel_mean]
        }        
    

    def compute_fa_roi_stats(self, fa_type: FAType, echo_time: Optional[float] = None) -> StatsDict:
        if self.wm_roi is None:
            raise RuntimeError("white matter ROI has not been initialized. Call init_wm_roi() first.")
            
        image_path = self.tmi.get_fa_path(fa_type, echo_time)
        fa_image = Nifti1Image.from_filename(image_path)
        fa_data = fa_image.get_fdata()
        fa_data_no_nan = np.nan_to_num(fa_data, nan=0.0, posinf=0.0, neginf=0.0)
        fa_image_no_nan = Nifti1Image(fa_data_no_nan, fa_image.affine, fa_image.header)

        if echo_time is None:
            wm_roi_for_stats = self.wm_roi
        else:
            # self.echo_to_wm_roi is not None if self.wm_roi is not None
            assert self.echo_to_wm_roi is not None
            wm_roi_for_stats = self.echo_to_wm_roi[echo_time]

        wm_mean, wm_std = compute_roi_mean_and_std(fa_image_no_nan, wm_roi_for_stats)
        
        roi1_mean, roi1_std = compute_roi_mean_and_std(fa_image_no_nan, self.roi1)
        roi2_mean, roi2_std = compute_roi_mean_and_std(fa_image_no_nan, self.roi2)
        voxel_mean, _ = compute_roi_mean_and_std(fa_image_no_nan, self.single_voxel)

        return {
            "wm": [wm_mean, wm_std],
            "roi1": [roi1_mean, roi1_std],
            "roi2": [roi2_mean, roi2_std],
            "voxel": [voxel_mean]
        } 

    
    def save_benchmark(self, save_path: Path, stages: List[DWIStage], fa_types: List[FAType], valid_echo_times: Optional[List[float]] = None):
        if self.wm_roi is None:
            raise RuntimeError("white matter ROI has not been initialized. Call init_wm_roi() first.")
            
        b0_stats: Dict[DWIStage, StatsDict] = {}
        fa_stats: Dict[str, StatsDict] = {}

        wm_voxel_cnt = self.count_white_matter_voxels()

        for stage in stages:
            stats = self.compute_b0_roi_stats(stage)
            b0_stats[stage] = stats
        
        for fa_type in fa_types:
            if valid_echo_times is None:
                stats = self.compute_fa_roi_stats(fa_type)
                fa_stats[fa_type] = stats
            else:
                for echo_time in valid_echo_times:
                    stats = self.compute_fa_roi_stats(fa_type, echo_time)
                    fa_stats[f"{fa_type}_te{echo_time}"] = stats

        benchmark_data = {
            "white_matter_voxel_count": wm_voxel_cnt,
            "b0_stats": b0_stats,
            "fa_stats": fa_stats
        }

        with open(save_path, "w") as f:
            json.dump(benchmark_data, f, indent=2)


class E2ERunner:

    def __init__(self, designer: DesignerRunner, tmi: TMIRunner, stats_computer: StatsComputer, *, scratch_dir: Path):
        self.designer = designer
        self.tmi = tmi
        self.stats_computer = stats_computer

        self.scratch_dir = scratch_dir


    def run(self, input_image_paths: List[Path]):
        dwi_designer_path, sigma_path, brain_mask_path = self.designer.run(input_image_paths)

        self.tmi.run(dwi_designer_path.nifti, sigma_path, brain_mask_path)
        
        # Initialize ROIs after FA images have been generated
        self.stats_computer.init_wm_roi()


    def cleanup(self):
        if self.scratch_dir.exists():
            shutil.rmtree(self.scratch_dir)


    def __getattr__(self, name: str):
        """Forward any unknown attribute/method access to stats_computer instance."""
        return getattr(self.stats_computer, name)


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
    e2e_runner = E2ERunner(designer, tmi, stats_computer, scratch_dir=scratch_dir)

    return e2e_runner


def prepare_complex_data_e2e_runner(scratch_dir: Path, data_dir: Path) -> E2ERunner:
    phase_image_paths=[
        data_dir / "phase1_ds.nii.gz", 
        data_dir / "phase2_ds.nii.gz", 
        data_dir / "phase3_ds.nii.gz"
    ]

    rpe_image_path=data_dir / "rpe_ds.nii.gz"

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
    stats_computer = StatsComputer(designer, tmi, data_dir, valid_echo_times=valid_echo_times)
    e2e_runner = E2ERunner(designer, tmi, stats_computer, scratch_dir=scratch_dir)

    return e2e_runner

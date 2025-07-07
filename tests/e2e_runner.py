from pathlib import Path
import shutil
import subprocess
from functools import reduce
import json
from typing import List, Tuple, Optional, Dict

import numpy as np
import nibabel as nib
from nibabel.nifti1 import Nifti1Image
from nipreps.synthstrip.wrappers.nipype import SynthStrip

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
        self.scratch_dir = scratch_dir
        self.params_dir = self.scratch_dir / "params"
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


    DIFFUSION_MODELS = {
        "SMI": "-SMI",    # Simple Model of diffusion Imaging
        "DKI": "-DKI",    # Diffusion Kurtosis Imaging
        "WDKI": "-WDKI",  # White matter Diffusion Kurtosis Imaging
        "DTI": "-DTI"     # Diffusion Tensor Imaging
    }


    @classmethod
    def create_with_all_models(cls, scratch_dir: Path, *, with_sigma: bool = True, with_mask: bool = True) -> "TMIRunner":
        return cls(scratch_dir, list(cls.DIFFUSION_MODELS.values()), with_sigma=with_sigma, with_mask=with_mask)


    @classmethod
    def create_dti_only(cls, scratch_dir: Path, *, with_sigma: bool = True, with_mask: bool = True) -> "TMIRunner":
        return cls(scratch_dir, [cls.DIFFUSION_MODELS["DTI"]], with_sigma=with_sigma, with_mask=with_mask)


class StatsComputer:

    def __init__(
        self, 
        designer: DesignerRunner, 
        tmi: TMIRunner, 
        roi_dir: Path, 
        *, 
        valid_echo_times: Optional[List[float]] = None, 
        with_skull_stripping: bool = False
    ):  
        self.designer = designer
        self.tmi = tmi
        self.valid_echo_times = valid_echo_times
        self.with_skull_stripping = with_skull_stripping

        self.processing_dir = self.designer.processing_dir
        self.params_dir = self.tmi.params_dir
        self.scratch_dir = self.designer.scratch_dir
        
        self.roi1 = Nifti1Image.from_filename(roi_dir / "roi1.nii.gz")
        self.roi2 = Nifti1Image.from_filename(roi_dir / "roi2.nii.gz")
        self.single_voxel = Nifti1Image.from_filename(roi_dir / "voxel.nii.gz")
        
        self.wm_roi: Optional[Nifti1Image] = None
        self.echo_to_wm_roi: Optional[Dict[float, Nifti1Image]] = None


    def init_wm_roi(self):
        """Initialize white matter ROI(s) after FA images have been generated."""
        if self.valid_echo_times is None and not self.with_skull_stripping:
            self._init_single_wm_roi()
        elif self.valid_echo_times is not None and not self.with_skull_stripping:
            self._init_multi_echo_wm_roi()
        elif self.valid_echo_times is None and self.with_skull_stripping:
            self._init_skull_stripped_wm_roi()
        else:
            raise ValueError("Invalid combination of echo times and skull stripping. Currently not supported.")


    def _init_single_wm_roi(self):
        self.wm_roi = create_binary_mask_from_fa(self.tmi.get_fa_path("dki"))


    def _init_multi_echo_wm_roi(self):
        assert self.valid_echo_times is not None

        self.echo_to_wm_roi = {
                te: create_binary_mask_from_fa(self.tmi.get_fa_path("dki", te))
                for te in self.valid_echo_times
            }
        merged_wm_roi = reduce(
            np.bitwise_or,
            [roi.get_fdata().astype(bool) for roi in self.echo_to_wm_roi.values()]
        )
        first_echo_time = self.valid_echo_times[0]
        self.wm_roi = Nifti1Image(
            merged_wm_roi,
            self.echo_to_wm_roi[first_echo_time].affine,
            self.echo_to_wm_roi[first_echo_time].header
        )


    # TODO: improve (same threshold for creating wm roi is used in multiple places)
    def _init_skull_stripped_wm_roi(self):
        synth = SynthStrip()

        synth.inputs.in_file = str(self.processing_dir / "b0bc.nii")
        synth.inputs.model = "tests/models/synthstrip_v7.4.1_.1.pt"
        synth.inputs.out_file = str(self.scratch_dir / "brain.nii.gz")  # this is not used
        synth.inputs.out_mask = str(self.scratch_dir / "brain_mask.nii")
        synth.inputs.use_gpu = False
        result = synth.run()

        assert Path(result.outputs.out_mask).exists()

        dti_path = self.tmi.get_fa_path("dti")
        dti_image = Nifti1Image.from_filename(dti_path)
        brain_mask = Nifti1Image.from_filename(result.outputs.out_mask)

        fa_dti_can = nib.as_closest_canonical(dti_image)
        mask_can = nib.as_closest_canonical(brain_mask)

        brain = fa_dti_can.get_fdata() * mask_can.get_fdata()
        brain_no_nan = np.nan_to_num(brain, nan=0.0, posinf=0.0, neginf=0.0)

        threshold = 0.3
        wm_roi = (brain_no_nan >= threshold).astype(np.uint8)
        self.wm_roi = Nifti1Image(wm_roi, fa_dti_can.affine, fa_dti_can.header)


    def count_white_matter_voxels(self) -> int:
        if self.wm_roi is None:
            raise RuntimeError("white matter ROI has not been initialized. Call init_wm_roi() first.")
        return np.count_nonzero(self.wm_roi.get_fdata())


    def compute_b0_roi_stats(self, stage: DWIStage) -> StatsDict:
        if self.wm_roi is None:
            raise RuntimeError("white matter ROI has not been initialized. Call init_wm_roi() first.")
            
        image_path: DWIImagePath = self.designer.get_dwi_path(stage)
        b0_image: Nifti1Image = extract_mean_b0(image_path.nifti, image_path.bval)

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

    
    def save_benchmark(
        self, 
        save_path: Path, 
        stages: List[DWIStage], 
        fa_types: List[FAType], 
        valid_echo_times: Optional[List[float]] = None
    ):
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

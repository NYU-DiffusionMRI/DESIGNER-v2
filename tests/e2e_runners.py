from pathlib import Path
import shutil
import subprocess
from functools import reduce
import json
from typing import List, Tuple, Optional, Dict

import numpy as np
from nibabel.nifti1 import Nifti1Image
from nipreps.synthstrip.wrappers.nipype import SynthStrip

from tests.types import DWIStage, DiffusionModelType, DWIImagePath, StatsDict
from tests.utils import compute_roi_mean_and_std, create_binary_mask_from_fa, extract_mean_b0


class DesignerRunner:

    def __init__(self, scratch_dir: Path, input_dwi_paths: List[Path], cli_options: List[str]):
        self.scratch_dir = scratch_dir
        self.processing_dir = self.scratch_dir / "processing"
        self.input_dwi_paths = input_dwi_paths

        self.cli_options = cli_options

        self.out_path = DWIImagePath(
            nifti=self.scratch_dir / "dwi_designer.nii",
            bval=self.scratch_dir / "dwi_designer.bval",
            bvec=self.scratch_dir / "dwi_designer.bvec"
        )

        self.brain_mask_path: Path = self.processing_dir / "brain_mask.nii"


    def run(self) -> Tuple[DWIImagePath, Path, Path]:
        dwi_args=",".join([str(path) for path in self.input_dwi_paths])

        if self.processing_dir.exists():
            shutil.rmtree(self.processing_dir)

        ret = subprocess.run([
            "designer", 
            *self.cli_options,
            "-mask", 
            "-scratch", str(self.processing_dir), 
            "-nocleanup", 
            dwi_args, 
            str(self.out_path["nifti"])
        ])
        assert ret.returncode == 0
        assert self.out_path["nifti"].exists()

        sigma_path = self.processing_dir / "sigma.nii"

        return self.out_path, sigma_path, self.brain_mask_path


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
    def __init__(self, scratch_dir: Path, input_dwi_path: Path, diffusion_models: List[DiffusionModelType], *, with_sigma: bool = True, with_mask: bool = True):
        self.scratch_dir = scratch_dir
        self.params_dir = self.scratch_dir / "params"
        self.input_dwi_path = input_dwi_path

        self.diffusion_models = diffusion_models

        self.with_sigma = with_sigma
        self.with_mask = with_mask
    

    def run(self, sigma_path: Path, brain_mask_path: Path):
        if self.params_dir.exists():
            shutil.rmtree(self.params_dir)

        diffusion_args = [f"-{model.upper()}" for model in self.diffusion_models]
        cmd = ["tmi", *diffusion_args]

        # if sigma_path is not None:
        if self.with_sigma:
            cmd.extend(["-sigma", str(sigma_path)])
        # if brain_mask_path is not None:
        if self.with_mask:
            cmd.extend(["-mask", str(brain_mask_path)])
        cmd.extend([str(self.input_dwi_path), str(self.params_dir)])

        ret = subprocess.run(cmd)
        assert ret.returncode == 0
        assert self.params_dir.exists()


    def get_fa_path(self, fa_model: DiffusionModelType, echo_time: Optional[float] = None) -> Path:
        if echo_time is None:
            return self.params_dir / f"fa_{fa_model}.nii"
        else:
            return self.params_dir / f"fa_{fa_model}_te{echo_time}.nii"


class E2ERunner:

    def __init__(
        self,
        scratch_dir: Path,
        input_dwi_paths: List[Path],
        designer_cli_options: List[str],
        diffusion_models: List[DiffusionModelType],
        *,
        with_sigma: bool = True,
        with_mask: bool = True
    ):
        self.scratch_dir = scratch_dir

        self.designer = DesignerRunner(scratch_dir, input_dwi_paths, designer_cli_options)

        tmi_input_path = self.designer.get_dwi_path("designer")["nifti"]
        self.tmi = TMIRunner(scratch_dir, tmi_input_path, diffusion_models, with_sigma=with_sigma, with_mask=with_mask)

        self.processing_dir = self.designer.processing_dir
        self.params_dir = self.tmi.params_dir

        self.brain_mask_path: Path = self.designer.brain_mask_path


    def run(self):
        _, sigma_path, brain_mask_path = self.designer.run()

        self.tmi.run(sigma_path, brain_mask_path)
        

    def cleanup(self):
        if self.scratch_dir.exists():
            shutil.rmtree(self.scratch_dir)


    def get_dwi_path(self, stage: DWIStage) -> DWIImagePath:
        return self.designer.get_dwi_path(stage)


    def get_fa_path(self, fa_model: DiffusionModelType, echo_time: Optional[float] = None) -> Path:
        return self.tmi.get_fa_path(fa_model, echo_time)


class StatsComputer:

    def __init__(
        self, 
        e2e_runner: E2ERunner,
        roi_dir: Path, 
        *, 
        valid_echo_times: Optional[List[float]] = None, 
        with_skull_stripping: bool = False
    ):  
        self.e2e_runner = e2e_runner
        self.valid_echo_times = valid_echo_times
        self.with_skull_stripping = with_skull_stripping

        self.processing_dir = self.e2e_runner.processing_dir
        self.params_dir = self.e2e_runner.params_dir
        self.scratch_dir = self.e2e_runner.scratch_dir
        
        self.roi1 = Nifti1Image.from_filename(roi_dir / "roi1.nii.gz")
        self.roi2 = Nifti1Image.from_filename(roi_dir / "roi2.nii.gz")
        self.single_voxel = Nifti1Image.from_filename(roi_dir / "voxel.nii.gz")
        
        self.echo_to_wm_roi: Optional[Dict[float, Nifti1Image]] = None
        self.brain_mask_path: Path = self.e2e_runner.brain_mask_path

        if self.valid_echo_times is None and not self.with_skull_stripping:
            self.wm_roi = self._generate_single_wm_roi()
        elif self.valid_echo_times is not None and not self.with_skull_stripping:
            self.wm_roi, self.echo_to_wm_roi = self._generate_multi_echo_wm_roi()
        elif self.valid_echo_times is None and self.with_skull_stripping:
            # using the brain mask from skull stripping
            self.wm_roi, self.brain_mask_path = self._generate_skull_stripped_wm_roi()
        else:
            raise ValueError("Invalid combination of echo times and skull stripping. Currently not supported.")


    def _generate_single_wm_roi(self) -> Nifti1Image:
        return create_binary_mask_from_fa(self.e2e_runner.get_fa_path("dki"))


    def _generate_multi_echo_wm_roi(self) -> Tuple[Nifti1Image, Dict[float, Nifti1Image]]:
        assert self.valid_echo_times is not None

        echo_to_wm_roi = {
                te: create_binary_mask_from_fa(self.e2e_runner.get_fa_path("dki", te))
                for te in self.valid_echo_times
            }
        merged_wm_roi = reduce(
            np.bitwise_or,
            [roi.get_fdata().astype(bool) for roi in echo_to_wm_roi.values()]
        )
        first_echo_time = self.valid_echo_times[0]
        merged_wm_image = Nifti1Image(
            merged_wm_roi,
            echo_to_wm_roi[first_echo_time].affine,
            echo_to_wm_roi[first_echo_time].header
        )

        return merged_wm_image, echo_to_wm_roi


    def _generate_skull_stripped_wm_roi(self) -> Tuple[Nifti1Image, Path]:
        synth = SynthStrip()

        synth.inputs.in_file = str(self.processing_dir / "b0bc.nii")
        synth.inputs.model = "tests/models/synthstrip_v7.3.2.pt"
        synth.inputs.out_file = str(self.scratch_dir / "brain.nii.gz")  # this is not used
        synth.inputs.out_mask = str(self.scratch_dir / "brain_mask.nii")
        synth.inputs.use_gpu = False
        result = synth.run()

        mask_path = Path(result.outputs.out_mask)
        assert mask_path.exists()

        fa_path = self.e2e_runner.get_fa_path("dti")

        return create_binary_mask_from_fa(fa_path, brain_mask=mask_path), mask_path


    def compute_wm_ratio(self) -> float:
        wm_voxel_cnt = np.count_nonzero(self.wm_roi.get_fdata())
        brain_mask = Nifti1Image.from_filename(self.brain_mask_path)
        total_voxel_cnt = np.count_nonzero(brain_mask.get_fdata())

        return wm_voxel_cnt / total_voxel_cnt


    def compute_b0_roi_stats(self, stage: DWIStage) -> StatsDict:
        image_path: DWIImagePath = self.e2e_runner.get_dwi_path(stage)
        b0_image: Nifti1Image = extract_mean_b0(image_path["nifti"], image_path["bval"])

        wm_mean, wm_std = compute_roi_mean_and_std(b0_image, self.wm_roi)
        roi1_mean, roi1_std = compute_roi_mean_and_std(b0_image, self.roi1)
        roi2_mean, roi2_std = compute_roi_mean_and_std(b0_image, self.roi2)
        voxel_mean, _ = compute_roi_mean_and_std(b0_image, self.single_voxel)

        return StatsDict(
            wm=[wm_mean, wm_std],
            roi1=[roi1_mean, roi1_std],
            roi2=[roi2_mean, roi2_std],
            voxel=[voxel_mean]
        )


    def compute_fa_roi_stats(self, fa_model: DiffusionModelType, echo_time: Optional[float] = None) -> StatsDict:
        image_path = self.e2e_runner.get_fa_path(fa_model, echo_time)
        fa_image = Nifti1Image.from_filename(image_path)
        fa_data = fa_image.get_fdata()
        fa_data_no_nan = np.nan_to_num(fa_data, nan=0.0, posinf=0.0, neginf=0.0)
        fa_image_no_nan = Nifti1Image(fa_data_no_nan, fa_image.affine, fa_image.header)

        if echo_time is None:
            wm_roi_for_stats = self.wm_roi
        else:
            assert self.valid_echo_times is not None
            assert self.echo_to_wm_roi is not None
            wm_roi_for_stats = self.echo_to_wm_roi[echo_time]

        wm_mean, wm_std = compute_roi_mean_and_std(fa_image_no_nan, wm_roi_for_stats)
        roi1_mean, roi1_std = compute_roi_mean_and_std(fa_image_no_nan, self.roi1)
        roi2_mean, roi2_std = compute_roi_mean_and_std(fa_image_no_nan, self.roi2)
        voxel_mean, _ = compute_roi_mean_and_std(fa_image_no_nan, self.single_voxel)

        return StatsDict(
            wm=[wm_mean, wm_std],
            roi1=[roi1_mean, roi1_std],
            roi2=[roi2_mean, roi2_std],
            voxel=[voxel_mean]
        )

    
    def save_benchmark(
        self, 
        save_path: Path, 
        stages: List[DWIStage], 
        fa_models: List[DiffusionModelType], 
        valid_echo_times: Optional[List[float]] = None
    ):
        b0_stats: Dict[DWIStage, StatsDict] = {}
        fa_stats: Dict[str, StatsDict] = {}

        wm_ratio = self.compute_wm_ratio()

        for stage in stages:
            stats = self.compute_b0_roi_stats(stage)
            b0_stats[stage] = stats
        
        for fa_model in fa_models:
            if valid_echo_times is None:
                stats = self.compute_fa_roi_stats(fa_model)
                fa_stats[fa_model] = stats
            else:
                for echo_time in valid_echo_times:
                    stats = self.compute_fa_roi_stats(fa_model, echo_time)
                    fa_stats[f"{fa_model}_te{echo_time}"] = stats

        benchmark_data = {
            "wm_ratio": wm_ratio,
            "b0_stats": b0_stats,
            "fa_stats": fa_stats
        }

        with open(save_path, "w") as f:
            json.dump(benchmark_data, f, indent=2)

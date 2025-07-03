from pathlib import Path
import shutil
import subprocess
import json
from typing import List, Tuple, Optional, Dict

import numpy as np
from nibabel.nifti1 import Nifti1Image

from tests.utils import compute_roi_mean_and_std, create_binary_mask_from_fa, extract_mean_b0
from tests.types import DWIStage, FAType, DWIImagePath, StatsDict


class DesignerRunner:

    def __init__(self, tmp_dir: Path, cli_options: List[str]):
        self.tmp_dir = tmp_dir
        self.processing_dir = self.tmp_dir / "processing"
        self.cli_options = cli_options

        self.out_path = DWIImagePath(
            nifti=self.tmp_dir / "dwi_designer.nii",
            bval=self.tmp_dir / "dwi_designer.bval",
            bvec=self.tmp_dir / "dwi_designer.bvec"
        )


    def run(self, input_image_paths: List[Path]) -> Tuple[DWIImagePath, Path, Path]:
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
                bval=self.tmp_dir / "dwi_designer.bval",    # same as designer bval
                bvec=self.tmp_dir / "dwi_designer.bvec"
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


    @classmethod
    def init_meso_nonsquare(cls, tmp_dir: Path):
        return cls(tmp_dir, [
            "-denoise", 
            "-shrinkage", "frob", 
            "-adaptive_patch", 
            "-rician", 
            "-degibbs", 
            "-b1correct", 
            "-normalize"
        ])


    @classmethod
    def init_meso_nonsquare_wo_bids(cls, tmp_dir: Path):
        return cls(tmp_dir, [
            "-denoise", 
            "-shrinkage", "frob", 
            "-adaptive_patch", 
            "-rician",
            "-pf", "6/8",
            "-pe_dir", "AP",
            "-degibbs", 
            "-b1correct", 
            "-normalize",
        ])


class TMIRunner:

    def __init__(self, tmp_dir: Path, cli_options: List[str], *, with_sigma: bool = True, with_mask: bool = True):
        self.params_dir = tmp_dir / "params"
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
            return self.params_dir / f"fa_{fa_type}_{echo_time}.nii"


    @classmethod
    def init_all_fa(cls, tmp_dir: Path, *, with_sigma: bool = True, with_mask: bool = True):
        return cls(tmp_dir, ["-SMI", "-DKI", "-WDKI", "-DTI"], with_sigma=with_sigma, with_mask=with_mask)


    @classmethod
    def init_dti(cls, tmp_dir: Path, *, with_sigma: bool = True, with_mask: bool = True):
        return cls(tmp_dir, ["-DTI"], with_sigma=with_sigma, with_mask=with_mask)    


class MESONonsquareRunner:

    def __init__(self, designer: DesignerRunner, tmi: TMIRunner, tmp_dir: Path):
        self.data_dir = Path("tests/data/D1/")
        self.tmp_dir = tmp_dir

        self.designer = designer
        self.tmi = tmi

        self._white_matter_roi = None
        self.roi1 = Nifti1Image.from_filename(self.data_dir / "roi1.nii.gz")
        self.roi2 = Nifti1Image.from_filename(self.data_dir / "roi2.nii.gz")
        self.single_voxel = Nifti1Image.from_filename(self.data_dir / "voxel.nii.gz")


    @property
    def white_matter_roi(self) -> Nifti1Image:
        if self._white_matter_roi is None:
            self._white_matter_roi = create_binary_mask_from_fa(self.tmi.get_fa_path("dki"))

        return self._white_matter_roi


    def run_pipeline(self):
        dwi_designer_path, sigma_path, brain_mask_path = self.designer.run([
            self.data_dir / "meso_slice_crop.nii.gz", 
            self.data_dir / "research_slice_crop.nii.gz"
        ])

        self.tmi.run(dwi_designer_path.nifti, sigma_path, brain_mask_path)


    def cleanup(self):
        if self.tmp_dir.exists():
            shutil.rmtree(self.tmp_dir)


    def count_white_matter_voxels(self) -> int:
        return np.count_nonzero(self.white_matter_roi.get_fdata())


    def compute_dwi_b0_stats(self, stage: DWIStage) -> StatsDict:
        image_path = self.designer.get_dwi_path(stage)
        b0_image = extract_mean_b0(image_path.nifti, image_path.bval)

        wm_mean, wm_std = compute_roi_mean_and_std(b0_image, self.white_matter_roi)
        roi1_mean, roi1_std = compute_roi_mean_and_std(b0_image, self.roi1)
        roi2_mean, roi2_std = compute_roi_mean_and_std(b0_image, self.roi2)
        voxel_mean, _ = compute_roi_mean_and_std(b0_image, self.single_voxel)
        
        return {
            "wm": [wm_mean, wm_std],
            "roi1": [roi1_mean, roi1_std],
            "roi2": [roi2_mean, roi2_std],
            "voxel": [voxel_mean]
        }


    def compute_fa_stats(self, fa_type: FAType) -> StatsDict:
        image_path = self.tmi.get_fa_path(fa_type)
        fa_image = Nifti1Image.from_filename(image_path)
        fa_data = fa_image.get_fdata()
        fa_data_no_nan = np.nan_to_num(fa_data, nan=0.0, posinf=0.0, neginf=0.0)
        fa_image_no_nan = Nifti1Image(fa_data_no_nan, fa_image.affine, fa_image.header)

        wm_mean, wm_std = compute_roi_mean_and_std(fa_image_no_nan, self.white_matter_roi)
        roi1_mean, roi1_std = compute_roi_mean_and_std(fa_image_no_nan, self.roi1)
        roi2_mean, roi2_std = compute_roi_mean_and_std(fa_image_no_nan, self.roi2)
        voxel_mean, _ = compute_roi_mean_and_std(fa_image_no_nan, self.single_voxel)

        return {
            "wm": [wm_mean, wm_std],
            "roi1": [roi1_mean, roi1_std],
            "roi2": [roi2_mean, roi2_std],
            "voxel": [voxel_mean]
        }


    def save_benchmark(self, save_path: Path, stages: List[DWIStage], fa_types: List[FAType]):
        b0_stats: Dict[DWIStage, StatsDict] = {}
        fa_stats: Dict[FAType, StatsDict] = {}

        wm_voxel_cnt = self.count_white_matter_voxels()

        for stage in stages:
            stats = self.compute_dwi_b0_stats(stage)
            b0_stats[stage] = stats

        for fa_type in fa_types:
            stats = self.compute_fa_stats(fa_type)
            fa_stats[fa_type] = stats

        benchmark_data = {
            "white_matter_voxel_count": wm_voxel_cnt,
            "b0_stats": b0_stats,
            "fa_stats": fa_stats
        }

        with open(save_path, "w") as f:
            json.dump(benchmark_data, f, indent=2)

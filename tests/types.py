from typing import Literal, List, TypedDict, get_args, NotRequired
from pathlib import Path


DWIStage = Literal["denoising", "degibbs","topup", "eddy", "b1correct", "rician", "designer"]
DiffusionModelType = Literal["smi", "dti", "dki", "wdki"]


def is_valid_dwi_stage(dwi_stage: str) -> bool:
    return dwi_stage in get_args(DWIStage)


def is_valid_diffusion_model_type(diffusion_model_type: str) -> bool:
    return diffusion_model_type in get_args(DiffusionModelType)


class DWIImagePath(TypedDict):
    nifti: Path
    bval: Path
    bvec: Path


class StatsDict(TypedDict):
    wm: List[float]
    roi1: List[float]
    roi2: List[float]
    voxel: List[float]


class ToleranceProfile(TypedDict):
    wm_ratio_atol: float
    b0_mean_rtol: float
    b0_std_rtol: float
    b0_wm_mean_rtol: NotRequired[float]
    b0_wm_std_rtol: NotRequired[float]
    fa_mean_atol: float
    fa_std_atol: float
    fa_wm_mean_rtol: NotRequired[float]
    fa_wm_std_rtol: NotRequired[float]

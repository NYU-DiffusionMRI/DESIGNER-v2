from typing import Literal, Dict, List, TypedDict, get_args, NotRequired
from pathlib import Path
from dataclasses import dataclass


DWIStage = Literal["denoising", "degibbs","topup", "eddy", "b1correct", "rician", "designer"]
DiffusionModelType = Literal["smi", "dti", "dki", "wdki"]


def is_valid_dwi_stage(dwi_stage: str) -> bool:
    return dwi_stage in get_args(DWIStage)


def is_valid_diffusion_model_type(diffusion_model_type: str) -> bool:
    return diffusion_model_type in get_args(DiffusionModelType)


@dataclass
class DWIImagePath:
    nifti: Path
    bval: Path
    bvec: Path


StatsDict = Dict[str, List[float]]


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

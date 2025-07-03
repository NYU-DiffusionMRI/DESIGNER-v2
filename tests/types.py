from typing import Literal, Dict, List
from pathlib import Path
from dataclasses import dataclass


DWIStage = Literal["denoising", "degibbs","topup", "eddy", "b1correct", "rician", "designer"]
FAType = Literal["dti", "dki", "wdki"]


@dataclass
class DWIImagePath:
    nifti: Path
    bval: Path
    bvec: Path

StatsDict = Dict[str, List[float]]

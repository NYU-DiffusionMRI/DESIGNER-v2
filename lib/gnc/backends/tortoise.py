"""GNC backend: TORTOISE's CreateGradientNonlinearityBMatrix."""

from pathlib import Path

import nibabel as nib
import numpy as np

from lib.designer_input_utils import splitext_
from lib.utils import force_remove


def resolve_is_ge(scanner: str) -> bool:
    """Map -gnc_scanner to TORTOISE's GE/non-GE distinction (--isGE flag).

    A plain equality check is safe here because -gnc_scanner is constrained
    to {'siemens', 'ge'} by argparse (designer2/designer.py), matching
    gradunwarp's own hard vendor constraint - no third value can reach here.
    """
    return scanner == 'ge'


def compute_warp(
    ref_volume_nii: str | Path,
    coeff_file: str | Path,
    scanner: str,
    out_prefix: str,
) -> str:
    """Runs TORTOISE's CreateGradientNonlinearityBMatrix on the reference
    volume - a single call that is a strict superset of
    CreateNonlinearityDisplacementMap, writing both an ITK displacement field
    and the per-voxel L tensor - then converts the displacement field to
    FSL's relative-warp convention (see _itk_field_to_fsl_rel_warp). The L
    tensor is deliberately left on disk for derive_graddev() to pick up right
    after this returns."""
    from mrtrix3 import run

    isge_flag = ' --isGE' if resolve_is_ge(scanner) else ''
    run.command(f'CreateGradientNonlinearityBMatrix -f {ref_volume_nii} -g {coeff_file}{isge_flag}')

    # Output names derive from ref_volume_nii's own basename, not out_prefix.
    ref_stem, _ = splitext_(str(ref_volume_nii))
    itk_field = f'{ref_stem}_gradwarp_field.nii'

    rel_warp = f'{out_prefix}_warp_rel.nii.gz'
    _itk_field_to_fsl_rel_warp(itk_field, ref_volume_nii, rel_warp)
    force_remove(itk_field)
    return rel_warp


def derive_graddev(
    ref_volume_nii: str | Path,
    coeff_file: str | Path,
    scanner: str,
    out_prefix: str,
) -> str:
    """Byproduct of compute_warp()'s single CreateGradientNonlinearityBMatrix
    call, left on disk for this."""
    ref_stem, _ = splitext_(str(ref_volume_nii))
    return f'{ref_stem}_graddev_c.nii'


def _lps_disp_to_fsl_rel(disp_lps, affine):
    """
    Convert per-voxel displacements (mm, components in the ITK LPS world
    frame, as TORTOISE writes them) into an FSL relative warp on the same grid
    (mm along the image's storage axes, x negated when det(affine) > 0, per
    FSL's internal radiological convention). Exact for oblique affines.

    disp_lps: (..., 3) array. affine: 4x4 voxel-to-RAS-world matrix.
    """
    M = np.asarray(affine, dtype=float)[:3, :3]
    d_ras = disp_lps * np.array([-1.0, -1.0, 1.0])
    d_vox = d_ras @ np.linalg.inv(M).T
    fsl = d_vox * np.linalg.norm(M, axis=0)
    if np.linalg.det(M) > 0:
        fsl[..., 0] *= -1
    return fsl


def _itk_field_to_fsl_rel_warp(
    itk_field_path: str | Path,
    ref_volume_nii: str | Path,
    out_rel_warp_path: str | Path,
) -> None:
    """
    TORTOISE's _gradwarp_field.nii is a pull displacement field (not absolute
    coordinates): corrected(r) = distorted(r + d(r)), d in mm in LPS
    (mk_displacement_siemens/_nonsiemens in TORTOISEV4's
    src/tools/gradnonlin/mk_displacementMaps.h). Re-express it as an FSL
    relative warp and canonicalize it through convertwarp --rel.
    """
    from mrtrix3 import MRtrixError, run

    ref_img = nib.load(str(ref_volume_nii))
    field_img = nib.load(str(itk_field_path))
    if field_img.shape[:3] != ref_img.shape[:3]:
        raise MRtrixError(f'TORTOISE displacement field grid {field_img.shape[:3]} does not '
                          f'match reference volume grid {ref_img.shape[:3]}')
    disp = np.asarray(field_img.dataobj, dtype=np.float64).reshape(*ref_img.shape[:3], 3)
    rel = _lps_disp_to_fsl_rel(disp, ref_img.affine).astype(np.float32)

    tmp_path = f'{Path(str(out_rel_warp_path)).name.removesuffix(".nii.gz")}_from_itk_tmp.nii.gz'
    nib.save(nib.Nifti1Image(rel, ref_img.affine, ref_img.header), tmp_path)
    run.command(
        f'convertwarp --rel --ref={ref_volume_nii} --warp1={tmp_path} '
        f'--relout --out={out_rel_warp_path}'
    )
    force_remove(tmp_path)

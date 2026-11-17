"""GNC backend: WashU gradunwarp (gradient_unwarp.py) plus FSL's calc_grad_perc_dev."""

from pathlib import Path

from lib.utils import force_remove


def compute_warp(
    ref_volume_nii: str | Path,
    coeff_file: str | Path,
    scanner: str,
    out_prefix: str,
) -> str:
    """
    Run gradient_unwarp.py on the reference volume and convert its absolute
    warp field to FSL's relative-warp convention.

    Arguments as in GncBackend.compute_warp. gradient_unwarp.py's own
    outputs ({out_prefix}_direct.nii.gz and fullWarp_abs.nii.gz in the CWD)
    are removed; nothing is left behind for derive_graddev().

    Returns:
        Path to {out_prefix}_warp_rel.nii.gz
    """
    from mrtrix3 import run

    run.command(
        f'gradient_unwarp.py {ref_volume_nii} {out_prefix}_direct.nii.gz {scanner} '
        f'-g {coeff_file} -n --interp_order 3'
    )
    force_remove(f'{out_prefix}_direct.nii.gz')

    # gradient_unwarp.py always writes fullWarp_abs.nii.gz into the CWD.
    abs_warp = 'fullWarp_abs.nii.gz'
    rel_warp = f'{out_prefix}_warp_rel.nii.gz'
    run.command(
        f'convertwarp --abs --ref={ref_volume_nii} --warp1={abs_warp} '
        f'--relout --out={rel_warp}'
    )
    force_remove(abs_warp)
    return rel_warp


def derive_graddev(
    ref_volume_nii: str | Path,
    coeff_file: str | Path,
    scanner: str,
    out_prefix: str,
) -> str:
    """
    Derive the L tensor from gradunwarp's absolute warp via FSL's
    calc_grad_perc_dev (gradunwarp has no native L tensor).

    The warp is recomputed on the reference grid padded by one voxel and
    cropped back afterwards: calc_grad_perc_dev zero-pads outside the grid,
    so differentiating absolute coordinates on the unpadded grid corrupts
    the whole outermost voxel layer (values ~50-150 instead of ~1).

    Arguments as in GncBackend.derive_graddev. All intermediates
    ({out_prefix}_graddev_* and fullWarp_abs.nii.gz) are removed.

    Returns:
        Path to {out_prefix}_graddev_raw.nii.gz, same grid as ref_volume_nii
    """
    from mrtrix3 import run

    padded_ref = f'{out_prefix}_graddev_padref.nii.gz'
    run.command(f'mrgrid {ref_volume_nii} pad -uniform 1 {padded_ref} -force -quiet')
    run.command(
        f'gradient_unwarp.py {padded_ref} {out_prefix}_graddev_paddirect.nii.gz {scanner} '
        f'-g {coeff_file} -n --interp_order 3'
    )

    abs_warp = 'fullWarp_abs.nii.gz'
    perc_prefix = f'{out_prefix}_graddev_perc'
    run.command(f'calc_grad_perc_dev --fullwarp={abs_warp} --out={perc_prefix}')

    merged = f'{out_prefix}_graddev_merged.nii.gz'
    run.command(
        f'fslmerge -t {merged} {perc_prefix}_x.nii.gz {perc_prefix}_y.nii.gz {perc_prefix}_z.nii.gz'
    )

    # calc_grad_perc_dev reports 100 * d(warp_i)/d(x_j) with no identity
    # subtraction (so the diagonal reads ~100, not ~1, at undistorted
    # voxels); dividing by 100 reproduces TORTOISE's native L convention
    # (within 1e-3 on every voxel of real data).
    padded_raw = f'{out_prefix}_graddev_padraw.nii.gz'
    run.command(f'mrcalc {merged} 100 -div {padded_raw} -force -quiet')
    raw = f'{out_prefix}_graddev_raw.nii.gz'
    run.command(f'mrgrid {padded_raw} crop -uniform 1 {raw} -force -quiet')
    force_remove(
        padded_ref, f'{out_prefix}_graddev_paddirect.nii.gz', abs_warp, merged, padded_raw,
        f'{perc_prefix}_x.nii.gz', f'{perc_prefix}_y.nii.gz', f'{perc_prefix}_z.nii.gz',
    )
    return raw

"""Gradient nonlinearity correction (GNC): spatial unwarping via gradunwarp.

Computes a single geometric warp (per scan) from the scanner's gradient
coefficient file and applies it either as a standalone pre-eddy correction
(for topup/eddy to operate on geometrically-correct data) or composed with
eddy's own per-volume displacement fields so the combined GNC+eddy correction
is resampled in a single interpolation pass.
"""

from pathlib import Path
from typing import Any


def compute_gnc_warp(
    ref_volume_nii: str | Path,
    coeff_file: str | Path,
    scanner: str,
    out_prefix: str,
) -> str:
    """
    Compute the relative FSL warp for gradient nonlinearity correction.

    Runs gradient_unwarp.py on a single reference volume (the warp is static
    in scanner space, so one volume is enough) and converts its absolute
    warp field to FSL's relative-warp convention for use with applywarp.

    Args:
        ref_volume_nii: Path to the single reference volume to unwarp
        coeff_file: Path to the scanner gradient coefficient file
        scanner: Scanner vendor passed to gradient_unwarp.py (e.g. 'siemens')
        out_prefix: Basename prefix for the intermediate/output warp files

    Returns:
        Path to the relative warp file (FSL convention, for use with applywarp)
    """
    from mrtrix3 import run

    run.command(
        f'gradient_unwarp.py {ref_volume_nii} {out_prefix}_direct.nii.gz {scanner} '
        f'-g {coeff_file} -n --interp_order 3'
    )

    # gradient_unwarp.py always writes fullWarp_abs.nii.gz into the CWD
    abs_warp = 'fullWarp_abs.nii.gz'
    rel_warp = f'{out_prefix}_warp_rel.nii.gz'
    run.command(
        f'convertwarp --abs --ref={ref_volume_nii} --warp1={abs_warp} '
        f'--relout --out={rel_warp}'
    )
    return rel_warp


def apply_warp_to_series(
    input_path: str | Path,
    rel_warp_path: str | Path,
    ref_volume_nii: str | Path,
    output_path: str | Path,
    split_prefix: str,
) -> None:
    """
    Apply a shared relative warp to every volume of a 3D/4D series.

    No Jacobian modulation is applied (matches gradient_unwarp.py -n), so this
    is purely a geometric correction. Splits/recombines per volume since
    applywarp operates on a single warp field regardless of input dimensionality,
    but we keep this explicit and consistent with the reference script.

    Args:
        input_path: Path to the input series (3D or 4D)
        rel_warp_path: Path to the shared relative warp (from compute_gnc_warp)
        ref_volume_nii: Path to the reference volume the warp was computed against
        output_path: Path to write the warped series to
        split_prefix: Basename prefix for the per-volume scratch files
    """
    from mrtrix3 import run, image

    n_volumes = int(image.Header(input_path).size()[3]) if len(image.Header(input_path).size()) == 4 else 1

    # fslsplit only understands NIfTI, but input_path may be a .mif (e.g. a
    # user-supplied -rpe_pair), so convert through a NIfTI intermediate first.
    split_input = f'{split_prefix}input.nii.gz'
    run.command(f'mrconvert -force {input_path} {split_input}')
    run.command(f'fslsplit {split_input} {split_prefix} -t')

    warped_files = []
    for i in range(n_volumes):
        idx = f'{i:04d}'
        vol_in = f'{split_prefix}{idx}.nii.gz'
        vol_out = f'{split_prefix}{idx}_gnc.nii.gz'
        run.command(
            f'applywarp --rel --interp=spline --datatype=float '
            f'-i {vol_in} -r {ref_volume_nii} -w {rel_warp_path} -o {vol_out}'
        )
        warped_files.append(vol_out)

    if n_volumes == 1:
        run.command(f'mrconvert -force {warped_files[0]} {output_path}')
    else:
        run.command(f'mrcat -force -axis 3 {" ".join(warped_files)} {output_path}')


def run_gnc_precorrection(dwi_metadata: dict[str, Any]) -> str | None:
    """
    Apply GNC to the pre-eddy working data (and rpe_pair, if present).

    Preserves the pre-GNC (denoised/degibbsed) data as working_predistort.mif
    for later single-interpolation recombination with eddy's displacement
    fields, and overwrites working.mif with the GNC-corrected series so
    topup/eddy operate on geometrically-correct data.

    If app.ARGS.rpe_pair is set, the same shared warp is applied to it and
    the corrected copy's path is returned (with matching .bval/.bvec/.json
    sidecars copied over from the original, if present) so the caller can
    pass it explicitly into run_eddy(). app.ARGS.rpe_pair itself is never
    modified.

    Args:
        dwi_metadata: DESIGNER's input metadata dict (used for 'stride')

    Returns:
        Path to the GNC-corrected rpe_pair image, or None if app.ARGS.rpe_pair
        is not set.
    """
    from mrtrix3 import app, run
    from lib.designer_input_utils import splitext_
    import shutil

    stride = dwi_metadata['stride']

    run.command('mrconvert -force working.mif working_predistort.mif')
    run.command(
        'mrconvert -force -export_grad_fsl working.bvec working.bval working.mif working_gnc_in.nii.gz'
    )
    run.command('mrconvert -force -coord 3 0 -axes 0,1,2 working_gnc_in.nii.gz gnc_ref_vol.nii.gz')

    # 'gnc' is a fixed scratch basename so finalize_single_interpolation() can
    # reuse the same relative warp (gnc_warp_rel.nii.gz) after eddy runs.
    compute_gnc_warp('gnc_ref_vol.nii.gz', app.ARGS.gnc_grad_coeff, app.ARGS.gnc_scanner, 'gnc')

    apply_warp_to_series(
        'working_gnc_in.nii.gz', 'gnc_warp_rel.nii.gz', 'gnc_ref_vol.nii.gz',
        'working_gnc.nii.gz', 'vol_working_')
    run.command(
        f'mrconvert -force -stride {stride} -fslgrad working.bvec working.bval '
        f'working_gnc.nii.gz working.mif'
    )

    if not getattr(app.ARGS, 'rpe_pair', None):
        return None

    rpe_path = app.ARGS.rpe_pair
    rpe_base, _ = splitext_(rpe_path)

    run.command('mrconvert -force gnc_ref_vol.nii.gz rpe_gnc_ref.nii.gz')
    apply_warp_to_series(
        rpe_path, 'gnc_warp_rel.nii.gz', 'rpe_gnc_ref.nii.gz',
        'rpe_pair_gnc.nii.gz', 'vol_rpe_')

    gnc_rpe_path = str(Path.cwd() / 'rpe_pair_gnc.nii.gz')
    gnc_rpe_base, _ = splitext_(gnc_rpe_path)
    # rpe_pair is always a NIfTI image, so its gradient table (if any) lives
    # in sidecar .bval/.bvec files next to it; carry them over by basename.
    for sidecar_ext in ('.json', '.bval', '.bvec'):
        src = rpe_base + sidecar_ext
        if Path(src).exists():
            shutil.copy(src, gnc_rpe_base + sidecar_ext)

    return gnc_rpe_path


def finalize_single_interpolation(
    eddy_out_prefix: str | Path,
    eddy_bvals_file: str | Path,
    n_volumes: int,
    dwi_metadata: dict[str, Any],
    mif_output: str | Path,
) -> str | Path:
    """
    Resample the pre-GNC data once with GNC+eddy warps composed.

    Combines the shared GNC relative warp with each volume's eddy
    displacement field, applies the composed warp in a single spline
    interpolation to working_predistort.mif, and rescales by eddy's own
    per-volume Jacobian only (GNC carries no Jacobian, matching
    gradient_unwarp.py -n).

    Args:
        eddy_out_prefix: Eddy's output basename (scratch_dir/dwi_post_eddy),
            used to locate its per-volume .eddy_displacement_fields.* and
            .eddy_rotated_bvecs
        eddy_bvals_file: Path to the bvals file eddy was run with
        n_volumes: Number of DWI volumes
        dwi_metadata: DESIGNER's input metadata dict (used for 'stride')
        mif_output: Path to write the single-interpolation composite mif to

    Returns:
        mif_output, unchanged, for convenience chaining
    """
    from mrtrix3 import run

    stride = dwi_metadata['stride']

    run.command('fslsplit working_predistort.mif vol_predistort_ -t')

    corrected_files = []
    for i in range(n_volumes):
        idx = f'{i:04d}'
        eddy_field = f'{eddy_out_prefix}.eddy_displacement_fields.{i + 1:03d}.nii.gz'

        full_warp_rel = f'vol_{idx}_full_warp_rel.nii.gz'
        run.command(
            f'convertwarp --ref=gnc_ref_vol.nii.gz --rel '
            f'--warp1=gnc_warp_rel.nii.gz --warp2={eddy_field} '
            f'--relout --out={full_warp_rel}'
        )

        eddy_warp_rel = f'vol_{idx}_eddy_warp_rel.nii.gz'
        eddy_jacobians = f'vol_{idx}_eddy_jacobians.nii.gz'
        run.command(
            f'convertwarp --ref=gnc_ref_vol.nii.gz --rel '
            f'--warp1={eddy_field} --relout --out={eddy_warp_rel} '
            f'--jacobian={eddy_jacobians}'
        )
        eddy_jacobian = f'vol_{idx}_eddy_jacobian.nii.gz'
        run.command(f'mrmath -axis 3 {eddy_jacobians} mean {eddy_jacobian} -force -quiet')

        transformed = f'vol_{idx}_full_transformed.nii.gz'
        run.command(
            f'applywarp --rel --interp=spline --datatype=float '
            f'-i vol_predistort_{idx}.nii.gz -r gnc_ref_vol.nii.gz '
            f'-w {full_warp_rel} -o {transformed}'
        )

        corrected = f'vol_{idx}_full_corrected.nii.gz'
        run.command(
            f'mrcalc {transformed} {eddy_jacobian} -multiply {corrected} -force -quiet'
        )
        corrected_files.append(corrected)

    run.command(f'mrcat -force -axis 3 {" ".join(corrected_files)} composite.nii.gz')
    run.command(
        f'mrconvert -force -stride {stride} '
        f'-fslgrad {eddy_out_prefix}.eddy_rotated_bvecs {eddy_bvals_file} '
        f'composite.nii.gz {mif_output}'
    )
    return mif_output

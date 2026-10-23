"""Gradient nonlinearity correction (GNC): spatial unwarping, with a choice of
backend ('gradunwarp' or 'tortoise').

Computes a single geometric warp (per scan) from the scanner's gradient
coefficient file and applies it either as a standalone pre-eddy correction
(for topup/eddy to operate on geometrically-correct data) or composed with
eddy's own per-volume displacement fields so the combined GNC+eddy correction
is resampled in a single interpolation pass.

Both backends also produce a per-voxel gradient-deviation ("L") tensor field
(see derive_gnc_graddev), consumed by designer2/tmi.py's -gnc option to
correct the b-matrix itself for gradient nonlinearity, not just image
geometry. Each backend only has to satisfy two contracts so the rest of this
module (and eddy composition in finalize_single_interpolation) never needs to
know which one ran: compute_gnc_warp() returns an FSL relative-warp file, and
derive_gnc_graddev() returns a NIfTI with last dimension 9 (a flattened 3x3
matrix per voxel, L ~= Identity at undistorted voxels).
"""

import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from lib.utils import force_remove


def _read_nonlinear_gradient_correction(sidecar: str | Path) -> bool | None:
    """Return the BIDS NonlinearGradientCorrection value, or None if unavailable."""
    if not Path(sidecar).is_file():
        return None
    with open(sidecar) as f:
        return json.load(f).get('NonlinearGradientCorrection')


def resolve_gnc_request(sidecars: Sequence[str | Path]) -> bool:
    """
    Decide whether a requested GNC should run, based on the BIDS
    NonlinearGradientCorrection field of every input that GNC would warp.

    All sidecars must agree: a sidecar that is missing (or lacks the field)
    is assumed to match the others, and if none report it the user's request
    is followed. GNC is skipped if the scanner already applied it.

    Args:
        sidecars: BIDS .json paths of the main DWI series (and the RPE image,
            when it is used by eddy)

    Returns:
        True if GNC should run, False if it should be skipped
    """
    import warnings
    from colorama import Fore
    from mrtrix3 import MRtrixError

    states = {str(s): _read_nonlinear_gradient_correction(s) for s in sidecars}
    known = {state for state in states.values() if state is not None}

    if len(known) > 1:
        details = ', '.join(f'{s}: {state}' for s, state in states.items() if state is not None)
        raise MRtrixError(
            'Inputs have inconsistent gradient nonlinearity correction states '
            f'(NonlinearGradientCorrection - {details}); cannot apply -gnc_grad_coeff.')

    if not known:
        warnings.warn(Fore.RED + '[WARNING] NonlinearGradientCorrection not found in BIDS sidecar(s); '
                      'cannot verify whether the scanner already applied gradient nonlinearity '
                      'correction. Proceeding with -gnc_grad_coeff.')
        return True

    already_corrected = known.pop()
    unknown = [s for s, state in states.items() if state is None]
    if unknown:
        warnings.warn(Fore.RED + f'[WARNING] NonlinearGradientCorrection not found for {", ".join(unknown)}; '
                      f'assuming {already_corrected} to match the other input(s).')

    if already_corrected:
        warnings.warn(Fore.RED + '[WARNING] BIDS sidecar(s) report NonlinearGradientCorrection=true '
                      '(already corrected on the scanner); ignoring -gnc_grad_coeff.')
    return not already_corrected


def compute_gnc_warp(
    ref_volume_nii: str | Path,
    coeff_file: str | Path,
    scanner: str,
    out_prefix: str,
    backend: str = 'gradunwarp',
) -> str:
    """
    Compute the relative FSL warp for gradient nonlinearity correction.

    Dispatches to the requested backend; both produce a warp in FSL's
    relative-warp convention (for use with applywarp), so nothing downstream
    (apply_warp_to_series, finalize_single_interpolation) needs to know which
    backend ran. Must be followed by a derive_gnc_graddev() call with the
    same ref_volume_nii/out_prefix/backend before any other GNC intermediates
    are produced - the tortoise backend deliberately leaves a byproduct file
    behind for that call to consume (the same kind of fixed-scratch-name
    coupling run_gnc_precorrection's 'gnc' basename already relies on).

    Args:
        ref_volume_nii: Path to the single reference volume to unwarp (the
            warp is static in scanner space, so one volume is enough)
        coeff_file: Path to the scanner gradient coefficient file
        scanner: Scanner vendor ('siemens' or 'ge')
        out_prefix: Basename prefix for the intermediate/output warp files
        backend: 'gradunwarp' (default) or 'tortoise'

    Returns:
        Path to the relative warp file (FSL convention, for use with applywarp)
    """
    from mrtrix3 import MRtrixError

    if backend == 'gradunwarp':
        return _compute_gnc_warp_gradunwarp(ref_volume_nii, coeff_file, scanner, out_prefix)
    if backend == 'tortoise':
        return _compute_gnc_warp_tortoise(ref_volume_nii, coeff_file, scanner, out_prefix)
    raise MRtrixError(f"Unknown GNC backend '{backend}'")


def _compute_gnc_warp_gradunwarp(
    ref_volume_nii: str | Path,
    coeff_file: str | Path,
    scanner: str,
    out_prefix: str,
) -> str:
    """Runs gradient_unwarp.py on the reference volume and converts its
    absolute warp field to FSL's relative-warp convention."""
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


def resolve_tortoise_is_ge(scanner: str) -> bool:
    """Map -gnc_scanner to TORTOISE's GE/non-GE distinction (--isGE flag).

    A plain equality check is safe here because -gnc_scanner is constrained
    to {'siemens', 'ge'} by argparse (designer2/designer.py), matching
    gradunwarp's own hard vendor constraint - no third value can reach here.
    """
    return scanner == 'ge'


def _compute_gnc_warp_tortoise(
    ref_volume_nii: str | Path,
    coeff_file: str | Path,
    scanner: str,
    out_prefix: str,
) -> str:
    """Runs TORTOISE's CreateGradientNonlinearityBMatrix on the reference
    volume - a single call that is a strict superset of
    CreateNonlinearityDisplacementMap, writing both an ITK displacement field
    and the per-voxel L tensor - then converts the displacement field to
    FSL's relative-warp convention (see _itk_field_to_fsl_rel_warp). The L tensor is deliberately left on disk
    for derive_gnc_graddev() to pick up right after this returns."""
    from mrtrix3 import run
    from lib.designer_input_utils import splitext_

    isge_flag = ' --isGE' if resolve_tortoise_is_ge(scanner) else ''
    run.command(f'CreateGradientNonlinearityBMatrix -f {ref_volume_nii} -g {coeff_file}{isge_flag}')

    # Output names derive from ref_volume_nii's own basename, not out_prefix.
    ref_stem, _ = splitext_(str(ref_volume_nii))
    itk_field = f'{ref_stem}_gradwarp_field.nii'

    rel_warp = f'{out_prefix}_warp_rel.nii.gz'
    _itk_field_to_fsl_rel_warp(itk_field, ref_volume_nii, rel_warp)
    force_remove(itk_field)
    return rel_warp


def _lps_disp_to_fsl_rel(disp_lps, affine):
    """
    Convert per-voxel displacements (mm, components in the ITK LPS world
    frame, as TORTOISE writes them) into an FSL relative warp on the same grid
    (mm along the image's storage axes, x negated when det(affine) > 0, per
    FSL's internal radiological convention). Exact for oblique affines.

    disp_lps: (..., 3) array. affine: 4x4 voxel-to-RAS-world matrix.
    """
    import numpy as np

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
    import numpy as np
    import nibabel as nib
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


def derive_gnc_graddev(
    ref_volume_nii: str | Path,
    coeff_file: str | Path,
    scanner: str,
    out_prefix: str,
    backend: str,
) -> str:
    """
    Produce {out_prefix}_graddev.nii.gz: the per-voxel gradient-deviation
    tensor L (L ~= Identity at undistorted voxels) consumed by
    designer2/tmi.py's -gnc option. Must be called immediately after
    compute_gnc_warp() with the same arguments.

    Args:
        ref_volume_nii: Same reference volume passed to compute_gnc_warp()
        coeff_file: Same coefficient file passed to compute_gnc_warp()
        scanner: Same scanner vendor passed to compute_gnc_warp()
        out_prefix: Same basename prefix passed to compute_gnc_warp()
        backend: 'gradunwarp' or 'tortoise' - same as compute_gnc_warp()

    Returns:
        Path to the sanitized graddev tensor field
    """
    from mrtrix3 import MRtrixError
    from lib.designer_input_utils import splitext_

    if backend == 'tortoise':
        # Byproduct of _compute_gnc_warp_tortoise()'s single
        # CreateGradientNonlinearityBMatrix call, left on disk for this.
        ref_stem, _ = splitext_(str(ref_volume_nii))
        raw = f'{ref_stem}_graddev_c.nii'
    elif backend == 'gradunwarp':
        raw = _derive_graddev_gradunwarp(ref_volume_nii, coeff_file, scanner, out_prefix)
    else:
        raise MRtrixError(f"Unknown GNC backend '{backend}'")

    out_path = f'{out_prefix}_graddev.nii.gz'
    _sanitize_graddev_nan(raw, out_path)
    force_remove(raw)
    return out_path


def _derive_graddev_gradunwarp(
    ref_volume_nii: str | Path,
    coeff_file: str | Path,
    scanner: str,
    out_prefix: str,
) -> str:
    """gradunwarp has no native L tensor; derive one with FSL's
    calc_grad_perc_dev from gradunwarp's absolute warp, computed on the
    reference grid padded by one voxel and cropped back afterwards.
    calc_grad_perc_dev zero-pads outside the grid, so differentiating
    absolute coordinates on the unpadded grid corrupts the whole outermost
    voxel layer (values ~50-150 instead of ~1)."""
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


def _sanitize_graddev_matrix(mat):
    """Replace any 3x3 gradient-deviation block containing a NaN with the
    identity matrix as a whole block (not entrywise). Works around a
    confirmed TORTOISE bug: CreateGradientNonlinearityBMatrix emits NaN for
    the entire 3x3 block exactly at the isocenter voxel, where the true
    value is Identity. Accepts (...,9) or (...,3,3); returns the same shape.
    """
    import numpy as np

    flat = mat.reshape(*mat.shape[:-1], 3, 3) if mat.shape[-1] == 9 else mat
    bad = np.isnan(flat).any(axis=(-2, -1))
    flat = np.where(bad[..., None, None], np.eye(3, dtype=flat.dtype), flat)
    return flat.reshape(mat.shape)


def _sanitize_graddev_nan(raw_path: str | Path, out_path: str | Path) -> None:
    import numpy as np
    import nibabel as nib

    img = nib.load(str(raw_path))
    data = _sanitize_graddev_matrix(np.asarray(img.dataobj)).astype(np.float32)
    nib.save(nib.Nifti1Image(data, img.affine), str(out_path))


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
    force_remove(split_input)

    warped_files = []
    for i in range(n_volumes):
        idx = f'{i:04d}'
        vol_in = f'{split_prefix}{idx}.nii.gz'
        vol_out = f'{split_prefix}{idx}_gnc.nii.gz'
        run.command(
            f'applywarp --rel --interp=spline --datatype=float '
            f'-i {vol_in} -r {ref_volume_nii} -w {rel_warp_path} -o {vol_out}'
        )
        force_remove(vol_in)
        warped_files.append(vol_out)

    if n_volumes == 1:
        run.command(f'mrconvert -force {warped_files[0]} {output_path}')
    else:
        run.command(f'mrcat -force -axis 3 {" ".join(warped_files)} {output_path}')
    force_remove(*warped_files)


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
    # reuse the same relative warp (gnc_warp_rel.nii.gz) after eddy runs, and
    # so the final output stage can export the L tensor (gnc_graddev.nii.gz).
    backend = app.ARGS.gnc_backend
    coeff, scanner = app.ARGS.gnc_grad_coeff, app.ARGS.gnc_scanner
    compute_gnc_warp('gnc_ref_vol.nii.gz', coeff, scanner, 'gnc', backend=backend)
    derive_gnc_graddev('gnc_ref_vol.nii.gz', coeff, scanner, 'gnc', backend)

    apply_warp_to_series(
        'working_gnc_in.nii.gz', 'gnc_warp_rel.nii.gz', 'gnc_ref_vol.nii.gz',
        'working_gnc.nii.gz', 'vol_working_')
    force_remove('working_gnc_in.nii.gz')
    run.command(
        f'mrconvert -force -stride {stride} -fslgrad working.bvec working.bval '
        f'working_gnc.nii.gz working.mif'
    )
    force_remove('working_gnc.nii.gz')

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

    # fslsplit only understands NIfTI, so convert through a NIfTI intermediate first.
    run.command('mrconvert -force working_predistort.mif working_predistort.nii.gz')
    run.command('fslsplit working_predistort.nii.gz vol_predistort_ -t')
    force_remove('working_predistort.nii.gz')

    corrected_files = []
    for i in range(n_volumes):
        idx = f'{i:04d}'
        eddy_field = f'{eddy_out_prefix}.eddy_displacement_fields.{i + 1:03d}.nii.gz'
        vol_predistort = f'vol_predistort_{idx}.nii.gz'

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
            f'-i {vol_predistort} -r gnc_ref_vol.nii.gz '
            f'-w {full_warp_rel} -o {transformed}'
        )
        force_remove(vol_predistort)

        corrected = f'vol_{idx}_full_corrected.nii.gz'
        run.command(
            f'mrcalc {transformed} {eddy_jacobian} -multiply {corrected} -force -quiet'
        )
        force_remove(full_warp_rel, eddy_warp_rel, eddy_jacobians, eddy_jacobian, transformed)
        corrected_files.append(corrected)

    run.command(f'mrcat -force -axis 3 {" ".join(corrected_files)} composite.nii.gz')
    run.command(
        f'mrconvert -force -stride {stride} '
        f'-fslgrad {eddy_out_prefix}.eddy_rotated_bvecs {eddy_bvals_file} '
        f'composite.nii.gz {mif_output}'
    )
    force_remove(*corrected_files, 'composite.nii.gz')
    return mif_output

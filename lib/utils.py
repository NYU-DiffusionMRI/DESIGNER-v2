from pathlib import Path
from typing import Optional
import numpy as np


def write_manual_pe_scheme(
    outfile: str | Path,
    pe_dir: str,
    n_volumes: int,
    readout_time: float = 0.1,
):
    """
    Reproduces dwifslpreproc behavior when:
    - pe_dir is provided
    - readout_time is NOT provided (defaults to 0.1)
    """
    from mrtrix3 import phaseencoding

    pe_vec = list(phaseencoding.direction(pe_dir))    # [dx, dy, dz]

    with open(outfile, "w") as f:
        for _ in range(n_volumes):
            f.write(f"{pe_vec[0]} {pe_vec[1]} {pe_vec[2]} {readout_time}\n")


def run_fsl_eddy(
    mif_input: str | Path,
    mif_output: str | Path,
    brain_mask: str | Path,
    scratch_dir: str | Path,
    *,
    eddy_opts: str = '',
    topup_prefix: Optional[str] = None,
    pe_dir: Optional[str] = None,
    grad_file: Optional[str] = None,
):
    from mrtrix3 import run, image

    # 1. convert the nifti brainmask in a format and orientation that eddy can use
    eddy_brain_mask = scratch_dir / 'eddy_mask.nii'
    run.command(f'mrconvert {brain_mask} {eddy_brain_mask} -datatype float32 -stride -1,+2,+3')

    # 2. convert the mif input to a nifti input
    nifti_input = scratch_dir / 'eddy_in.nii'
    cmd = f'mrconvert {mif_input} {nifti_input} -strides -1,+2,+3,+4 -export_grad_fsl {scratch_dir}/bvecs {scratch_dir}/bvals -export_pe_eddy {scratch_dir}/eddy_config.txt {scratch_dir}/eddy_indices.txt'

    if pe_dir is not None:
        manual_pe_scheme_path = scratch_dir / 'dwi_manual_pe_scheme.txt'
        n_volumes = int(image.Header(mif_input).size()[3])
        write_manual_pe_scheme(manual_pe_scheme_path, pe_dir, n_volumes)
        cmd += f' -import_pe_table {manual_pe_scheme_path}'

    if grad_file is not None:
        cmd += f' -grad {grad_file}'

    run.command(cmd)

    # 3. call FSL eddy
    eddy_out_wo_ext = scratch_dir / 'dwi_post_eddy'
    eddy_acqp_file = scratch_dir / 'eddy_config.txt'
    eddy_indices_file = scratch_dir / 'eddy_indices.txt'
    eddy_bvecs_file = scratch_dir / 'bvecs'
    eddy_bvals_file = scratch_dir / 'bvals'

    cmd = (
        f"eddy --imain={nifti_input} --mask={eddy_brain_mask} --acqp={eddy_acqp_file} "
        f"--index={eddy_indices_file} --bvecs={eddy_bvecs_file} --bvals={eddy_bvals_file} "
        f"--out={eddy_out_wo_ext} {eddy_opts} --verbose"
    )
    if topup_prefix is not None:
        cmd += f" --topup={topup_prefix}"
    
    run.command(cmd)

    # 4. convert the eddy output nifti to mif
    run.command(f'mrconvert {eddy_out_wo_ext}.nii.gz {mif_output} -strides -1,2,3,4 -fslgrad {eddy_out_wo_ext}.eddy_rotated_bvecs {eddy_bvals_file}')


def run_topup_and_prepare_for_eddy(
    b0_pair_file: str,
    pe_dir: str,
    output_prefix: str,
    fsl_suffix: str,
    scratch_dir: str | Path,
    *,
    acqp_file: Optional[str] = None,
) -> tuple[str, str]:
    """
    Run FSL topup and prepare brain mask for eddy.
    
    Args:
        b0_pair_file: Path to concatenated b0 pair for topup
        pe_dir: Phase encoding direction (e.g., 'j', 'j-', 'i', etc.)
        output_prefix: Prefix for output files (e.g., 'topup_results_1')
        fsl_suffix: FSL suffix from fsl.suffix()
        scratch_dir: Directory to store output files (mimics dwifslpreproc scratch dir)
        acqp_file: Optional path to existing acquisition parameters file
        
    Returns:
        tuple of (brain_mask_path, topup_prefix) for use with run_fsl_eddy
    """
    from mrtrix3 import run, image
    
    # Set up paths based on scratch directory
    scratch_dir = Path(scratch_dir)
    topup_prefix = scratch_dir / output_prefix
    topup_acqp = scratch_dir / f'{output_prefix}_acqp.txt'
    
    # 1. Create topup acquisition parameters if not provided
    if acqp_file is None:
        acqp_file = str(topup_acqp)
        acqp = np.zeros((2, 3))
        if 'i' in pe_dir: acqp[:, 0] = 1
        if 'j' in pe_dir: acqp[:, 1] = 1
        if 'k' in pe_dir: acqp[:, 2] = 1
        if '-' in pe_dir:
            acqp[0, :] = -acqp[0, :]
        else:
            acqp[1, :] = -acqp[1, :]
        
        acqp[acqp == -0] = 0
        acqp = np.hstack((acqp, np.array([0.1, 0.1])[..., None]))
        np.savetxt(acqp_file, acqp, fmt="%1.2f")
    
    # 2. Check for odd dimensions (affects topup subsampling)
    odd_dims = [int(s) for s in image.Header(b0_pair_file).size()[:3] if s % 2]
    subsamp_flag = '--subsamp=1' if np.any(np.array(odd_dims)) else ''
    
    # 3. Run topup
    run.command(
        f'topup --imain="{b0_pair_file}" --datain="{acqp_file}" '
        f'--config=b02b0.cnf {subsamp_flag} --scale=1 '
        f'--out="{topup_prefix}" --iout="{topup_prefix}{fsl_suffix}"'
    )
    
    # 4. Create brain mask from topup corrected image
    topup_mean = f'{topup_prefix}_mean.nii'
    topup_brain_prefix = f'{topup_prefix}_brain'
    topup_brain_mask = f'{topup_brain_prefix}_mask.nii.gz'
    
    run.command(f'mrmath "{topup_prefix}{fsl_suffix}" mean "{topup_mean}" -axis 3')
    run.command(f'bet "{topup_mean}" "{topup_brain_prefix}" -f 0.2 -m')
    
    return topup_brain_mask, str(topup_prefix)

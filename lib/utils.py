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
    Write PE scheme for single-direction acquisition.
    
    Args:
        outfile: Output file path
        pe_dir: Phase encoding direction (e.g., 'AP', 'PA', 'LR')
        n_volumes: Number of volumes
        readout_time: Total readout time (default: 0.1)
    """
    from mrtrix3 import phaseencoding

    pe_vec = list(phaseencoding.direction(pe_dir))    # [dx, dy, dz]

    with open(outfile, "w") as f:
        for _ in range(n_volumes):
            f.write(f"{pe_vec[0]} {pe_vec[1]} {pe_vec[2]} {readout_time}\n")


def write_rpe_all_pe_scheme(
    outfile: str | Path,
    pe_dir: str,
    n_volumes_forward: int,
    n_volumes_reverse: int,
    readout_time: float = 0.1,
):
    """
    Write PE scheme for -rpe_all acquisition (forward + reverse PE).
    
    Args:
        outfile: Output file path
        pe_dir: Phase encoding direction for forward acquisition (e.g., 'AP')
        n_volumes_forward: Number of forward PE volumes
        n_volumes_reverse: Number of reverse PE volumes
        readout_time: Total readout time (default: 0.1)
    """
    from mrtrix3 import phaseencoding

    pe_vec_fwd = list(phaseencoding.direction(pe_dir))
    pe_vec_rev = [-x for x in pe_vec_fwd]  # Reverse direction

    with open(outfile, "w") as f:
        # Forward PE volumes
        for _ in range(n_volumes_forward):
            f.write(f"{pe_vec_fwd[0]} {pe_vec_fwd[1]} {pe_vec_fwd[2]} {readout_time}\n")
        # Reverse PE volumes
        for _ in range(n_volumes_reverse):
            f.write(f"{pe_vec_rev[0]} {pe_vec_rev[1]} {pe_vec_rev[2]} {readout_time}\n")


def compute_jacobian_weight_for_rpe_all(
    field_map_path: str | Path,
    pe_config_row: np.ndarray,
    output_suffix: str,
    scratch_dir: str | Path = '.',
):
    """
    Compute Jacobian-based weight for -rpe_all volume recombination.
    
    Creates weight image based on distortion field Jacobian.
    Weight = J² where J = 1 + ∂(field)/∂(PE_direction)
    Higher weight indicates less distortion → higher contribution in recombination.
    
    Args:
        field_map_path: Path to topup field map image
        pe_config_row: Single row from eddy_config.txt [pe_x, pe_y, pe_z, readout_time]
        output_suffix: Suffix for output weight file (e.g., 'fwd', 'rev')
        scratch_dir: Directory for output file (default: current directory)
    """
    from mrtrix3 import run, app
    
    pe_axis = int(np.where(pe_config_row[:3] != 0)[0][0])   # 0=x,1=y,2=z
    pe_sign = float(pe_config_row[pe_axis])                 # +1 or -1
    readout = float(pe_config_row[3])

    scratch_dir = Path(scratch_dir)
    jac = scratch_dir / f"jac_{output_suffix}.mif"
    wgt = scratch_dir / f"weight_{output_suffix}.mif"

    sign_mult = " -1.0 -mult" if pe_sign < 0 else ""
    
    app.console(f"Computing Jacobian-based weight for {output_suffix} PE direction")

    # J = max(0, 1 + d( field*readout*sign ) / dPE)
    run.command(
        f"mrcalc {field_map_path} {readout} -mult{sign_mult} - | "
        f"mrfilter - gradient - | "
        f"mrconvert - -coord 3 {pe_axis} -axes 0,1,2 - | "
        f"mrcalc 1.0 - -add 0.0 -max {jac}",
    )

    # weight = J^2
    run.command(f"mrcalc {jac} {jac} -mult {wgt}")
    run.function(Path.unlink, jac, show=False)


def run_fsl_eddy(
    mif_input: str | Path,
    mif_output: str | Path,
    brain_mask: str | Path,
    scratch_dir: Path,
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
    include_field_map: bool = False,
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
    fout_flag = f'--fout={topup_prefix}_field_map.nii.gz' if include_field_map else ''

    run.command(
        f'topup --imain="{b0_pair_file}" --datain="{acqp_file}" '
        f'--config=b02b0.cnf {subsamp_flag} --scale=1 '
        f'--out="{topup_prefix}" --iout="{topup_prefix}{fsl_suffix}" {fout_flag}'
    )
    
    # 4. Create brain mask from topup corrected image
    topup_mean = f'{topup_prefix}_mean.nii'
    topup_brain_prefix = f'{topup_prefix}_brain'
    topup_brain_mask = f'{topup_brain_prefix}_mask.nii.gz'
    
    run.command(f'mrmath "{topup_prefix}{fsl_suffix}" mean "{topup_mean}" -axis 3')
    run.command(f'bet "{topup_mean}" "{topup_brain_prefix}" -f 0.2 -m')
    
    return topup_brain_mask, str(topup_prefix)

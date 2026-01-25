from pathlib import Path
from typing import Optional


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


def call_eddy(
    mif_input: str | Path,
    mif_output: str | Path,
    eddy_proc_dir: str | Path,
    brain_mask: str | Path,
    *,
    eddy_opts: str = '',
    pe_dir: Optional[str] = None,
    grad_file: Optional[str] = None,
):
    from mrtrix3 import run, image

    # 1. convert the nifti brainmask in a format and orientation that eddy can use
    eddy_brain_mask = eddy_proc_dir / 'eddy_mask.nii'
    run.command(f'mrconvert {brain_mask} {eddy_brain_mask} -datatype float32 -stride -1,+2,+3')

    # 2. convert the mif input to a nifti input
    nifti_input = eddy_proc_dir / 'eddy_in.nii'
    cmd = f'mrconvert {mif_input} {nifti_input} -strides -1,+2,+3,+4 -export_grad_fsl {eddy_proc_dir}/bvecs {eddy_proc_dir}/bvals -export_pe_eddy {eddy_proc_dir}/eddy_config.txt {eddy_proc_dir}/eddy_indices.txt'

    if pe_dir is not None:
        manual_pe_scheme_path = eddy_proc_dir / 'dwi_manual_pe_scheme.txt'
        n_volumes = int(image.Header(mif_input).size()[3])
        write_manual_pe_scheme(manual_pe_scheme_path, pe_dir, n_volumes)
        cmd += f' -import_pe_table {manual_pe_scheme_path}'

    if grad_file is not None:
        cmd += f' -grad {grad_file}'

    run.command(cmd)

    # 3. call FSL eddy
    eddy_out_wo_ext = eddy_proc_dir / 'dwi_post_eddy'
    eddy_acqp_file = eddy_proc_dir / 'eddy_config.txt'
    eddy_indices_file = eddy_proc_dir / 'eddy_indices.txt'
    eddy_bvecs_file = eddy_proc_dir / 'bvecs'
    eddy_bvals_file = eddy_proc_dir / 'bvals'

    run.command(
        f"eddy --imain={nifti_input} --mask={eddy_brain_mask} --acqp={eddy_acqp_file} "
        f"--index={eddy_indices_file} --bvecs={eddy_bvecs_file} --bvals={eddy_bvals_file} "
        f"--topup=topup_results --out={eddy_out_wo_ext} --verbose {eddy_opts}"
    )

    # 4. convert the eddy output nifti to mif
    run.command(f'mrconvert {eddy_out_wo_ext}.nii.gz {mif_output} -strides -1,2,3,4 -fslgrad {eddy_out_wo_ext}.eddy_rotated_bvecs {eddy_bvals_file}')

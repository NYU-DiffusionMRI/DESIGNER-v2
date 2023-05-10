#!/usr/bin/env python3

import os, sys
import numpy as np, gzip, shutil
__version__ = "0.0.1"

from lib.designer_input_utils import *
from lib.designer_func_wrappers import *

def usage(cmdline): #pylint: disable=unused-variable
    from mrtrix3 import app #pylint: disable=no-name-in-module, import-outside-toplevel
    cmdline.set_author('Benjamin Ades-Aron (benjamin.ades-aron@nyulangone.org)')
    cmdline.set_synopsis("""1. pre-check: concatenate all dwi series and make sure the all diffusion AP images and PA image have the same matrix dimensions/orientations 
 						
                        2. denoise the complete dataset\n
						
                        3. gibbs ringing correction on complete dataset\n
 						 						
                        4. topup + eddy\n
 						
                        5. b0 normalization\n
 						                        
                        6. Rician bias correction\n
                         						
                        7. Fitting and or outlier correction\n
 						
                        """)

    cmdline.add_citation('Veraart, J.; Novikov, D.S.; Christiaens, D.; Ades-aron, B.; Sijbers, J. & Fieremans, E. Denoising of diffusion MRI using random matrix theory. NeuroImage, 2016, 142, 394-406, doi: 10.1016/j.neuroimage.2016.08.016',is_external=True)
    cmdline.add_citation('Veraart, J.; Fieremans, E. & Novikov, D.S. Diffusion MRI noise mapping using random matrix theory. Magn. Res. Med., 2016, 76(5), 1582-1593, doi:10.1002/mrm.26059',is_external=True)
    cmdline.add_citation('Kellner, E., et al., Gibbs-Ringing Artifact Removal Based on Local Subvoxel-Shifts. Magnetic Resonance in Medicine, 2016. 76(5): p. 1574-1581.',is_external=True)
    cmdline.add_citation('Koay, C.G. and P.J. Basser, Analytically exact correction scheme for signal extraction from noisy magnitude MR signals. Journal of Magnetic Resonance, 2006. 179(2): p. 317-322.',is_external=True)
    cmdline.add_citation('Andersson, J. L. & Sotiropoulos, S. N. An integrated approach to correction for off-resonance effects and subject movement in diffusion MR imaging. NeuroImage, 2015, 125, 1063-1078', is_external=True)
    cmdline.add_citation('Smith, S. M.; Jenkinson, M.; Woolrich, M. W.; Beckmann, C. F.; Behrens, T. E.; Johansen-Berg, H.; Bannister, P. R.; De Luca, M.; Drobnjak, I.; Flitney, D. E.; Niazy, R. K.; Saunders, J.; Vickers, J.; Zhang, Y.; De Stefano, N.; Brady, J. M. & Matthews, P. M. Advances in functional and structural MR image analysis and implementation as FSL. NeuroImage, 2004, 23, S208-S219', is_external=True)
    cmdline.add_citation('Skare, S. & Bammer, R. Jacobian weighting of distortion corrected EPI data. Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 5063', is_external=True)
    cmdline.add_citation('Andersson, J. L.; Skare, S. & Ashburner, J. How to correct susceptibility distortions in spin-echo echo-planar images: application to diffusion tensor imaging. NeuroImage, 2003, 20, 870-888', is_external=True)
    cmdline.add_citation('Zhang, Y.; Brady, M. & Smith, S. Segmentation of brain MR images through a hidden Markov random field model and the expectation-maximization algorithm. IEEE Transactions on Medical Imaging, 2001, 20, 45-57', is_external=True)
    cmdline.add_citation('Smith, S. M.; Jenkinson, M.; Woolrich, M. W.; Beckmann, C. F.; Behrens, T. E.; Johansen-Berg, H.; Bannister P. R.; De Luca, M.; Drobnjak, I.; Flitney, D. E.; Niazy, R. K.; Saunders, J.; Vickers, J.; Zhang, Y.; DeStefano, N.; Brady, J. M. & Matthews, P. M. Advances in functional and structural MR image analysis and implementation as FSL. NeuroImage, 2004,23, S208-S219', is_external=True)
    cmdline.add_citation('Collier, Q., et al., Iterative reweighted linear least squares for accurate, fast, and robust estimation of diffusion magnetic resonance parameters. Magn Reson Med, 2015. 73(6): p. 2174-84.', is_external=True)
    
    cmdline.add_argument('input',  help='The input DWI series. For multiple input series, separate file names with commas (i.e. dwi1.nii,dwi2.nii,...)')
    cmdline.add_argument('output', help='The output basename')
    
    options = cmdline.add_argument_group('Other options for the DESIGNER script')
    options.add_argument('-degibbs', action='store_true', help='Perform (RPG) Gibbs artifact correction. Must include PF factor with -pf (e.g. 6/8, 7/8) and PF dimension -dim (1, 2 or 3 for i, j  or k respectively).')
    options.add_argument('-rician', action='store_true', help='Perform Rician bias correction')
    options.add_argument('-eddy', action='store_true', help='run fsl eddy (note that if you choose this command you must also choose a phase encoding option')
    options.add_argument('-b1correct', action='store_true', help='Include a bias correction step in dwi preprocessing', default=False)
    options.add_argument('-normalize', action='store_true', help='normalize the dwi volume to median b0 CSF intensity of 1000 (useful for multiple dwi acquisitions)', default=False)
    options.add_argument('-mask', action='store_true',help='compute a brain mask prior to tensor fitting to stip skull and improve efficientcy')
    
    options.add_argument('-datatype', metavar=('<spec>'), help='If using the "-processing_only" option, you can specify the output datatype. Valid options are float32, float32le, float32be, float64, float64le, float64be, int64, uint64, int64le, uint64le, int64be, uint64be, int32, uint32, int32le, uint32le, int32be, uint32be, int16, uint16, int16le, uint16le, int16be, uint16be, cfloat32, cfloat32le, cfloat32be, cfloat64, cfloat64le, cfloat64be, int8, uint8, bit')
    options.add_argument('-fslbvec',metavar=('<bvecs>'),help='specify bvec path if path is different from the path to the dwi or the file has an unusual extention')
    options.add_argument('-fslbval',metavar=('<bvals>'),help='specify bval path if path is different from the path to the dwi or the file has an unusual extention')
    options.add_argument('-bids',metavar=('<bids1,bids2,...>'),help='specify bids.json path if path is different from the path to the dwi or the file has an unusual extention')
    options.add_argument('-n_cores',metavar=('<ncores>'),help='specify the number of cores to use in parallel tasts, by default designer will use available cores - 2', default=-3)
    options.add_argument('-echo_time',metavar=('<TE1,TE2,...>'),help='specify the echo time used in the acquisition (comma separated list the same length as number of inputs)')
    options.add_argument('-bshape',metavar=('<beta1,beta2,...>'),help='specify the b-shape used in the acquisition (comma separated list the same length as number of inputs)')
    options.add_argument('-pre_align',action='store_true',help='rigidly align each input series to correct for large motion between series')
    options.add_argument('-ants_motion_correction',action='store_true',help='perform rigid motion correction using ANTs (useful for cases where eddy breaks down)')

    options.add_argument('-pe_dir', metavar=('<phase encoding direction>'), help='Specify the phase encoding direction of the input series (required if using the eddy option). Can be a signed axis number (e.g. -0, 1, +2), an axis designator (e.g. RL, PA, IS), or NIfTI axis codes (e.g. i-, j, k)')
    options.add_argument('-pf', metavar=('<PF factor>'), help='Specify the partial fourier factor (e.g. 7/8, 6/8)')

    mp_options = cmdline.add_argument_group('Options for specifying denoising parameters')
    mp_options.add_argument('-patch2self', action='store_true', help='Perform patch2self')
    mp_options.add_argument('-denoise', action='store_true', help='Perform MPPCA')
    mp_options.add_argument('-shrinkage',metavar=('<shrink>'),help='specify shrinkage type for MPPCA. Options are "threshold" or "frob"', default='frob')
    mp_options.add_argument('-algorithm',metavar=('<alg>'),help='specify MP algorithm. Options are "veraart","cordero-grande","jespersen"',default='cordero-grande')
    mp_options.add_argument('-extent', metavar=('<size>'), help='MPPCA Denoising extent. Default is 5,5,5')
    mp_options.add_argument('-phase', metavar=('<image>'), help='Diffusion phases', default=None)
    mp_options.add_argument('-adaptive_patch', action='store_true', help='Run MPPCA with adaptive patching')
    
    rpe_options = cmdline.add_argument_group('Options for eddy and to specify the acquisition phase-encoding design')
    rpe_options.add_argument('-eddy_groups',metavar=('<index1,index2,...'),help='specify how input series should be grouped when running eddy, for use with variable TE or b-shape data. (Comma separared list of intergers beginning with 1, i.e. 1,1,1,2,2,2 for 6 series where the first 3 series and second 3 series have different echo times).', default=None)
    rpe_options.add_argument('-rpe_none', action='store_true', help='Specify that no reversed phase-encoding image data is being provided; eddy will perform eddy current and motion correction only')
    rpe_options.add_argument('-rpe_pair', metavar=('<reverse PE b=0 image>'), help='Specify the reverse phase encoding image')
    rpe_options.add_argument('-rpe_all', metavar=('<reverse PE dwi volume>'), help='Specify that ALL DWIs have been acquired with opposing phase-encoding; this information will be used to perform a recombination of image volumes (each pair of volumes with the same b-vector but different phase encoding directions will be combined together into a single volume). The argument to this option is the set of volumes with reverse phase encoding but the same b-vectors as the input image')
    rpe_options.add_argument('-rpe_header', action='store_true', help='Specify that the phase-encoding information can be found in the image header(s), and that this is the information that the script should use')
    rpe_options.add_argument('-rpe_te', metavar=('<echo time (s)>'), help='Specify the echo time of the reverse phase encoded image, if it is not accompanied by a bids .json sidecar.')

def execute(): #pylint: disable=unused-variable
    from mrtrix3 import app, fsl, run, path #pylint: disable=no-name-in-module, import-outside-toplevel
    import pandas as pd

    # create a temparary directory to store processing files
    app.make_scratch_dir()
    
    # grab the fsl suffix stored in $FSLOUTPUTTYPLE
    fsl_suffix = fsl.suffix()

    # create dict containing metadata from either user input args or bids .json sidecars
    dwi_metadata = get_input_info(
        app.ARGS.input, app.ARGS.fslbval, app.ARGS.fslbvec, app.ARGS.bids)
    
    # ensure inputs are reasonable for subsequent dwi processing
    assert_inputs(dwi_metadata, app.ARGS.pe_dir, app.ARGS.pf)

    # convert input data to .mif format and concatenate
    convert_input_data(dwi_metadata)

    # get a table of b-shells, echo times, and b-shapes
    shell_table = create_shell_table(dwi_metadata)
    shell_rows = ['b-value', 'b-shape', 'n volumes', 'echo time']
    shell_df = pd.DataFrame(data = shell_table, 
                  index = shell_rows)
    print('input DWI data has properties:')
    print(shell_df)

    app.goto_scratch_dir()
    
    # begin pipeline

    # denoising
    if app.ARGS.denoise:
        if app.ARGS.phase:
            phasepath = path.from_user(app.ARGS.phase)
        else:
            phasepath = None

        # by default perform mppca denoising with optional phase denoising for complex data
        run_mppca(
            app.ARGS.extent, phasepath, app.ARGS.shrinkage, app.ARGS.algorithm)
    
    # patch2self can be run either in addition or alternatively to mppca
    if app.ARGS.patch2self:
        run_patch2self()

    # rpg gibbs artifact correction
    if app.ARGS.degibbs:
        run_degibbs(dwi_metadata['pf'], dwi_metadata['pe_dir'])

    # rigid alignment of b0s from separate input series
    if app.ARGS.pre_align:
        run_pre_align(dwi_metadata)

    # rigid alignment of each dwi volume to the prior volume
    if app.ARGS.ants_motion_correction:
        run_ants_moco(dwi_metadata)

    # eddy current, succeptibility, motion correction
    if app.ARGS.eddy:
        run_eddy(shell_table, dwi_metadata)

    # b1 bias correction
    if app.ARGS.b1correct:
        run_b1correct(dwi_metadata)

    # generate a final brainmask
    if app.ARGS.mask or app.ARGS.normalize:
        create_brainmask(fsl_suffix)    

    # rician bias correction
    if app.ARGS.phase and app.ARGS.rician:
        print('phase included, skipping approximated bias correction')
    elif app.ARGS.rician and not app.ARGS.phase:
        run_rice_bias_correct()

    # normalize separate input series (voxelwise) such that b0 images are properly scaled
    if app.ARGS.normalize:
        run_normalization(dwi_metadata)

    run.command('mrinfo -export_grad_fsl dwi_designer.bvec dwi_designer.bval working.mif', show=False)
    run.command('mrconvert -force -datatype float32le working.mif dwi_designer.nii', show=False)

    outpath = path.from_user(app.ARGS.output, True)
    outfname = path.from_user(os.path.splitext(os.path.basename(outpath))[0], True)

    if app.ARGS.datatype:
        run.command('mrconvert -force -datatype %s %s %s' % 
            (app.ARGS.datatype, path.to_scratch('working.mif', True), outpath))
    else:
        run.command('mrconvert -force %s %s' % 
            (path.to_scratch('working.mif', True), outpath))

    shutil.copyfile(path.to_scratch('dwi_designer.bvec',True), outfname + '.bvec')
    shutil.copyfile(path.to_scratch('dwi_designer.bval',True), outfname + '.bval')

def main():
    import mrtrix3
    mrtrix3.execute() #pylint: disable=no-member



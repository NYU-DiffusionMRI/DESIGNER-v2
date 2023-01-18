#!/usr/bin/env python3

import os, sys
import numpy as np, gzip, shutil
__version__ = "0.0.1"

from lib.designer_input_utils import get_input_info, convert_input_data, assert_inputs
from lib.designer_func_wrappers import run_denoising, run_degibbs, run_eddy, run_b1correct, create_brainmask, run_rice_bias_correct, run_normalization

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

    options.add_argument('-pe_dir', metavar=('<phase encoding direction>'), help='Specify the phase encoding direction of the input series (required if using the eddy option). Can be a signed axis number (e.g. -0, 1, +2), an axis designator (e.g. RL, PA, IS), or NIfTI axis codes (e.g. i-, j, k)')
    options.add_argument('-pf', metavar=('<PF factor>'), help='Specify the partial fourier factor (e.g. 7/8, 6/8)')

    mp_options = cmdline.add_argument_group('Options for specifying mp parameters')
    mp_options.add_argument('-p2c', action='store_true', help='Perform dwidenoise')
    mp_options.add_argument('-denoise', action='store_true', help='Perform dwidenoise')
    mp_options.add_argument('-shrinkage',metavar=('<shrink>'),help='specify shrinkage type. Options are "threshold" or "frob"', default='frob')
    mp_options.add_argument('-algorithm',metavar=('<alg>'),help='specify MP algorithm. Options are "veraart","cordero-grande","jespersen"',default='cordero-grande')
    mp_options.add_argument('-extent', metavar=('<size>'), help='Denoising extent. Default is 5,5,5')
    mp_options.add_argument('-phase', metavar=('<image>'), help='Diffusion phases', default=None)
    
    rpe_options = cmdline.add_argument_group('Options for specifying the acquisition phase-encoding design')
    rpe_options.add_argument('-rpe_none', action='store_true', help='Specify that no reversed phase-encoding image data is being provided; eddy will perform eddy current and motion correction only')
    rpe_options.add_argument('-rpe_pair', metavar=('<reverse PE b=0 image>'), help='Specify the reverse phase encoding image')
    rpe_options.add_argument('-rpe_all', metavar=('<reverse PE dwi volume>'), help='Specify that ALL DWIs have been acquired with opposing phase-encoding; this information will be used to perform a recombination of image volumes (each pair of volumes with the same b-vector but different phase encoding directions will be combined together into a single volume). The argument to this option is the set of volumes with reverse phase encoding but the same b-vectors as the input image')
    rpe_options.add_argument('-rpe_header', action='store_true', help='Specify that the phase-encoding information can be found in the image header(s), and that this is the information that the script should use')

def execute(): #pylint: disable=unused-variable
    from mrtrix3 import app, fsl, run, path #pylint: disable=no-name-in-module, import-outside-toplevel
    
    app.make_scratch_dir()
    
    fsl_suffix = fsl.suffix()

    (isdicom, DWIext, bveclist, bvallist, DWInlist, bidslist) = get_input_info(
        app.ARGS.input, app.ARGS.fslbval, app.ARGS.fslbvec, app.ARGS.bids)
    
    (pe_dir, pf, TE, bshape) = assert_inputs(bidslist, app.ARGS.pe_dir, app.ARGS.pf)

    (idxlist) = convert_input_data(
        isdicom, DWIext, bveclist, bvallist, DWInlist)

    print(TE)
    print(bshape)
    import pdb; pdb.set_trace()

    app.goto_scratch_dir()
    
    # denoising
    if app.ARGS.denoise:
        if app.ARGS.phase:
            phasepath = path.from_user(app.ARGS.phase)
        else:
            phasepath = None

        run_denoising(
            app.ARGS.extent, phasepath, app.ARGS.shrinkage, app.ARGS.algorithm)

        if app.ARGS.p2c:
            run_patch2self()


    # rpg gibbs artifact correction
    if app.ARGS.degibbs:
        run_degibbs(pf, pe_dir)

    if app.ARGS.eddy:
        run_eddy(pe_dir, app.ARGS.rpe_pair)

    if app.ARGS.b1correct:
        run_b1correct(DWInlist, idxlist)

    # generate a final brainmask
    if app.ARGS.mask or app.ARGS.normalize:
        create_brainmask(fsl_suffix)    

    # rician bias correction
    if app.ARGS.phase and app.ARGS.rician:
        print('phase included, skipping approximated bias correction')
    elif app.ARGS.rician and not app.ARGS.phase:
        run_rice_bias_correct()

    if app.ARGS.normalize:
        run_normalization(DWInlist, idxlist)

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



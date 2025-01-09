"""
Functions to run designer preprocessing steps

Includes:
run_denoising:

run_degibbs:

run_eddy: 

run_b1correct:

create_brainmask:

run_rice_bias_correct:

run_normalization:
"""

import time
import shutil

def run_mppca(args_extent, args_phase, args_shrinkage, args_algorithm):
    """
    wrapper for complex or magnitude, adatpive or local mppca
    
    """
    from mrtrix3 import run, app
    from ants import image_read, image_write, from_numpy
    import numpy as np

    if app.ARGS.adaptive_patch:
        import lib.mpdenoise_adaptive as mp
    else:
        import lib.mpcomplex as mp

    if args_extent:
        extent = args_extent
        if type(extent) == str:
            extent = extent.split(",")
            extent = [int(i) for i in extent]
    else: 
        extent = [5,5,5]

    run.command('mrconvert dwi.mif -export_grad_mrtrix grad.txt tmp_dwi.nii', show=False)
    nii = image_read('tmp_dwi.nii')
    dwi = nii.numpy()

    # note adaptive patch does not include mpcomplex
    n_cores = app.ARGS.n_cores
    terminal_width = shutil.get_terminal_size().columns
    separator = "=" * terminal_width

    print("\n" + separator)
    print('...denoising...')
    start_time = time.time()

    if app.ARGS.adaptive_patch:
        
        if app.ARGS.adaptive_patch_length:
            adapt_length = int(app.ARGS.adaptive_patch_length)
        else:
            adapt_length = None

        Signal, Sigma, Nparameters = mp.denoise(
            dwi, phase=args_phase, kernel=extent, patchtype='nonlocal', patchsize=adapt_length, shrinkage=args_shrinkage, algorithm=args_algorithm
            )
    else:    
        Signal, Sigma, Nparameters = mp.denoise(
            dwi, phase=args_phase, kernel=extent, shrinkage=args_shrinkage, algorithm=args_algorithm, n_cores=n_cores
            )
    Sigma[np.isnan(Sigma)] = 0

    out = from_numpy(
        Signal, origin=nii.origin, spacing=nii.spacing, direction=nii.direction)
    image_write(out, 'tmp_dwidn.nii')
    out = from_numpy(
        Sigma, origin=nii.origin[:-1], spacing=nii.spacing[:-1], direction=nii.direction[:-1,:])
    image_write(out, 'sigma.nii')
    out = from_numpy(
        Nparameters, origin=nii.origin[:-1], spacing=nii.spacing[:-1], direction=nii.direction[:-1,:])
    image_write(out, 'Npars.nii')

    run.command('mrconvert -grad grad.txt tmp_dwidn.nii dwidn.mif', show=False)
    run.command('mrconvert sigma.nii noisemap.mif', show=False)
    app.cleanup('tmp_dwi.nii')
    app.cleanup('tmp_dwidn.nii')
    app.cleanup('grad.txt')
    
    #run.command('dwidenoise -noise fullnoisemap.mif -estimator Exp2 working.mif dwidn.mif')
    run.command('mrconvert -force dwidn.mif working.mif', show=False)

    # End timer
    end_time = time.time()
    elapsed_time = end_time - start_time

    # Convert elapsed time to hours and minutes
    hours, rem = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(rem, 60)

    print(f"Denoising completed in {int(hours):02} hours, {int(minutes):02} minutes, {int(seconds):02} seconds.")
    print(separator + "\n")

def run_patch2self():
    """
    wrapper for patch2self
    """
    from mrtrix3 import run, app
    from dipy.denoise.patch2self import patch2self
    from ants import image_read, image_write, from_numpy
    import numpy as np

    terminal_width = shutil.get_terminal_size().columns
    separator = "=" * terminal_width

    print("\n" + separator)
    print('...denoising...')
    start_time = time.time()

    run.command('mrconvert -force -export_grad_fsl working.bvec working.bval working.mif working.nii', show=False)
    nii = image_read('working.nii')
    dwi = nii.numpy()

    bvals = np.loadtxt('working.bval')
    dwi_dn = patch2self(dwi, bvals)
    out = from_numpy(
        dwi_dn, origin=nii.origin, spacing=nii.spacing, direction=nii.direction)
    image_write(out, 'working_p2s.nii')
    run.command('mrconvert -force -fslgrad working.bvec working.bval working_p2s.nii working.mif', show=False)

    # End timer
    end_time = time.time()
    elapsed_time = end_time - start_time

    # Convert elapsed time to hours and minutes
    hours, rem = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(rem, 60)

    print(f"Denoising completed in {int(hours):02} hours, {int(minutes):02} minutes, {int(seconds):02} seconds.")
    print(separator + "\n")

def get_phase_encoding_dimension(pe_dir, strides):
    # Define mappings for 'i', 'j', 'k' to their respective anatomical axes
    axis_map = {
        'i': 0,  # Left-Right (X-axis)
        'j': 1,  # Anterior-Posterior (Y-axis)
        'k': 2   # Superior-Inferior (Z-axis)
    }
    
    # Check for valid phase encoding direction
    if pe_dir not in axis_map:
        raise ValueError("Invalid phase encoding direction. Choose from 'i', 'j', or 'k'.")
    
    # Map 'i', 'j', or 'k' to the target anatomical axis
    target_axis = axis_map[pe_dir]
    
    # Determine which dimension in strides corresponds to the target axis
    try:
        dimension = strides.index(target_axis + 1)  # +1 since strides use 1-based indexing
    except ValueError:
        raise ValueError("Strides configuration does not contain the expected axis order.")
    
    return dimension

def run_degibbs(pf, pe_dir):
    """
    wrapper for rpg degibbs
    """

    import lib.rpg as rpg
    from mrtrix3 import run, app, MRtrixError, image
    from ants import image_read, image_write, from_numpy
    import os
    import numpy as np

    # rpg_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'rpg_cpp')

    # convert working.mif to nii
    run.command('mrconvert -force -export_grad_fsl working.bvec working.bval working.mif working.nii', show=False)
    nii = image_read('working.nii')
    dwi = nii.numpy()

    terminal_width = shutil.get_terminal_size().columns
    separator = "=" * terminal_width

    print("\n" + separator)
    print('...RPG degibbsing...')
    start_time = time.time()
    
    n_cores = app.ARGS.n_cores

    print('Gibbs correction with parameters:')
    #print('phase-encoding direction = %s' % (pe_dir))
    #print('partial-fourier factor   = %s' % (pf))

    # grab original strides
    UserCpath = app.ARGS.input.rsplit(',')
    dwi0 = [os.path.realpath(i) for i in UserCpath][0]
    strides_orig = image.Header(dwi0).strides()
    orient_orig = [abs(i) for i in strides_orig]

    if 'i' in pe_dir:
        pe_dir = 0
    elif 'j' in pe_dir:
        pe_dir = 1
    elif 'k' in pe_dir:
        pe_dir = 2

    pe_dir_orig = (np.array(orient_orig)-1)[pe_dir]
    if pe_dir_orig == 2:
        raise ValueError("PE dir=k should not be possible. Phase encoding direction must be along the first or second axis of the image.")

    # we need to under the -1,2,3,4 striding and use the original PE direction here
    # unring expects a transposed image along x,y (y along dim 0)
    transpose_order = np.argsort(orient_orig)
    dwi_t = np.ascontiguousarray(dwi.transpose(transpose_order).transpose(3,2,1,0))
    dwi_dg_t = rpg.unring(dwi_t, minW=1, maxW=3, nsh=20, pfv=float(pf), pfdimf=pe_dir_orig, phase_flag=False)
    dwi_dg = dwi_dg_t[0].copy().transpose(3,2,1,0).transpose(np.argsort(transpose_order))
    
    out = from_numpy(
        dwi_dg, origin=nii.origin, spacing=nii.spacing, direction=nii.direction)
    image_write(out, 'working_rpg.nii')

    #convert gibbs corrected nii to .mif
    run.command('mrconvert -force -fslgrad working.bvec working.bval working_rpg.nii working.mif', show=False)

    # End timer
    end_time = time.time()
    elapsed_time = end_time - start_time

    # Convert elapsed time to hours and minutes
    hours, rem = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(rem, 60)

    print(f"RPG gibbs correction completed in {int(hours):02} hours, {int(minutes):02} minutes, {int(seconds):02} seconds.")
    print(separator + "\n")

def group_alignment(group_list):
    from mrtrix3 import app, run, path, image, fsl, MRtrixError
    import numpy as np

    fsl_suffix = fsl.suffix()

    terminal_width = shutil.get_terminal_size().columns
    separator = "=" * terminal_width

    print("\n" + separator)
    print('... Rigidly aligning groups...')
    start_time = time.time()

    # 1) find the longest series and extract the b0
    group_len = []
    for i in range(len(group_list)):
        std = run.command('mrinfo -size "%s"' %
                          (group_list[i]
                           ), show=False)
        group_len.append([int(j) for j in std[0].strip().split(' ')][-1])
    longest_series = np.argmax(group_len)

    # idxlist = dwi_metadata['idxlist']
    # series_inds = [i.split(",") for i in idxlist]
    # series_lens = [len(i) for i in series_inds]
    # longest_series = np.argmax(series_lens)

    # For each group: Register b0 ref PE and b0 RPE to the b0 and run topup + eddy
    for i in range(len(group_list)):
        run.command('dwiextract -force -bzero "%s" - | mrmath - mean "%s" -axis 3 -force' %
                    (group_list[i], 
                    'b0_align_series_' + str(i) + '.nii'),
                    show=False)
        run.command('bet "%s" "%s" -f 0.2 -m' %
                    ('b0_align_series_' + str(i) + '.nii', 
                    'b0_align_series_' + str(i) + '_brain'),
                    show=False)

    # 2) extract b0 images from all other series
    for i in range(len(group_list)):                
        if i == longest_series:
            continue
        else:
            # 3) rigidly register all b0 images to the b0 from the longest series
            run.command('flirt -in "%s" -ref "%s" -omat "%s" -dof 6' %
                    ('b0_align_series_' + str(i) + '_brain',
                    'b0_align_series_' + str(longest_series) + '_brain',
                    'b0_align_series_to_longest_series_' + str(i) + '.mat'),
                    show=False)
            run.command('transformconvert -force "%s" "%s" "%s" flirt_import "%s"' %
                ('b0_align_series_to_longest_series_' + str(i) + '.mat',
                'b0_align_series_' + str(i) + '_brain' + fsl_suffix,
                'b0_align_series_' + str(longest_series) + '_brain' + fsl_suffix,
                'b0_align_series_scan' + str(i) + 'to_longest_series_mrtrix.txt'),
                show=False)
            run.command('mrtransform -force -linear "%s" -interp cubic "%s" "%s"' %
                ('b0_align_series_scan' + str(i) + 'to_longest_series_mrtrix.txt',
                group_list[i],
                'dwi_align_series_' + str(i) + '_to_ref.mif'),
                show=False)
    
    # concatenate data in correct order
    series_cat_all = []
    for i in range(len(group_list)):
        if i == longest_series:
            series_cat_all.append(group_list[i])
        else:
            series_to_cat = 'dwi_align_series_' + str(i) + '_to_ref.mif'
            series_cat_all.append(series_to_cat)
    SList = ' '.join(series_cat_all)
    run.command('mrcat -force %s %s' %
                (SList,
                'dwi_align_series_aligned.mif'),
                show=False)
    run.command('mrconvert -force dwi_align_series_aligned.mif working.mif', show=False)

    # End timer
    end_time = time.time()
    elapsed_time = end_time - start_time

    # Convert elapsed time to hours and minutes
    hours, rem = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(rem, 60)

    print(f"rigid alignment completed in {int(hours):02} hours, {int(minutes):02} minutes, {int(seconds):02} seconds.")
    print(separator + "\n")


def convert_ants_xform(mat, i):
    import scipy.io as sio
    import numpy as np
    from mrtrix3 import app, run, path, image, fsl, MRtrixError

    ants_dict = sio.loadmat(mat)
    lps2ras = np.diag([-1, -1, 1])
    rot = ants_dict['AffineTransform_float_3_3'][0:9].reshape((3, 3))
    trans = ants_dict['AffineTransform_float_3_3'][9:12]
    offset = ants_dict['fixed']
    r_trans = (np.dot(rot, offset) - offset - trans).T * [1, 1, -1]
    
    data = np.eye(4)
    data[0:3, 3] = r_trans
    data[:3, :3] = np.dot(np.dot(lps2ras, rot), lps2ras)

    #np.savetxt('ants_textfile_' + str(i) + '_np.txt', data, fmt='%10.10f', delimiter=' ')
    return data

def pre_eddy_ants_moco(dwi_metadata):
    import ants
    from mrtrix3 import app, run, path, image, fsl, MRtrixError
    import numpy as np

    terminal_width = shutil.get_terminal_size().columns
    separator = "=" * terminal_width

    print("\n" + separator)
    print('... ANTS motion correction...')
    start_time = time.time()
    
    run.command('mrconvert -force -export_grad_fsl working.bvec working.bval working.mif working.nii', show=False)
    nii = ants.image_read('working.nii')
    ants_moco_params = ants.motion_correction(nii, moreaccurate=2)
    ants.image_write(ants_moco_params['motion_corrected'], 'working_antsmoco.nii')

    # ants_interp = ants.apply_transforms(fixed=nii, moving=nii, 
    #                                     transformlist=ants_moco_params['motion_parameters'],
    #                                     interpolator='bSpline', imagetype=3)

    dirs = dwi_metadata['grad'][:,:3]
    dirs_rot = np.zeros_like(dirs)
    # rotate gradients
    #img_cat_all = ''
    for i in range(dwi_metadata['grad'].shape[0]):

        antsaffine = ants_moco_params['motion_parameters'][i][0]
        aff = convert_ants_xform(antsaffine, i)
        diri = np.hstack((dirs[i,:],0))
        dirs_rot[i,:] = (aff @ diri.T)[:3]

    # run.command('mrconvert working.mif -coord 3 %s -axes 0,1,2 dwi_antsmoco_%s.mif' % 
    #             (str(i), str(i)), show=False
    #             )
    #     run.command('mrtransform -linear %s %s %s' %
    #                 ('ants_textfile_' + str(i) + '_np.txt',
    #                 'dwi_antsmoco_' + str(i) + '.mif',
    #                 'dwi_antsmoco_' + str(i) + '_withgrad.mif'),
    #                 show=False)
    #     series_to_cat = path.to_scratch('dwi_antsmoco_' + str(i) + '_withgrad.mif')
    #     img_cat_all = img_cat_all + series_to_cat + ' '
    # run.command('mrcat %s %s' %
    #             (img_cat_all,
    #             path.to_scratch('dwi_pre_eddy_antsmoco.mif')))

    np.savetxt('working_antsmoco.bvec', dirs_rot.T, delimiter=' ', fmt='%4.10f')
    run.command('mrconvert -force -fslgrad working_antsmoco.bvec working.bval working_antsmoco.nii working.mif', show=False)

    # End timer
    end_time = time.time()
    elapsed_time = end_time - start_time

    # Convert elapsed time to hours and minutes
    hours, rem = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(rem, 60)

    print(f"motion correction completed in {int(hours):02} hours, {int(minutes):02} minutes, {int(seconds):02} seconds.")
    print(separator + "\n")

def run_pre_align(dwi_metadata):
    from mrtrix3 import app, run, path, image, fsl, MRtrixError

    if len(dwi_metadata['dwi_list']) == 1:
        print('Warning -pre_align option was used but only 1 input series was provided')
    else:
        series_list = []
        for i in range(len(dwi_metadata['dwi_list'])):
            run.command('mrconvert -coord 3 "%s" working.mif "%s"' % 
                    (dwi_metadata['idxlist'][i], 'dwi_align_series_' + str(i) + '.mif'),
                    show=False)
            series_list.append('dwi_align_series_' + str(i) + '.mif')
            
        group_alignment(series_list)

def run_ants_moco(dwi_metadata):
    from mrtrix3 import app, run, path, image, fsl, MRtrixError

    pre_eddy_ants_moco(dwi_metadata)

def run_eddy(shell_table, dwi_metadata):
    from mrtrix3 import app, run, path, image, fsl, MRtrixError
    from lib.designer_input_utils import splitext_
    import numpy as np
    import json
    import os
    import glob

    terminal_width = shutil.get_terminal_size().columns
    separator = "=" * terminal_width

    print("\n" + separator)
    print('... Eddy current, EPI, motion correction...')
    start_time = time.time()

    run.command('mrconvert -force -export_grad_fsl working.bvec working.bval working.mif working.nii', show=False)
    run.command('mrconvert -force -fslgrad working.bvec working.bval -json_import working.json working.nii working.mif', show=False)

    eddyopts = '" --cnr_maps --repol --data_is_shelled "'

    fsl_suffix = fsl.suffix()

    pe_dir = dwi_metadata['pe_dir']

    echo_times_per_series = dwi_metadata['TE']
    bshape_per_series = dwi_metadata['bshape']
    if  echo_times_per_series.count(echo_times_per_series[0]) == len(echo_times_per_series):
        flag_variable_TE = False
    else:
        flag_variable_TE = True

    if  bshape_per_series.count(bshape_per_series[0]) == len(bshape_per_series):
        flag_variable_bshape = False
    else:
        flag_variable_bshape = True

    if flag_variable_TE or flag_variable_bshape:
        if app.ARGS.eddy_groups is None and app.ARGS.eddy_fakeb is None:
            raise MRtrixError('for variable TE data, the user must supply the \
                            -eddy_groups or -eddy_fakeb command line argument')

        # scale b-value by the fake_b scaling factors
        if app.ARGS.eddy_fakeb is not None:
            eddy_fakeb = [float(i) for i in app.ARGS.eddy_fakeb.rsplit(',')]
            if len(dwi_metadata['bvals']) != len(eddy_fakeb):
                raise ValueError("The number of bval directories does not match the length of the eddy_fakeb array.")

            fake_bvals = np.array([])
            for i, bval_path in enumerate(dwi_metadata['bvals']):
                bvals = np.loadtxt(bval_path)
                scaled_bvals = bvals * eddy_fakeb[i]
                fake_bvals = np.concatenate((fake_bvals, scaled_bvals))

            all_concatenated_bvecs = np.empty((0, 3))
            for bvec_path in dwi_metadata['bvecs']:
                bvecs = np.transpose(np.loadtxt(bvec_path))
                all_concatenated_bvecs = np.vstack((all_concatenated_bvecs, bvecs))

            fakeb_grad = np.hstack((all_concatenated_bvecs, fake_bvals[:, np.newaxis]))
            np.savetxt('fakeb_grad.txt', fakeb_grad, fmt='%0.8f')
            print("Scaled b-values and concatenated b-vectors have been saved to fakeb_grad.txt")
        
         # get TE of the PA
        if app.ARGS.rpe_pair:
            if app.ARGS.rpe_te:
                te_pa = float(app.ARGS.rpe_te)
            else:
                pa_fpath = splitext_(app.ARGS.rpe_pair)[0]
                pa_bids_path = pa_fpath + '.json'
                if not os.path.exists(pa_bids_path):
                    raise MRtrixError('for variable TE data the RPE \
                        image must be accompanies by a bids .json or use\
                        must use the -rpe_te argument to specify the PA echo time')

                pa_bids = json.load(open(pa_bids_path))
                te_pa = float(pa_bids['echo_time'])

            if (te_pa > 1) and (min(echo_times_per_series) < 1):
                te_pa /= 1000
            elif (te_pa < 1) and (min(echo_times_per_series) > 1):
                te_pa *= 1000

            id_dwi_match_pa = np.where(te_pa == np.array(echo_times_per_series))[0][0]
            if id_dwi_match_pa.size == 0:
                raise MRtrixError('echo time of reverse phase encoding image does not\
                                  match any of the input echo times, please check.')

            bidslist = dwi_metadata['bidslist']
            rpe_fpath = splitext_(app.ARGS.rpe_pair)[0]
            rpe_bids_path = rpe_fpath + '.json'
            rpe_bvals_path = rpe_fpath + '.bval'
            rpe_bvec_path = rpe_fpath + '.bvec'

            if os.path.exists(bidslist[0]) and os.path.exists(rpe_bids_path):
                # extract b0s of dwi matching te of PA image
                run.command('mrconvert -coord 3 %s working.mif - | dwiextract -bzero - - | mrmath - mean "%s" -axis 3' %
                    (dwi_metadata['idxlist'][id_dwi_match_pa],
                    'pe_original_meanb0.mif'),
                    show=False)

                rpe_size = [ int(s) for s in image.Header(app.ARGS.rpe_pair).size() ]
                if len(rpe_size) == 4:
                    run.command('mrconvert -coord 3 0 -strides -1,+2,+3 -json_import "%s" "%s" "%s"' % 
                        (rpe_bids_path, app.ARGS.rpe_pair, 'rpe_b0.mif'))
                else: 
                    run.command('mrconvert -strides -1,+2,+3 -json_import "%s" "%s" "%s"' % 
                        (rpe_bids_path, app.ARGS.rpe_pair, 'rpe_b0.mif'))
                run.command('mrconvert pe_original_meanb0.mif pe_original_meanb0.nii')
                run.command('mrconvert rpe_b0.mif rpe_b0.nii')
            else:
                # extract b0s of dwi matching te of PA image
                run.command('mrconvert -coord 3 %s working.mif - | dwiextract -bzero - - | mrmath - mean "%s" -axis 3' %
                    (dwi_metadata['idxlist'][id_dwi_match_pa],
                    'pe_original_meanb0.nii'),
                    show=False)

                rpe_size = [ int(s) for s in image.Header(app.ARGS.rpe_pair).size() ]
                if len(rpe_size) == 4:
                    run.command('mrconvert -coord 3 0 -strides -1,+2,+3 "%s" "%s"' % 
                        (app.ARGS.rpe_pair, 'rpe_b0.nii'))
                else: 
                    run.command('mrconvert -strides -1,+2,+3 "%s" "%s"' % 
                        (app.ARGS.rpe_pair, 'rpe_b0.nii'))
            
            # extract brain from mean b0
            run.command('bet "%s" "%s" -f 0.2 -m' %
                ('pe_original_meanb0.nii', 
                'pe_original_brain'))
            
            run.command('bet "%s" "%s" -f 0.20' % 
                ('rpe_b0.nii', 'rpe_b0_brain'))


    if app.ARGS.eddy_groups:
        eddy_groups = [int(i) for i in app.ARGS.eddy_groups.rsplit(',')]
        ngroups = max(eddy_groups)
        group_ids = np.arange(1,ngroups+1,1)

        for i in iter(group_ids):
            group_idx = np.where(np.array(eddy_groups) == i)[0]
            volume_idx_list = [dwi_metadata['idxlist'][j] for j in group_idx]
            volume_idx = ','.join(volume_idx_list)

            # get the indices for each series in the group
            group_vol_inds = [list(map(int, j.split(','))) for j in volume_idx_list]
            start = 0
            ginds = []
            for j in range(len(group_vol_inds)):
                gindsi = np.arange(start, start + len(group_vol_inds[j]))
                gindsi_str = ','.join([str(k) for k in gindsi])
                ginds.append(gindsi_str)
                start += len(group_vol_inds[j])            
           
            run.command('mrconvert -coord 3 "%s" working.mif "%s"' % 
                (volume_idx,
                'dwi_pre_eddy_' + str(i) + '.mif'),
                show=False)
            
            if app.ARGS.rpe_pair:
                run.command('dwiextract -bzero "%s" - | mrmath - mean "%s" -axis 3' %
                    ('dwi_pre_eddy_' + str(i) + '.mif', 
                    'b0_pre_eddy_' + str(i) + '.nii'))
                run.command('bet "%s" "%s" -f 0.2 -m' %
                    ('b0_pre_eddy_' + str(i) + '.nii', 
                    'b0_pre_eddy_' + str(i) + '_brain'))

                run.command('flirt -in "%s" -ref "%s" -out "%s" -dof 6' %
                    ('rpe_b0_brain',
                    'b0_pre_eddy_' + str(i) + '_brain',
                    'rpe_to_ref_' + str(i)))

                run.command('flirt -in "%s" -ref "%s" -out "%s" -dof 6' %
                    ('pe_original_brain' + fsl_suffix, 
                    'b0_pre_eddy_' + str(i) + '_brain', 
                    'pe_to_ref_' + str(i)))

                run.command('mrcat -force -axis 3 "%s" "%s" "%s"' %
                    ('pe_to_ref_' + str(i) + fsl_suffix,
                    'rpe_to_ref_' + str(i) + fsl_suffix,
                    'b0_pair_topup_' + str(i) + '.nii'))

                if os.path.exists(bidslist[0]) and os.path.exists(rpe_bids_path):
                    run.command('mrinfo pe_original_meanb0.mif -export_pe_eddy topup_config_1.txt topup_indicies_1.txt')
                    run.command('mrinfo rpe_b0.mif -export_pe_eddy topup_config_2.txt topup_indicies_2.txt')
                    filenames = ['topup_config_1.txt', 'topup_config_2.txt']
                    with open('topup_acqp.txt', 'w') as outfile:
                        for fname in filenames:
                            with open(fname) as infile:
                                outfile.write(infile.read())
                else:
                    acqp = np.zeros((2,3))
                    if 'i' in pe_dir: acqp[:,0] = 1
                    if 'j' in pe_dir: acqp[:,1] = 1
                    if 'k' in pe_dir: acqp[:,2] = 1
                    if '-' in pe_dir:
                        acqp[0,:] = -acqp[0,:]
                    else:
                        acqp[1,:] = -acqp[1,:]
                    
                    acqp[acqp==-0] = 0
                    acqp = np.hstack((acqp, np.array([0.1,0.1])[...,None]))
                    np.savetxt('topup_acqp.txt', acqp, fmt="%1.2f")

                # if any of the image dims are odd dont subsample during topup
                odd_dims = [ int(s) for s in image.Header('pe_to_ref_' + str(i) + fsl_suffix).size()[:3] if s % 2 ]
                if np.any(np.array(odd_dims)):
                    flag_no_subsampling = True
                else:
                    flag_no_subsampling = False

                if flag_no_subsampling:
                    run.command('topup --imain="%s" --datain="%s" --config=b02b0.cnf --subsamp=1 --scale=1 --out="%s" --iout="%s"' %
                        ('b0_pair_topup_' + str(i) + '.nii',
                        'topup_acqp.txt',
                        'topup_results_' + str(i),
                        'topup_results_' + str(i) + fsl_suffix))
                else:
                    run.command('topup --imain="%s" --datain="%s" --config=b02b0.cnf --scale=1 --out="%s" --iout="%s"' %
                        ('b0_pair_topup_' + str(i) + '.nii',
                        'topup_acqp.txt',
                        'topup_results_' + str(i),
                        'topup_results_' + str(i) + fsl_suffix))
                
                # mask the topup corrected image
                run.command('mrmath "%s" mean "%s" -axis 3' %
                            ('topup_results_' + str(i) + fsl_suffix,
                             'topup_corrected_' + str(i) + '_mean.nii'
                            ))
                 
                run.command('bet "%s" "%s" -f 0.2 -m' %
                    ('topup_corrected_' + str(i) + '_mean.nii', 
                    'topup_corrected_' + str(i) + '_brain'))
                
                run.command('dwifslpreproc -nocleanup -scratch "%s" -eddy_options "%s" -rpe_none -eddy_mask "%s" -topup_files "%s" -pe_dir "%s" "%s" "%s"' % 
                    ('eddy_processing_' + str(i), 
                    eddyopts, 
                    'topup_corrected_' + str(i) + '_brain_mask' + fsl_suffix,
                    'topup_results_' + str(i),
                    pe_dir,
                    'dwi_pre_eddy_' + str(i) + '.mif',
                    'dwi_post_eddy_' + str(i) + '.mif'))
                
            elif app.ARGS.rpe_none:

                run.command('dwifslpreproc -nocleanup -scratch "%s" -eddy_options "%s" -rpe_none -pe_dir "%s" "%s" "%s"' % 
                    ('eddy_processing_' + str(i), 
                    eddyopts, 
                    pe_dir,
                    'dwi_pre_eddy_' + str(i) + '.mif',
                    'dwi_post_eddy_' + str(i) + '.mif'))
                
            elif app.ARGS.rpe_all:
                ### this needs to be changed for realistic compatiblity to eddy_groups
                run.command('mrconvert -export_grad_mrtrix grad.txt working.mif tmp.mif', show=False)
                run.command('mrconvert -strides -1,+2,+3,+4 -grad grad.txt ' + app.ARGS.rpe_all + ' dwirpe.mif', show=False)
                run.command('mrcat -axis 3 working.mif dwirpe.mif dwipe_rpe.mif')
                run.command('dwifslpreproc -nocleanup -scratch "%s" -eddy_options "%s" -rpe_all -pe_dir "%s" "%s" "%s"' %
                    ('eddy_processing',
                     eddyopts, 
                     pe_dir, 
                     'dwi_pre_eddy_' + str(i) + '.mif',
                     'dwi_post_eddy_' + str(i) + '.mif'))
                run.function(os.remove,'tmp.mif')

            elif app.ARGS.rpe_header:

                run.command('dwifslpreproc -nocleanup -scratch "%s" -eddy_options "%s" -rpe_header "%s" "%s"' % 
                ('eddy_processing',
                 eddyopts,
                'dwi_pre_eddy_' + str(i) + '.mif',
                'dwi_post_eddy_' + str(i) + '.mif'))

            run.command('dwiextract -bzero "%s" - | mrmath - mean "%s" -axis 3' %
                    ('dwi_post_eddy_' + str(i) + '.mif', 
                    'b0_post_eddy_' + str(i) + '.nii'))

            # undo grouping:
            for j in range(len(group_idx)):
                run.command('mrconvert -coord 3 "%s" "%s" "%s"' %
                            (ginds[j],
                            'dwi_post_eddy_' + str(i) + '.mif',
                            'dwi_post_eddy_ungrouped_' + str(group_idx[j]) + '.mif'),
                            show=False)
                
        
        registered_data = glob.glob(os.path.join(path.to_scratch(''),'dwi_post_eddy_ungrouped*.mif'))
        registered_data.sort()

        group_alignment(registered_data)
        run.command('mrconvert -force working.mif dwiec.mif', show=False)
        

    # if not variable TE (not using eddy_groups - only run eddy once)
    else:

        if app.ARGS.rpe_pair:
            bidslist = dwi_metadata['bidslist']
            rpe_fpath = splitext_(app.ARGS.rpe_pair)[0]
            rpe_bids_path = rpe_fpath + '.json'
            rpe_bvals_path = rpe_fpath + '.bval'
            rpe_bvec_path = rpe_fpath + '.bvec'
            
            if os.path.exists(bidslist[0]) and os.path.exists(rpe_bids_path):
                run.command('dwiextract -bzero dwi.mif - | mrconvert -coord 3 0 - b0pe.mif')
                rpe_size = [ int(s) for s in image.Header(app.ARGS.rpe_pair).size() ]
                if len(rpe_size) == 4:
                    run.command('mrconvert "%s" -coord 3 0 -strides -1,+2,+3 -json_import "%s" b0rpe.mif' % 
                            (app.ARGS.rpe_pair, rpe_bids_path))
                else: 
                    run.command('mrconvert "%s" -strides -1,+2,+3 -json_import "%s" b0rpe.mif' % 
                            (app.ARGS.rpe_pair, rpe_bids_path))
                
                run.command('mrinfo b0pe.mif -export_pe_eddy topup_config_1.txt topup_indicies_1.txt')
                run.command('mrinfo b0rpe.mif -export_pe_eddy topup_config_2.txt topup_indicies_2.txt')
                filenames = ['topup_config_1.txt', 'topup_config_2.txt']
                with open('topup_acqp.txt', 'w') as outfile:
                    for fname in filenames:
                        with open(fname) as infile:
                            outfile.write(infile.read())

                run.command('mrconvert b0rpe.mif b0rpe.nii')
                run.command('mrconvert b0pe.mif b0pe.nii')
            else:
                run.command('dwiextract -bzero dwi.mif - | mrconvert -coord 3 0 - b0pe.nii')
                rpe_size = [ int(s) for s in image.Header(app.ARGS.rpe_pair).size() ]
                if len(rpe_size) == 4:
                    run.command('mrconvert "%s" -strides -1,+2,+3 -coord 3 0 b0rpe.nii' % (app.ARGS.rpe_pair))
                else: 
                    run.command('mrconvert -strides -1,+2,+3 "%s" b0rpe.nii' % (app.ARGS.rpe_pair))

                acqp = np.zeros((2,3))
                if 'i' in pe_dir: acqp[:,0] = 1
                if 'j' in pe_dir: acqp[:,1] = 1
                if 'k' in pe_dir: acqp[:,2] = 1
                if '-' in pe_dir:
                    acqp[0,:] = -acqp[0,:]
                else:
                    acqp[1,:] = -acqp[1,:]
                    
                acqp[acqp==-0] = 0
                acqp = np.hstack((acqp, np.array([0.1,0.1])[...,None]))
                np.savetxt('topup_acqp.txt', acqp, fmt="%1.2f")

            run.command('flirt -in b0rpe.nii -ref b0pe.nii -dof 6 -out b0rpe2pe.nii.gz')
            run.command('mrcat -axis 3 b0pe.nii b0rpe2pe.nii.gz b0_pair_topup.nii')

            # if any of the image dims are odd dont subsample during topup. might be better off changing this to padding so topup doesnt take forever
            odd_dims = [ int(s) for s in image.Header('b0pe.nii').size()[:3] if s % 2 ]
            if np.any(np.array(odd_dims)):
                flag_no_subsampling = True
            else:
                flag_no_subsampling = False

            if flag_no_subsampling:
                run.command('topup --imain="%s" --datain="%s" --config=b02b0.cnf --subsamp=1 --scale=1 --out="%s" --iout="%s"' %
                    ('b0_pair_topup.nii',
                    'topup_acqp.txt',
                    'topup_results',
                    'topup_results' + fsl_suffix))
            else:
                run.command('topup --imain="%s" --datain="%s" --config=b02b0.cnf --scale=1 --out="%s" --iout="%s"' %
                    ('b0_pair_topup.nii',
                    'topup_acqp.txt',
                    'topup_results',
                    'topup_results' + fsl_suffix))
                
            # mask the topup corrected image
            run.command('mrmath %s mean %s -axis 3' %
                        ('topup_results' + fsl_suffix,
                            'topup_corrected_mean.nii'
                        ))
                 
            run.command('bet %s %s -f 0.2 -m' %
                ('topup_corrected_mean.nii', 
                'topup_corrected_brain'))
            
            if app.ARGS.eddy_fakeb is None:
                run.command('dwifslpreproc -nocleanup -scratch "%s" -eddy_options %s -rpe_none -eddy_mask "%s" -topup_files "%s" -pe_dir "%s" "%s" "%s"' % 
                        ('eddy_processing', 
                        eddyopts, 
                        'topup_corrected_brain_mask' + fsl_suffix,
                        'topup_results',
                        pe_dir,
                        'working.mif',
                        'dwiec.mif'))
            else:
                run.command('dwifslpreproc -nocleanup -scratch "%s" -grad "%s" -eddy_options %s -rpe_none -eddy_mask "%s" -topup_files "%s" -pe_dir "%s" "%s" "%s"' % 
                        ('eddy_processing',
                        'fakeb_grad.txt', 
                        eddyopts, 
                        'topup_corrected_brain_mask' + fsl_suffix,
                        'topup_results',
                        pe_dir,
                        'working.mif',
                        'dwiec.mif'))

        elif app.ARGS.rpe_none:

            run.command('dwiextract -bzero dwi.mif - | mrconvert -coord 3 0 - b0pe.nii')
            run.command('bet %s %s -f 0.2 -m' %
                ('b0pe.nii', 
                'b0_pe_brain'))
            
            if app.ARGS.eddy_fakeb is None:
                run.command('dwifslpreproc -nocleanup -scratch "%s" -eddy_options %s -rpe_none -eddy_mask "%s" -pe_dir "%s" "%s" "%s"' % 
                        ('eddy_processing', 
                        eddyopts, 
                        'b0_pe_brain_mask' + fsl_suffix,
                        pe_dir,
                        'working.mif',
                        'dwiec.mif'))
            else:
                run.command('dwifslpreproc -nocleanup -scratch "%s" -grad "%s" -eddy_options %s -rpe_none -eddy_mask "%s" -pe_dir "%s" "%s" "%s"' % 
                        ('eddy_processing', 
                        'fakeb_grad.txt', 
                        eddyopts, 
                        'b0_pe_brain_mask' + fsl_suffix,
                        pe_dir,
                        'working.mif',
                        'dwiec.mif'))
            
        elif app.ARGS.rpe_all:
            # run an initial topup to create a brain mask
            run.command('mrconvert -export_grad_mrtrix grad.txt dwi.mif tmp.mif', show=False)
            run.command('dwiextract -bzero working.mif - | mrconvert -coord 3 0 - b0pe.mif')

            bidslist = dwi_metadata['bidslist']
            rpe_fpath = splitext_(app.ARGS.rpe_all)[0]
            rpe_bids_path = rpe_fpath + '.json'
            
            if os.path.exists(bidslist[0]) and os.path.exists(rpe_bids_path):
                run.command('mrconvert -grad grad.txt -strides -1,+2,+3,+4 -json_import "%s" "%s" dwirpe.mif' % (rpe_bids_path,app.ARGS.rpe_all))
                run.command('dwiextract -bzero dwirpe.mif - | mrconvert -coord 3 0 - b0rpe.mif')
                run.command('mrinfo b0pe.mif -export_pe_eddy topup_config_1.txt topup_indicies_1.txt')
                run.command('mrinfo b0rpe.mif -export_pe_eddy topup_config_2.txt topup_indicies_2.txt')
                filenames = ['topup_config_1.txt', 'topup_config_2.txt']
                with open('topup_acqp.txt', 'w') as outfile:
                    for fname in filenames:
                        with open(fname) as infile:
                            outfile.write(infile.read())
                run.command('mrconvert b0rpe.mif b0rpe.nii')
            else:
                acqp = np.zeros((2,3))
                if 'i' in pe_dir: acqp[:,0] = 1
                if 'j' in pe_dir: acqp[:,1] = 1
                if 'k' in pe_dir: acqp[:,2] = 1
                if '-' in pe_dir:
                    acqp[0,:] = -acqp[0,:]
                else:
                    acqp[1,:] = -acqp[1,:]
                    
                acqp[acqp==-0] = 0
                acqp = np.hstack((acqp, np.array([0.1,0.1])[...,None]))
                np.savetxt('topup_acqp.txt', acqp, fmt="%1.2f")

            run.command('mrconvert -strides -1,+2,+3,+4 -grad grad.txt "%s" dwirpe.mif' % (app.ARGS.rpe_all))
            run.command('dwiextract -bzero dwirpe.mif - | mrconvert -coord 3 0 - b0rpe.nii')
            run.command('mrconvert b0pe.mif b0pe.nii')
            run.command('flirt -in b0rpe.nii -ref b0pe.nii -dof 6 -out b0rpe2pe.nii.gz')
            run.command('mrcat -axis 3 b0pe.mif b0rpe2pe.nii.gz b0_pair_topup.nii')

            # if any of the image dims are odd dont subsample during topup. might be better off changing this to padding so topup doesnt take forever
            odd_dims = [ int(s) for s in image.Header('b0pe.nii').size()[:3] if s % 2 ]
            if np.any(np.array(odd_dims)):
                flag_no_subsampling = True
            else:
                flag_no_subsampling = False

            if flag_no_subsampling:
                run.command('topup --imain="%s" --datain="%s" --config=b02b0.cnf --subsamp=1 --scale=1 --out="%s" --iout="%s"' %
                    ('b0_pair_topup.nii',
                    'topup_acqp.txt',
                    'topup_results',
                    'topup_results' + fsl_suffix))
            else:
                run.command('topup --imain="%s" --datain="%s" --config=b02b0.cnf --scale=1 --out="%s" --iout="%s"' %
                    ('b0_pair_topup.nii',
                    'topup_acqp.txt',
                    'topup_results',
                    'topup_results' + fsl_suffix))
                
            # mask the topup corrected image
            run.command('mrmath "%s" mean "%s" -axis 3' %
                        ('topup_results' + fsl_suffix,
                            'topup_corrected_mean.nii'
                        ))
                 
            run.command('bet "%s" "%s" -f 0.2 -m' %
                ('topup_corrected_mean.nii', 
                'topup_corrected_brain'))

            run.command('mrcat -axis 3 working.mif dwirpe.mif dwipe_rpe.mif')
            
            if app.ARGS.eddy_fakeb is None:
                run.command('dwifslpreproc -nocleanup -scratch "%s" -eddy_options %s -rpe_all -pe_dir "%s" -eddy_mask "%s" dwipe_rpe.mif dwiec.mif' %
                            ('eddy_processing',
                            eddyopts, 
                            pe_dir,
                            'topup_corrected_brain_mask' + fsl_suffix
                            ))
            else:
                run.command('dwifslpreproc -nocleanup -scratch "%s" -grad "%s" -eddy_options %s -rpe_all -pe_dir "%s" -eddy_mask "%s" dwipe_rpe.mif dwiec.mif' %
                            ('eddy_processing',
                            'fakeb_grad.txt',
                            eddyopts, 
                            pe_dir,
                            'topup_corrected_brain_mask' + fsl_suffix
                            ))
            run.function(os.remove,'tmp.mif')

        elif app.ARGS.rpe_header:
            if app.ARGS.eddy_fakeb is None: 
                cmd = ('dwifslpreproc -nocleanup -scratch "%s" -eddy_options "%s" -rpe_header working.mif dwiec.mif' % 
                    ('eddy_processing',eddyopts))
            else:
                cmd = ('dwifslpreproc -nocleanup -scratch "%s" -grad "%s" -eddy_options %s -rpe_header working.mif dwiec.mif' % 
                    ('eddy_processing','fakeb_grad.txt', eddyopts))
            run.command(cmd)

        elif not app.ARGS.rpe_header and not app.ARGS.rpe_all and not app.ARGS.rpe_pair:
            raise MRtrixError("the eddy option must run alongside -rpe_header, -rpe_all, or -rpe_pair option")
    if app.ARGS.eddy_quad_output:
        run.command('eddy_quad "%s" -idx "%s" -par "%s" -m "%s" -b "%s" -g "%s" -o "%s"' %
                    ('eddy_processing/dwi_post_eddy',
                    'eddy_processing/eddy_indices.txt',
                    'eddy_processing/eddy_config.txt',
                    'eddy_processing/eddy_mask.nii',
                    'eddy_processing/bvals',
                    'eddy_processing/bvecs',
                    app.ARGS.eddy_quad_output
                    ))
    else:
        run.command('eddy_quad "%s" -idx "%s" -par "%s" -m "%s" -b "%s" -g "%s"' %
                    ('eddy_processing/dwi_post_eddy',
                    'eddy_processing/eddy_indices.txt',
                    'eddy_processing/eddy_config.txt',
                    'eddy_processing/eddy_mask.nii',
                    'eddy_processing/bvals',
                    'eddy_processing/bvecs'
                    ))

    run.command('mrconvert -force -fslgrad working.bvec working.bval dwiec.mif working.mif', show=False)
    # End timer
    end_time = time.time()
    elapsed_time = end_time - start_time

    # Convert elapsed time to hours and minutes
    hours, rem = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(rem, 60)

    print(f"Eddy current, EPI, motion correction completed in {int(hours):02} hours, {int(minutes):02} minutes, {int(seconds):02} seconds.")
    print(separator + "\n")

def run_b1correct(dwi_metadata):
    from mrtrix3 import run

    DWInlist = dwi_metadata['dwi_list']
    idxlist = dwi_metadata['idxlist']

    terminal_width = shutil.get_terminal_size().columns
    separator = "=" * terminal_width

    print("\n" + separator)
    print('... B1 correction...')
    start_time = time.time()

    # b1 bias field correction
    if len(DWInlist) == 1:
        run.command('dwibiascorrect ants -bias biasfield.mif working.mif dwibc.mif')
    else:
        # b1 correction may still need to be done individually for each diffusion series ...
        miflist = []
        for idx,i in enumerate(DWInlist):
            cmd = ('mrconvert -coord 3 "%s" working.mif dwiprebc%s.mif' % 
                (idxlist[idx], str(idx)))
            run.command(cmd)
            cmd = ('dwibiascorrect ants -bias biasfield%s.mif dwiprebc%s.mif dwibc%s.mif' % 
                (str(idx), str(idx), str(idx)))
            run.command(cmd)
            miflist.append('dwibc' + str(idx) + '.mif')
            DWImif = ' '.join(miflist)
        run.command('mrcat -axis 3 ' + DWImif + ' dwibc.mif')
    run.command('mrconvert -force dwibc.mif working.mif', show=False)

    # End timer
    end_time = time.time()
    elapsed_time = end_time - start_time

    # Convert elapsed time to hours and minutes
    hours, rem = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(rem, 60)

    print(f"Eddy current, EPI, motion correction completed in {int(hours):02} hours, {int(minutes):02} minutes, {int(seconds):02} seconds.")
    print(separator + "\n")

def create_brainmask(fsl_suffix):
    import os, gzip, shutil
    from mrtrix3 import run

    print("...Computing brain mask...")
    run.command('dwiextract -bzero working.mif - | mrmath -axis 3 - mean b0bc.nii')
    # run.command('dwi2mask dwibc.mif - | maskfilter - dilate brain_mask.nii')
    # run.command('fslmaths b0bc.nii -mas brain_mask.nii brain')
    run.command('bet b0bc.nii brain' + fsl_suffix + ' -m -f 0.25')
    if os.path.isfile('brain_mask.nii.gz'):
        with gzip.open('brain_mask' + fsl_suffix, 'rb') as f_in, open('brain_mask.nii', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        run.function(os.remove,'brain_mask' + fsl_suffix, show=False)

def run_rice_bias_correct():
    from mrtrix3 import run, app

    print("...Rician Bias correction...")
    if app.ARGS.denoise:
        run.command('mrcalc noisemap.mif -finite noisemap.mif 0 -if lowbnoisemap.mif', show=False)
        run.command('mrcalc working.mif 2 -pow lowbnoisemap.mif 2 -pow -sub -abs -sqrt - | mrcalc - -finite - 0 -if dwirc.mif')
    else:
        run.command('dwidenoise -noise lowbnoisemap.mif -estimator Exp2 dwi.mif dwitmp.mif', show=False)
        app.cleanup('dwitmp.mif')
        run.command('mrcalc working.mif 2 -pow lowbnoisemap.mif 2 -pow -sub -abs -sqrt - | mrcalc - -finite - 0 -if dwirc.mif')
        if not app.ARGS.degibbs:
            run.command('mrconvert -force -export_grad_fsl working.bvec working.bval working.mif working.nii', show=False)
        run.command('mrconvert -force -fslgrad working.bvec working.bval dwirc.mif working.mif', show=False)

def run_normalization(dwi_metadata):
    from mrtrix3 import app, run

    DWInlist = dwi_metadata['dwi_list']
    idxlist = dwi_metadata['idxlist']

    print("...B0 Normalization...")
    if app.ARGS.normalize and not app.ARGS.b1correct:
        if not len(DWInlist) == 1:
            miflist = []
            for idx,i in enumerate(DWInlist):
                run.command('mrconvert -coord 3 ' + idxlist[idx] + ' working.mif dwiprenm' + str(idx) + '.mif')
                run.command('dwiextract dwiprenm' + str(idx) + '.mif - -bzero | mrmath - mean mean_bzero_prenm' + str(idx) + '.mif -axis 3')
                run.command('mrfilter -stdev 3 mean_bzero_prenm' + str(idx) + '.mif smooth mean_bzero_sm' + str(idx) + '.mif')
                miflist.append('mean_bzero_sm' + str(idx) + '.mif')
            nmlist = []
            nmlist.append('dwiprenm0.mif')
            for idx,i in enumerate(miflist[1:]):
                run.command('mrcalc ' + miflist[0] + ' ' + i + ' -div ratio0to' + str(idx+1) + '.mif')
                run.command('mrcalc dwiprenm' + str(idx+1) + '.mif ratio0to' + str(idx+1) + '.mif -mult dwinm' + str(idx+1) + '.mif')
                nmlist.append('dwinm' + str(idx+1) + '.mif')
                DWImif = ' '.join(nmlist)
            run.command('mrcat -axis 3 ' + DWImif + ' dwinm.mif')
            run.command('mrconvert -force dwinm.mif working.mif')
    elif app.ARGS.normalize and app.ARGS.b1correct:
        if not len(DWInlist) == 1:
            for idx,i in enumerate(DWInlist):
                run.command('mrconvert -force -coord 3 ' + idxlist[idx] + ' working.mif dwiprenm' + str(idx) + '.mif')
            nmlist = []
            nmlist.append('dwiprenm0.mif')
            for idx,i in enumerate(DWInlist[1:]):
                run.command('mrcalc -force biasfield0.mif biasfield' + str(idx+1) + '.mif' + ' -div ratio0to' + str(idx+1) + '.mif')
                run.command('mrcalc -force dwiprenm' + str(idx+1) + '.mif ratio0to' + str(idx+1) + '.mif -mult dwinm' + str(idx+1) + '.mif')
                nmlist.append('dwinm' + str(idx+1) + '.mif')
                DWImif = ' '.join(nmlist)
            run.command('mrcat -axis 3 ' + DWImif + ' dwinm.mif')
            run.command('mrconvert -force dwinm.mif working.mif')
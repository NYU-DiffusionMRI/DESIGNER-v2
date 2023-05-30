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

def run_mppca(args_extent, args_phase, args_shrinkage, args_algorithm):
    """
    wrapper for complex or magnitude, adatpive or local mppca
    
    """
    from mrtrix3 import run, app
    from ants import image_read, image_write, from_numpy
    import numpy as np

    if app.ARGS.adaptive_patch:
        import lib.mpdenoise_sj as mp
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
    print('...denoising...')
    if app.ARGS.adaptive_patch:
        Signal, Sigma, Nparameters = mp.denoise(
            dwi, patchtype='nonlocal', shrinkage=args_shrinkage, algorithm=args_algorithm)
        Sigma[np.isnan(Sigma)] = 0
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

def run_patch2self():
    """
    wrapper for patch2self
    """
    from mrtrix3 import run, app
    from dipy.denoise.patch2self import patch2self
    from ants import image_read, image_write, from_numpy
    import numpy as np

    run.command('mrconvert -force -export_grad_fsl working.bvec working.bval working.mif working.nii', show=False)
    nii = image_read('working.nii')
    dwi = nii.numpy()

    bvals = np.loadtxt('working.bval')
    dwi_dn = patch2self(dwi, bvals)
    out = from_numpy(
        dwi_dn, origin=nii.origin, spacing=nii.spacing, direction=nii.direction)
    image_write(out, 'working_p2s.nii')
    run.command('mrconvert -force -fslgrad working.bvec working.bval working_p2s.nii working.mif', show=False)


def run_degibbs(pf, pe_dir):
    """
    wrapper for rpg degibbs
    """

    import lib.gibbs_removal_rpg as rpg
    from mrtrix3 import run, app
    from ants import image_read, image_write, from_numpy
    import os

    rpg_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'rpg_cpp')

    # convert working.mif to nii
    run.command('mrconvert -force -export_grad_fsl working.bvec working.bval working.mif working.nii', show=False)
    nii = image_read('working.nii')
    dwi = nii.numpy()
    print("...RPG degibbsing...")

    if pe_dir[0] == 'i':
        pe_rpg = 0
    elif pe_dir[0] == 'j':
        pe_rpg = 1
    elif pe_dir[0] == 'k':
        pe_rpg = 2
    
    n_cores = app.ARGS.n_cores

    print('Gibbs correction with parameters:')
    print('phase-encoding direction = %s' % (pe_dir))
    print('partial-fourier factor   = %s' % (pf))


    if 'i' in pe_dir:
        pe_dir = 'x'
    elif 'j' in pe_dir:
        pe_dir = 'y'
    elif 'k' in pe_dir:
        pe_dir = 'z'
    else:
        pe_dir = 'y'

    run.command('%s -pf %s -dim %s %s %s' %
                (os.path.join(rpg_path,'rpg'),
                pf, 
                pe_dir,
                'working.nii',
                'working_rpg.nii'
                ), show=False)

    # img_dg = rpg.gibbs_removal(dwi, pf_fact=pf, pe_dim=pe_rpg, nproc=n_cores)
    # out = from_numpy(
    #     img_dg, origin=nii.origin, spacing=nii.spacing, direction=nii.direction)
    # image_write(out, 'working_rpg.nii')

    #convert gibbs corrected nii to .mif
    run.command('mrconvert -force -fslgrad working.bvec working.bval working_rpg.nii working.mif', show=False)


def group_alignment(group_list):
    from mrtrix3 import app, run, path, image, fsl, MRtrixError
    import numpy as np

    fsl_suffix = fsl.suffix()

    print('... Rigidly aligning groups...')

    # 1) find the longest series and extract the b0
    group_len = []
    for i in range(len(group_list)):
        std = run.command('mrinfo -size %s' %
                          (path.to_scratch(group_list[i])
                           ), show=False)
        group_len.append([int(j) for j in std[0].strip().split(' ')][-1])
    longest_series = np.argmax(group_len)

    # idxlist = dwi_metadata['idxlist']
    # series_inds = [i.split(",") for i in idxlist]
    # series_lens = [len(i) for i in series_inds]
    # longest_series = np.argmax(series_lens)

    # For each group: Register b0 ref PE and b0 RPE to the b0 and run topup + eddy
    for i in range(len(group_list)):
        run.command('dwiextract -force -bzero %s - | mrmath - mean %s -axis 3 -force' %
                    (path.to_scratch(group_list[i]), 
                    path.to_scratch('b0_align_series_' + str(i) + '.nii')),
                    show=False)
        run.command('bet %s %s -f 0.2 -m' %
                    (path.to_scratch('b0_align_series_' + str(i) + '.nii'), 
                    path.to_scratch('b0_align_series_' + str(i) + '_brain')),
                    show=False)

    # 2) extract b0 images from all other series
    for i in range(len(group_list)):                
        if i == longest_series:
            continue
        else:
            # 3) rigidly register all b0 images to the b0 from the longest series
            run.command('flirt -in %s -ref %s -omat %s -dof 6' %
                    (path.to_scratch('b0_align_series_' + str(i) + '_brain'),
                    path.to_scratch('b0_align_series_' + str(longest_series) + '_brain'),
                    path.to_scratch('b0_align_series_to_longest_series_' + str(i) + '.mat')),
                    show=False)
            run.command('transformconvert -force %s %s %s flirt_import %s' %
                (path.to_scratch('b0_align_series_to_longest_series_' + str(i) + '.mat'),
                path.to_scratch('b0_align_series_' + str(i) + '_brain' + fsl_suffix),
                path.to_scratch('b0_align_series_' + str(longest_series) + '_brain' + fsl_suffix),
                path.to_scratch('b0_align_series_scan' + str(i) + 'to_longest_series_mrtrix.txt')),
                show=False)
            run.command('mrtransform -force -linear %s -interp cubic %s %s' %
                (path.to_scratch('b0_align_series_scan' + str(i) + 'to_longest_series_mrtrix.txt'),
                path.to_scratch(group_list[i]),
                path.to_scratch('dwi_align_series_' + str(i) + '_to_ref.mif')),
                show=False)
    
    # concatenate data in correct order
    series_cat_all = ''
    for i in range(len(group_list)):
        if i == longest_series:
            series_to_cat = path.to_scratch(group_list[i])
        else:
            series_to_cat = path.to_scratch('dwi_align_series_' + str(i) + '_to_ref.mif')
        series_cat_all = series_cat_all + series_to_cat + ' '
    run.command('mrcat -force %s %s' %
                (series_cat_all,
                path.to_scratch('dwi_align_series_aligned.mif')),
                show=False)
    run.command('mrconvert -force dwi_align_series_aligned.mif working.mif', show=False)


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

def run_pre_align(dwi_metadata):
    from mrtrix3 import app, run, path, image, fsl, MRtrixError

    if len(dwi_metadata['dwi_list']) == 1:
        print('Warning -pre_align option was used but only 1 input series was provided')
    else:
        series_list = []
        for i in range(len(dwi_metadata['dwi_list'])):
            run.command('mrconvert -coord 3 %s working.mif %s' % 
                    (dwi_metadata['idxlist'][i], path.to_scratch('dwi_align_series_' + str(i) + '.mif')),
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

    print("...Eddy current, EPI, motion correction...")

    eddyopts = '" --cnr_maps --repol --data_is_shelled "'

    fsl_suffix = fsl.suffix()

    pe_dir = dwi_metadata['pe_dir']

    echo_times_per_series = dwi_metadata['TE']
    if  echo_times_per_series.count(echo_times_per_series[0]) == len(echo_times_per_series):
        flag_variable_TE = False
    else:
        flag_variable_TE = True

    if flag_variable_TE:
        if app.ARGS.eddy_groups is None:
            raise MRtrixError('for variable TE data, the user must supply the \
                            -eddy_groups command line argument')
        
         # get TE of the PA
        if app.ARGS.rpe_pair:
            if app.ARGS.rpe_te:
                te_pa = float(app.ARGS.rpe_te)
            else:
                pa_fpath = splitext_(path.from_user(app.ARGS.rpe_pair))[0]
                pa_bids_path = pa_fpath + '.json'
                if not os.path.exists(pa_bids_path):
                    raise MRtrixError('for variable TE data the RPE \
                        image must be accompanies by a bids .json or use\
                        must use the -rpe_te argument to specify the PA echo time')

                pa_bids = json.load(open(pa_bids_path))
                te_pa = float(pa_bids['echo_time'])


            id_dwi_match_pa = np.where(te_pa == np.array(echo_times_per_series))[0][0]
            if id_dwi_match_pa.size == 0:
                raise MRtrixError('echo time of reverse phase encoding image does not\
                                  match any of the input echo times, please check.')

            # extract b0s of dwi matching te of PA image
            run.command('mrconvert -coord 3 %s working.mif - | dwiextract -bzero - - | mrmath - mean %s -axis 3' %
                (dwi_metadata['idxlist'][id_dwi_match_pa],
                path.to_scratch('pe_original_meanb0.nii')))
            # extract brain from mean b0
            run.command('bet %s %s -f 0.2 -m' %
                (path.to_scratch('pe_original_meanb0.nii'), 
                path.to_scratch('pe_original_brain')))

            rpe_size = [ int(s) for s in image.Header(path.from_user(app.ARGS.rpe_pair)).size() ]
            if len(rpe_size) == 4:
                run.command('mrconvert -coord 3 0 %s %s' % 
                    (path.from_user(app.ARGS.rpe_pair), path.to_scratch('rpe_b0.nii')))
            else: 
                run.command('mrconvert %s %s' % 
                    (path.from_user(app.ARGS.rpe_pair), path.to_scratch('rpe_b0.nii')))
            run.command('bet %s %s -f 0.20' % 
                (path.to_scratch('rpe_b0.nii'), path.to_scratch('rpe_b0_brain')))


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
           
            run.command('mrconvert -coord 3 %s working.mif %s' % 
                (volume_idx,
                path.to_scratch('dwi_pre_eddy_' + str(i) + '.mif')),
                show=False)
            
            if app.ARGS.rpe_pair:
                run.command('dwiextract -bzero %s - | mrmath - mean %s -axis 3' %
                    (path.to_scratch('dwi_pre_eddy_' + str(i) + '.mif'), 
                    path.to_scratch('b0_pre_eddy_' + str(i) + '.nii')))
                run.command('bet %s %s -f 0.2 -m' %
                    (path.to_scratch('b0_pre_eddy_' + str(i) + '.nii'), 
                    path.to_scratch('b0_pre_eddy_' + str(i) + '_brain')))

                run.command('flirt -in %s -ref %s -out %s -dof 6' %
                    (path.to_scratch('rpe_b0_brain'),
                    path.to_scratch('b0_pre_eddy_' + str(i) + '_brain'),
                    path.to_scratch('rpe_to_ref_' + str(i))))

                run.command('flirt -in %s -ref %s -out %s -dof 6' %
                    (path.to_scratch('pe_original_brain' + fsl_suffix), 
                    path.to_scratch('b0_pre_eddy_' + str(i) + '_brain'), 
                    path.to_scratch('pe_to_ref_' + str(i))))

                run.command('mrcat -force -axis 3 %s %s %s' %
                    (path.to_scratch('pe_to_ref_' + str(i) + fsl_suffix),
                    path.to_scratch('rpe_to_ref_' + str(i) + fsl_suffix),
                    path.to_scratch('b0_pair_topup_' + str(i) + '.nii')))

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
                np.savetxt(path.to_scratch('topup_acqp.txt'), acqp, fmt="%1.2f")

                run.command('topup --imain=%s --datain=%s --config=b02b0.cnf --scale=1 --out=%s --iout=%s' %
                    (path.to_scratch('b0_pair_topup_' + str(i) + '.nii'),
                    path.to_scratch('topup_acqp.txt'),
                    path.to_scratch('topup_results_' + str(i)),
                    path.to_scratch('topup_results_' + str(i) + fsl_suffix)))
                
                # mask the topup corrected image
                run.command('mrmath %s mean %s -axis 3' %
                            (path.to_scratch('topup_results_' + str(i) + fsl_suffix),
                             path.to_scratch('topup_corrected_' + str(i) + '_mean.nii')
                            ))
                 
                run.command('bet %s %s -f 0.2 -m' %
                    (path.to_scratch('topup_corrected_' + str(i) + '_mean.nii'), 
                    path.to_scratch('topup_corrected_' + str(i) + '_brain')))
                
                run.command('dwifslpreproc -nocleanup -scratch %s -eddy_options %s -rpe_none -eddy_mask %s -topup_files %s -pe_dir %s %s %s' % 
                    (path.to_scratch('eddy_processing_' + str(i)), 
                    eddyopts, 
                    path.to_scratch('topup_corrected_' + str(i) + '_brain_mask' + fsl_suffix),
                    path.to_scratch('topup_results_' + str(i)),
                    pe_dir,
                    path.to_scratch('dwi_pre_eddy_' + str(i) + '.mif'),
                    path.to_scratch('dwi_post_eddy_' + str(i) + '.mif')))
                
            elif app.ARGS.rpe_none:

                run.command('dwifslpreproc -nocleanup -scratch %s -eddy_options %s -pe_dir %s %s %s' % 
                    (path.to_scratch('eddy_processing_' + str(i)), 
                    eddyopts, 
                    pe_dir,
                    path.to_scratch('dwi_pre_eddy_' + str(i) + '.mif'),
                    path.to_scratch('dwi_post_eddy_' + str(i) + '.mif')))
                
            elif app.ARGS.rpe_all:
                ### this needs to be changed for realistic compatiblity to eddy_groups
                run.command('mrconvert -export_grad_mrtrix grad.txt working.mif tmp.mif', show=False)
                run.command('mrconvert -grad grad.txt ' + path.from_user(app.ARGS.rpe_all) + ' dwirpe.mif', show=False)
                run.command('mrcat -axis 3 working.mif dwirpe.mif dwipe_rpe.mif')
                run.command('dwifslpreproc -nocleanup -eddy_options %s -rpe_all -pe_dir %s %s %s' %
                    (eddyopts, 
                     pe_dir, 
                     path.to_scratch('dwi_pre_eddy_' + str(i) + '.mif'),
                     path.to_scratch('dwi_post_eddy_' + str(i) + '.mif')))
                run.function(os.remove,'tmp.mif')

            elif app.ARGS.rpe_header:

                run.command('dwifslpreproc -nocleanup -eddy_options %s -rpe_header %s %s' % 
                (eddyopts,
                path.to_scratch('dwi_pre_eddy_' + str(i) + '.mif'),
                path.to_scratch('dwi_post_eddy_' + str(i) + '.mif')))

            run.command('dwiextract -bzero %s - | mrmath - mean %s -axis 3' %
                    (path.to_scratch('dwi_post_eddy_' + str(i) + '.mif'), 
                    path.to_scratch('b0_post_eddy_' + str(i) + '.nii')))

            # undo grouping:
            for j in range(len(group_idx)):
                run.command('mrconvert -coord 3 %s %s %s' %
                            (ginds[j],
                            path.to_scratch('dwi_post_eddy_' + str(i) + '.mif'),
                            path.to_scratch('dwi_post_eddy_ungrouped_' + str(group_idx[j]) + '.mif')),
                            show=False)
                
        
        registered_data = glob.glob(os.path.join(path.to_scratch(''),'dwi_post_eddy_ungrouped*.mif'))
        registered_data.sort()

        group_alignment(registered_data)
        run.command('mrconvert -force working.mif dwiec.mif', show=False)
        

    # if not variable TE (not using eddy_groups - only run eddy once)
    else:

        if app.ARGS.rpe_pair:

            run.command('dwiextract -bzero dwi.mif - | mrconvert -coord 3 0 - b0pe.nii')
            rpe_size = [ int(s) for s in image.Header(path.from_user(app.ARGS.rpe_pair)).size() ]
            if len(rpe_size) == 4:
                run.command('mrconvert -coord 3 0 ' + path.from_user(app.ARGS.rpe_pair) + ' b0rpe.nii')
            else: 
                run.command('mrconvert ' + path.from_user(app.ARGS.rpe_pair) + ' b0rpe.nii')
            run.command('flirt -in b0rpe.nii -ref b0pe.nii -dof 6 -out b0rpe2pe.nii.gz')
            run.command('mrcat -axis 3 b0pe.nii b0rpe2pe.nii.gz b0_pair_topup.nii')

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
            np.savetxt(path.to_scratch('topup_acqp.txt'), acqp, fmt="%1.2f")

            run.command('topup --imain=%s --datain=%s --config=b02b0.cnf --scale=1 --out=%s --iout=%s' %
                (path.to_scratch('b0_pair_topup.nii'),
                path.to_scratch('topup_acqp.txt'),
                path.to_scratch('topup_results'),
                path.to_scratch('topup_results' + fsl_suffix)))
                
            # mask the topup corrected image
            run.command('mrmath %s mean %s -axis 3' %
                        (path.to_scratch('topup_results' + fsl_suffix),
                            path.to_scratch('topup_corrected_mean.nii')
                        ))
                 
            run.command('bet %s %s -f 0.2 -m' %
                (path.to_scratch('topup_corrected_mean.nii'), 
                path.to_scratch('topup_corrected_brain')))
            
            run.command('dwifslpreproc -nocleanup -scratch %s -eddy_options %s -rpe_none -eddy_mask %s -topup_files %s -pe_dir %s %s %s' % 
                    (path.to_scratch('eddy_processing'), 
                    eddyopts, 
                    path.to_scratch('topup_corrected_brain_mask' + fsl_suffix),
                    path.to_scratch('topup_results'),
                    pe_dir,
                    path.to_scratch('working.mif'),
                    path.to_scratch('dwiec.mif')))

        elif app.ARGS.rpe_none:

            run.command('dwiextract -bzero dwi.mif - | mrconvert -coord 3 0 - b0pe.nii')
            run.command('bet %s %s -f 0.2 -m' %
                (path.to_scratch('b0pe.nii'), 
                path.to_scratch('b0_pe_brain')))
            
            run.command('dwifslpreproc -nocleanup -scratch %s -eddy_options %s -rpe_none -eddy_mask %s -pe_dir %s %s %s' % 
                    (path.to_scratch('eddy_processing'), 
                    eddyopts, 
                    path.to_scratch('b0_pe_brain_mask' + fsl_suffix),
                    pe_dir,
                    path.to_scratch('working.mif'),
                    path.to_scratch('dwiec.mif')))
            
        elif app.ARGS.rpe_all:

            run.command('mrconvert -export_grad_mrtrix grad.txt dwi.mif tmp.mif', show=False)
            run.command('mrconvert -grad grad.txt ' + path.from_user(app.ARGS.rpe_all) + ' dwirpe.mif')
            run.command('mrcat -axis 3 working.mif dwirpe.mif dwipe_rpe.mif')
            cmd = ('dwifslpreproc -nocleanup -eddy_options %s -rpe_all -pe_dir %s dwipe_rpe.mif dwiec.mif' %
                (eddyopts, pe_dir))
            run.command('cmd')
            run.function(os.remove,'tmp.mif')

        elif app.ARGS.rpe_header:

            cmd = ('dwifslpreproc -nocleanup -eddy_options %s -rpe_header working.mif dwiec.mif' % 
                (eddyopts))
            run.command(cmd)

        elif not app.ARGS.rpe_header and not app.ARGS.rpe_all and not app.ARGS.rpe_pair:
            raise MRtrixError("the eddy option must run alongside -rpe_header, -rpe_all, or -rpe_pair option")

    run.command('mrconvert -force dwiec.mif working.mif', show=False)

def run_b1correct(dwi_metadata):
    from mrtrix3 import run

    DWInlist = dwi_metadata['dwi_list']
    idxlist = dwi_metadata['idxlist']

    # b1 bias field correction
    print("...B1 correction...")
    if len(DWInlist) == 1:
        run.command('dwibiascorrect ants -bias biasfield.mif working.mif dwibc.mif')
    else:
        # b1 correction may still need to be done individually for each diffusion series ...
        miflist = []
        for idx,i in enumerate(DWInlist):
            cmd = ('mrconvert -coord 3 %s working.mif dwiprebc%s.mif' % 
                (idxlist[idx], str(idx)))
            run.command(cmd)
            cmd = ('dwibiascorrect ants -bias biasfield%s.mif dwiprebc%s.mif dwibc%s.mif' % 
                (str(idx), str(idx), str(idx)))
            run.command(cmd)
            miflist.append('dwibc' + str(idx) + '.mif')
            DWImif = ' '.join(miflist)
        run.command('mrcat -axis 3 ' + DWImif + ' dwibc.mif')
    run.command('mrconvert -force dwibc.mif working.mif', show=False)

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
    run.command('mrconvert -force dwirc.mif working.mif')

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
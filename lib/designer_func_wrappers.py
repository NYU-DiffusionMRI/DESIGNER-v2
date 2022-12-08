"""
Utilities for checking compatibility of inputs to designer

Includes:
run_denoising:

run_degibbs:

run_eddy: 

run_b1correct:

create_brainmask:

run_rice_bias_correct:

run_normalization:
"""

def run_denoising(args_extent, args_phase, args_shrinkage, args_algorithm):
    import mpcomplex as mp #pylint: disable=import-error
    from mrtrix3 import run, app
    from ants import image_read, image_write, from_numpy
    import numpy as np

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
    #sx,sy,sz,N = np.shape(dwi)

    print('...denoising...')
    # nonlocal denoising, no shrinkage
    Signal, Sigma, Nparameters = mp.denoise(
        dwi, phase=args_phase, kernel=extent, shrinkage=args_shrinkage, algorithm=args_algorithm
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

def run_degibbs(pf, pe_dir):
    import gibbs_removal_rpg as rpg
    from mrtrix3 import run
    from ants import image_read, image_write, from_numpy

    # convert working.mif to nii
    run.command('mrconvert -export_grad_fsl working.bvec working.bval working.mif working.nii', show=False)
    nii = image_read('working.nii')
    dwi = nii.numpy()
    print("...RPG degibbsing...")

    if pe_dir[0] == 'i':
        pe_rpg = 0
    elif pe_dir[0] == 'j':
        pe_rpg = 1
    elif pe_dir[0] == 'k':
        pe_rpg = 2

    print('Gibbs correction with parameters:')
    print('phase-encoding direction = %s' % (pe_dir))
    print('partial-fourier factor   = %s' % (pf))
    
    img_dg = rpg.gibbs_removal(dwi, pf_fact=pf, pe_dim=pe_rpg, nproc=8)
    out = from_numpy(
        img_dg, origin=nii.origin, spacing=nii.spacing, direction=nii.direction)
    image_write(out, 'working_rpg.nii')

    #convert gibbs corrected nii to .mif
    run.command('mrconvert -force -fslgrad working.bvec working.bval working_rpg.nii working.mif', show=False)

def run_eddy(pe_dir, rpe):
    from mrtrix3 import app, run, path, image
     # epi + eddy current and motion correction
    # if number of input volumes is greater than 1, make a new acqp and index file.
    # eddyindexlist = [(i+1)*np.ones((1,nvols[(i+1)])) for i in range(len(DWInlist))]
    # eddyindex = np.hstack(eddyindexlist)
    # np.savetxt('eddyindex.txt',eddyindex,fmt="%d")
    # idxpath = path.to_scratch('eddyindex.txt')

    # acqp = np.zeros((len(DWInlist),4))
    # for i in range(len(DWInlist)):
    #     if app.ARGS.pe_dir == 'AP' or app.ARGS.pe_dir == '-j':
    #         acqp[i,:] = (0,-1,0,0.1)
    #     if app.ARGS.pe_dir == 'PA' or app.ARGS.pe_dir == 'j':
    #         acqp[i,:] = (0,1,0,0.1)
    #     if app.ARGS.pe_dir == 'LR' or app.ARGS.pe_dir == 'i':
    #         acqp[i,:] = (1,0,0,0.1)
    #     if app.ARGS.pe_dir == 'RL' or app.ARGS.pe_dir == '-i':
    #         acqp[i,:] = (-1,0,0,0.1)
    #     if app.ARGS.pe_dir == 'IS' or app.ARGS.pe_dir == 'k':
    #         acqp[i,:] = (0,0,1,0.1)
    #     if app.ARGS.pe_dir == 'SI' or app.ARGS.pe_dir == '-k':
    #         acqp[i,:] = (0,0,-1,0.1)
    # np.savetxt('eddyacqp.txt',acqp,fmt="%1.1f")
    # acqpath = path.to_scratch('eddyacqp.txt')

    print("...Eddy current, EPI, motion correction...")
    eddyopts = '" --cnr_maps --repol --data_is_shelled "'
    #eddyopts = '" --cnr_maps --repol --data_is_shelled --ol_type=both --mb=' + app.ARGS.mb + ' "'
    #print('dwifslpreproc -nocleanup -scratch ' + path.to_scratch('',True) + '/eddy_processing' + ' -eddy_options ' + eddyopts + ' -rpe_none -pe_dir ' + app.ARGS.pe_dir + ' working.mif dwiec.mif')
    if app.ARGS.rpe_none:
        #run.command('dwifslpreproc -nocleanup -eddy_options " --index=' + idxpath + ' --repol --data_is_shelled" -rpe_none -pe_dir ' + app.ARGS.pe_dir + ' working.mif dwiec.mif')
        cmd = ('dwifslpreproc -nocleanup -scratch %s/eddy_processing -eddy_options %s -rpe_none -pe_dir %s working.mif dwiex.mif' % 
            (path.to_scratch('',True), eddyopts, pe_dir))
        run.command(cmd)
        #run.command('dwifslpreproc -nocleanup -scratch ' + path.to_scratch('',True) + '/eddy_processing' + ' -eddy_options ' + eddyopts + ' -rpe_none -pe_dir ' + app.ARGS.pe_dir + ' working.mif dwiec.mif')
    elif app.ARGS.rpe_pair:
        run.command('dwiextract -bzero dwi.mif - | mrconvert -coord 3 0 - b0pe.nii')
        rpe_size = [ int(s) for s in image.Header(path.from_user(rpe)).size() ]
        if len(rpe_size) == 4:
            run.command('mrconvert -coord 3 0 ' + path.from_user(rpe) + ' b0rpe.nii')
        else: 
            run.command('mrconvert ' + path.from_user(rpe) + ' b0rpe.nii')
        run.command('flirt -in b0rpe.nii -ref b0pe.nii -dof 6 -out b0rpe2pe.nii.gz')
        run.command('mrcat -axis 3 b0pe.nii b0rpe2pe.nii.gz rpepair.mif')
        #run.command('dwifslpreproc -nocleanup -eddy_options " --acqp=' + acqpath + ' --index=' + idxpath + ' --repol --data_is_shelled" -rpe_pair -se_epi rpepair.mif -align_seepi -pe_dir ' + app.ARGS.pe_dir + ' working.mif dwiec.mif')
        cmd = ('dwifslpreproc -nocleanup -scratch %s/eddy_processing -eddy_options %s -rpe_pair -se_epi rpepair.mif -align_seepi -pe_dir %s working.mif dwiec.mif' % 
            (path.to_scratch('',True), eddyopts, pe_dir))
        run.command(cmd)
    elif app.ARGS.rpe_all:
        run.command('mrconvert -export_grad_mrtrix grad.txt dwi.mif tmp.mif', show=False)
        run.command('mrconvert -grad grad.txt ' + path.from_user(rpe) + ' dwirpe.mif')
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
        print("the eddy option must run alongside -rpe_header, -rpe_all, or -rpe_pair option")
        quit()
    run.command('mrconvert -force dwiec.mif working.mif')

def run_b1correct(DWInlist, idxlist):
    from mrtrix3 import run

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
        run.function(os.remove,'brain_mask' + fsl_suffix)

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

def run_normalization(DWInlist, idxlist):
    from mrtrix3 import app, run

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
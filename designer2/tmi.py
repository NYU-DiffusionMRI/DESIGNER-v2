#!/usr/bin/env python3
import os

from lib.designer_input_utils import get_input_info, convert_input_data, create_shell_table, assert_inputs
from lib.designer_fit_wrappers import refit_or_smooth, save_params

def usage(cmdline): #pylint: disable=unused-variable
    from mrtrix3 import app #pylint: disable=no-name-in-module, import-outside-toplevel

    cmdline.set_copyright("""Copyright (c) 2016 New York University.\n\n
                    
                        
                    Permission is hereby granted, free of charge, to any non-commercial entity (\'Recipient\') obtaining a copy of this software and associated documentation files (the \'Software\'), to the Software solely for non-commercial research, including the rights to use, copy and modify the Software, subject to the following conditions\n
                    
                    1. The above copyright notice and this permission notice shall be included by Recipient in all copies or substantial portions of the Software.\n

                    2. THE SOFTWARE IS PROVIDED \'AS IS\', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIESvOF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF ORIN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.\n
                    
                    3. In no event shall NYU be liable for direct, indirect, special, incidental or consequential damages in connection with the Software. Recipient will defend, indemnify and hold NYU harmless from any claims or liability resulting from the use of the Software by recipient.\n
                    
                    4. Neither anything contained herein nor the delivery of the Software to recipient shall be deemed to grant the Recipient any right or licenses under any patents or patent application owned by NYU.\n
                                            
                    5. The Software may only be used for non-commercial research and may not be used for clinical care.\n
                    
                    6. Any publication by Recipient of research involving the Software shall cite the references listed below.\n
                                                                
                    """)


    cmdline.set_author('Benjamin Ades-Aron (benjamin.ades-aron@nyulangone.org)')
    cmdline.set_synopsis("""1. pre-check: concatenate all dwi series and make sure the all diffusion AP images and PA image have the same matrix dimensions/orientations \n 
 						
                        2. Run fitting: Available options include DTI (TMI will look for b0 and b1000 images by default.\n
                            DKI (constrained fitting is available using the -fit_constraints option).\n
                            W-DKI (DKI parameters using the W definition instead of K).\n
                            WMTI (using dipy implimentation from Rafael Henriques).\n
                            SMI (From Santiago Coelho).\n
                        
                        """)

    cmdline.add_citation('Veraart, J.; Novikov, D.S.; Christiaens, D.; Ades-aron, B.; Sijbers, J. & Fieremans, E. Denoising of diffusion MRI using random matrix theory. NeuroImage, 2016, 142, 394-406, doi: 10.1016/j.neuroimage.2016.08.016',is_external=True)

    cmdline.add_argument('input',  help='The input DWI series. For multiple input series, separate file names with commas (i.e. dwi1.nii,dwi2.nii,...)')
    cmdline.add_argument('output', help='The output directory (includes diffusion parameters, kurtosis parameters and processed dwi) unless option -processing_only is used, in which case this is the output basename')
    
    options = cmdline.add_argument_group('Other options for the TMI script')    
    options.add_argument('-mask', metavar=('<path>'),help='path to mask previously computed by designer script')
    options.add_argument('-datatype', metavar=('<spec>'), help='Valid options are float32, float32le, float32be, float64, float64le, float64be, int64, uint64, int64le, uint64le, int64be, uint64be, int32, uint32, int32le, uint32le, int32be, uint32be, int16, uint16, int16le, uint16le, int16be, uint16be, cfloat32, cfloat32le, cfloat32be, cfloat64, cfloat64le, cfloat64be, int8, uint8, bit')
    options.add_argument('-fslbvec',metavar=('<bvecs>'),help='specify bvec path if path is different from the path to the dwi or the file has an unusual extention')
    options.add_argument('-fslbval',metavar=('<bvals>'),help='specify bval path if path is different from the path to the dwi or the file has an unusual extention')
    options.add_argument('-bids',metavar=('<bids>'),help='specify bids.json path if path is different from the path to the dwi or the file has an unusual extention')
    options.add_argument('-n_cores',metavar=('<ncores>'),help='specify the number of cores to use in parallel tasts, by default designer will use available cores - 2', default=-3)
    options.add_argument('-echo_time',metavar=('<TE1,TE2,...>'),help='specify the echo time used in the acquisition (comma separated list the same length as number of inputs)')
    options.add_argument('-bshape',metavar=('<beta1,beta2,...>'),help='specify the b-shape used in the acquisition (comma separated list the same length as number of inputs)')

    dki_options = cmdline.add_argument_group('tensor options for the TMI script')
    dki_options.add_argument('-DKI', action='store_true', help='Include DKI parameters in output folder (mk,ak,rk)')
    dki_options.add_argument('-DTI', action='store_true', help='Include DTI parameters in output folder (md,ad,rd,fa,eigenvalues, eigenvectors')
    dki_options.add_argument('-WDKI', action='store_true',help='output the kurtosis tensor with W cumulant rather than K')    
    dki_options.add_argument('-akc_outliers', action='store_true', help='brute force K tensor outlier rejection')
    dki_options.add_argument('-fit_constraints',metavar=('<string>'),help='constrain the wlls fit (default 0,1,0)')
    dki_options.add_argument('-fit_smoothing',metavar=('<percentile>'),help='NLM smoothing on wlls fit')
    dki_options.add_argument('-polyreg',action='store_true',help='polynomial regression based DKI estimation')

    smi_options = cmdline.add_argument_group('tensor options for the TMI script')
    smi_options.add_argument('-SMI', action='store_true',help='Perform estimation of SMI (standard model of Diffusion in White Matter). Please use in conjunction with the -bshape, -echo_time, -sigma, and -compartments options.')  
    smi_options.add_argument('-compartments', metavar=('<compartments>'),help='SMI compartments (IAS, EAS, and FW), default=IAS,EAS')  
    smi_options.add_argument('-sigma', metavar=('<noisemap>'),help='path to noise map for SMI parameter estimation')

    #wmti_options = cmdline.add_argument_group('tensor options for the TMI script')
    #wmti_options.add_argument('-WMTI', action='store_true', help='Include WMTI parameters in output folder (awf,ias_params,eas_params)')


def execute(): #pylint: disable=unused-variable
    from mrtrix3 import app, path, run, MRtrixError #pylint: disable=no-name-in-module, import-outside-toplevel
    import lib.tensor as tensor
    import numpy as np
    from ants import image_read
    import pandas as pd

    app.make_scratch_dir()

    dwi_metadata = get_input_info(
        app.ARGS.input, app.ARGS.fslbval, app.ARGS.fslbvec, app.ARGS.bids)

    assert_inputs(dwi_metadata, None, None)

    convert_input_data(dwi_metadata)

    shell_table = create_shell_table(dwi_metadata)
    shell_rows = ['b-value', 'b-shape', 'n volumes', 'echo time']
    shell_df = pd.DataFrame(data = shell_table, 
                  index = shell_rows)
    print('input DWI data has properties:')
    print(shell_df)

    app.goto_scratch_dir()

    run.command('mrconvert dwi.mif -export_grad_fsl dwi.bvec dwi.bval dwi.nii', show=False)
    nii = image_read('dwi.nii')
    dwi = nii.numpy()

    outdir = path.from_user(app.ARGS.output, True)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    nii = image_read('dwi.nii')
    dwi = nii.numpy()
    bvec = np.loadtxt('dwi.bvec')
    bval = np.loadtxt('dwi.bval')

    order = np.floor(np.log(abs(np.max(bval)+1)) / np.log(10))
    if order >= 2:
        bval = bval / 1000

    if app.ARGS.mask:
        mask = image_read(path.from_user(app.ARGS.mask)).numpy()
    else:
        mask = np.ones(dwi.shape[:-1])
    mask = mask.astype(bool)

    if app.ARGS.fit_constraints:
        constraints = app.ARGS.fit_constraints
        if type(constraints) == str:
                constraints = constraints.split(",")
                constraints = [int(i) for i in constraints]
        else:
            raise MRtrixError("Constraints must be a 3 element comma separated string (i.e. 0,1,0)")
    else:
        constraints = [0,0,0]

    if (len(set(dwi_metadata['bshape_per_volume'])) > 1) or (len(set(dwi_metadata['echo_time_per_volume'])) > 1):
        multi_te_beta = True
        dwi_orig = dwi.copy()
        bval_orig = bval.copy()
        bvec_orig = bvec.copy()
    else:
        multi_te_beta = False


    if len(set(dwi_metadata['bshape_per_volume'])) > 1:
        if (app.ARGS.DTI) or (app.ARGS.DKI) or (app.ARGS.WDKI):
            print('For variable b-shape data DTI/DKI are run only on LTE part')
            lte_idx = (dwi_metadata['bshape_per_volume'] == 1)

            if np.sum(lte_idx) < 6:
                raise MRtrixError("Not enough LTE measurements for DTI/DKI")
            dwi = dwi_orig[:,:,:,lte_idx]
            bval = bval_orig[lte_idx]
            bvec = bvec_orig[:,lte_idx]


    if len(set(dwi_metadata['echo_time_per_volume'])) == 1:

        if app.ARGS.DTI:
            print('...Single shell DTI fit...')
            dtishell = (bval <= 0.1) | ((bval > .5) & (bval <= 1.5))
            dwi_dti = dwi[:,:,:,dtishell]
            bval_dti = bval[dtishell]
            bvec_dti = bvec[:,dtishell]
            grad_dti = np.hstack((bvec_dti.T, bval_dti[None,...].T))
            dti = tensor.TensorFitting(grad_dti, int(app.ARGS.n_cores))
            dt_dti, s0_dti, b_dti = dti.dti_fit(dwi_dti, mask)

        if app.ARGS.DKI or app.ARGS.WDKI:
            print('...Multi shell DKI fit with constraints = ' + str(constraints))

            maxb = 2.5
            dwi_dki = dwi[:,:,:,bval < maxb]
            bvec_dki = bvec[:,bval < maxb]
            bval_dki = bval[bval < maxb]

            grad = np.hstack((bvec_dki.T, bval_dki[None,...].T))
            dki = tensor.TensorFitting(grad, int(app.ARGS.n_cores))
            dt_dki, s0_dki, b_dki = dki.dki_fit(dwi_dki, mask, constraints=constraints)

        if app.ARGS.polyreg and app.ARGS.DKI:
            dt_poly = dki.train_rotated_bayes_fit(dwi_dki, dt_dki, s0_dki, b_dki, mask)

        if app.ARGS.akc_outliers:
            from lib.mpunits import vectorize
            import scipy.io as sio

            dwd = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            mat = sio.loadmat(os.path.join(dwd,'constant','dirs10000.mat'))
            dir = mat['dir']

            print('...Outlier detection...')
            if not (app.ARGS.DKI or app.ARGS.WDKI):
                raise MRtrixError("AKC Outlier detection must be accompanied by DKI option")
            else:
                akc_mask = dki.outlierdetection(dt_dki, mask, dir)
                
            akc_mask = vectorize(akc_mask, mask).astype(bool)
            print('N outliers = %s' % (np.sum(akc_mask)))

            dwi_new = refit_or_smooth(akc_mask, dwi_dki, n_cores=int(app.ARGS.n_cores))
            dt_new,_,_ = dki.dki_fit(dwi_new, akc_mask)

            x,y,z = np.where(akc_mask == 1)
            DT = vectorize(dt_dki, mask)
            DT[x,y,z,:] = dt_new.T
            dt_dki = vectorize(DT, mask)
            akc_mask = dki.outlierdetection(dt_dki, mask, dir)
            akc_mask = vectorize(akc_mask, mask).astype(bool)
            print('N outliers = %s' % (np.sum(akc_mask)))
        else:
            akc_mask = np.zeros_like(mask)

        if app.ARGS.fit_smoothing:
            print('...Nonlocal smoothing...')
            if app.ARGS.DTI:
                dwi_new = refit_or_smooth(akc_mask, dwi_dti, mask=mask, smoothlevel=int(app.ARGS.fit_smoothing))
                dt_dti,_,_ = dti.dti_fit(dwi_new, mask)
            if (app.ARGS.DKI or app.ARGS.WDKI):
                dwi_new = refit_or_smooth(akc_mask, dwi_dki, mask=mask, smoothlevel=int(app.ARGS.fit_smoothing))
                dt_dki,_,_ = dki.dki_fit(dwi_new, mask)

        if app.ARGS.DTI:
            print('...extracting and saving DTI maps...')
            params_dti = dti.extract_parameters(dt_dti, b_dti, mask, extract_dti=True, extract_dki=False, fit_w=False)
            save_params(params_dti, nii, model='dti', outdir=outdir)

        if app.ARGS.DKI:
            print('...extracting and saving DKI maps...')
            params_dki = dki.extract_parameters(dt_dki, b_dki, mask, extract_dti=True, extract_dki=True, fit_w=False)
            save_params(params_dki, nii, model='dki', outdir=outdir)

        if app.ARGS.polyreg:
            params_dki_poly = dki.extract_parameters(dt_poly, b_dki, mask, extract_dti=True, extract_dki=True, fit_w=False)
            save_params(params_dki_poly, nii, model='dki_poly', outdir=outdir)

        if app.ARGS.WDKI:
            print('...extracting and saving WDKI maps...')
            params_dwi = dki.extract_parameters(dt_dki, b_dki, mask, extract_dti=False, extract_dki=True, fit_w=True)
            save_params(params_dwi, nii, model='wdki', outdir=outdir)

    else:

        for te in set(dwi_metadata['echo_time_per_volume']):

            te_idx = (dwi_metadata['echo_time_per_volume'] == te) & (dwi_metadata['bshape_per_volume'] == 1)
            dwi = dwi_orig[:,:,:,te_idx]
            bval = bval_orig[te_idx]
            bvec = bvec_orig[:,te_idx]

            if app.ARGS.DTI:

                if np.sum(te_idx) < 6:
                    continue

                print('...Single shell DTI fit for TE=' + str(te) + '...')
                dtishell = (bval <= 0.1) | ((bval > .5) & (bval <= 1.5))
                dwi_dti = dwi[:,:,:,dtishell]
                bval_dti = bval[dtishell]
                bvec_dti = bvec[:,dtishell]
                grad_dti = np.hstack((bvec_dti.T, bval_dti[None,...].T))
                dti = tensor.TensorFitting(grad_dti, int(app.ARGS.n_cores))
                dt_dti, s0_dti, b_dti = dti.dti_fit(dwi_dti, mask)

            if app.ARGS.DKI or app.ARGS.WDKI:

                if np.sum(te_idx) < 21:
                    continue

                print('...Multi shell DKI fit  for TE=' + str(te) +  ' with constraints = ' + str(constraints))

                maxb = 2.5
                dwi_dki = dwi[:,:,:,bval < maxb]
                bvec_dki = bvec[:,bval < maxb]
                bval_dki = bval[bval < maxb]

                grad = np.hstack((bvec_dki.T, bval_dki[None,...].T))
                dki = tensor.TensorFitting(grad, int(app.ARGS.n_cores))
                dt_dki, s0_dki, b_dki = dki.dki_fit(dwi_dki, mask, constraints=constraints)

            if app.ARGS.polyreg and app.ARGS.DKI:
                dt_poly = dki.train_rotated_bayes_fit(dwi_dki, dt_dki, s0_dki, b_dki, mask)

            if app.ARGS.akc_outliers:
                from lib.mpunits import vectorize
                import scipy.io as sio

                dwd = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
                mat = sio.loadmat(os.path.join(dwd,'constant','dirs10000.mat'))
                dir = mat['dir']

                print('...Outlier detection...')
                if not (app.ARGS.DKI or app.ARGS.WDKI):
                    raise MRtrixError("AKC Outlier detection must be accompanied by DKI option")
                else:
                    akc_mask = dki.outlierdetection(dt_dki, mask, dir)
                    
                akc_mask = vectorize(akc_mask, mask).astype(bool)
                print('N outliers = %s' % (np.sum(akc_mask)))

                dwi_new = refit_or_smooth(akc_mask, dwi_dki, n_cores=int(app.ARGS.n_cores))
                dt_new,_,_ = dki.dki_fit(dwi_new, akc_mask)

                x,y,z = np.where(akc_mask == 1)
                DT = vectorize(dt_dki, mask)
                DT[x,y,z,:] = dt_new.T
                dt_dki = vectorize(DT, mask)
                akc_mask = dki.outlierdetection(dt_dki, mask, dir)
                akc_mask = vectorize(akc_mask, mask).astype(bool)
                print('N outliers = %s' % (np.sum(akc_mask)))
            else:
                akc_mask = np.zeros_like(mask)

            if app.ARGS.fit_smoothing:
                print('...Nonlocal smoothing...')
                if app.ARGS.DTI:
                    dwi_new = refit_or_smooth(akc_mask, dwi_dti, mask=mask, smoothlevel=int(app.ARGS.fit_smoothing))
                    dt_dti,_,_ = dti.dti_fit(dwi_new, mask)
                if (app.ARGS.DKI or app.ARGS.WDKI):
                    dwi_new = refit_or_smooth(akc_mask, dwi_dki, mask=mask, smoothlevel=int(app.ARGS.fit_smoothing))
                    dt_dki,_,_ = dki.dki_fit(dwi_new, mask)

            if app.ARGS.DTI:
                print('...extracting and saving DTI maps...')
                params_dti = dti.extract_parameters(dt_dti, b_dti, mask, extract_dti=True, extract_dki=False, fit_w=False)
                save_params(params_dti, nii, model='dti_te'+str(te), outdir=outdir)

            if app.ARGS.DKI:
                print('...extracting and saving DKI maps...')
                params_dki = dki.extract_parameters(dt_dki, b_dki, mask, extract_dti=True, extract_dki=True, fit_w=False)
                save_params(params_dki, nii, model='dki_te'+str(te), outdir=outdir)

            if app.ARGS.polyreg:
                params_dki_poly = dki.extract_parameters(dt_poly, b_dki, mask, extract_dti=True, extract_dki=True, fit_w=False)
                save_params(params_dki_poly, nii, model='dki_poly_te'+str(te), outdir=outdir)

            if app.ARGS.WDKI:
                print('...extracting and saving WDKI maps...')
                params_dwi = dki.extract_parameters(dt_dki, b_dki, mask, extract_dti=False, extract_dki=True, fit_w=True)
                save_params(params_dwi, nii, model='wdki_te'+str(te), outdir=outdir)


    # if app.ARGS.WMTI:
    #     import dipy.reconst.dki as dki
    #     import dipy.reconst.dki_micro as dki_micro
    #     from dipy.core.gradients import gradient_table

    #     print('...WMTI fit...')
    #     gtab = gradient_table(bval*1000, bvec)
    #     dki_micro_model = dki_micro.KurtosisMicrostructureModel(gtab)
    #     dki_micro_fit = dki_micro_model.fit(dwi, mask=mask)
        
    #     params_wmti = {}
    #     params_wmti['awf'] = dki_micro_fit.awf
    #     params_wmti['tort'] = dki_micro_fit.tortuosity
    #     params_wmti['had'] = dki_micro_fit.hindered_ad
    #     params_wmti['hrd'] = dki_micro_fit.hindered_rd
    #     params_wmti['axd'] = dki_micro_fit.axonal_diffusivity
    #     save_params(params_wmti, nii, model='wmti', outdir=outdir)

    if app.ARGS.SMI:
        from lib.smi import SMI
        import warnings
        warnings.simplefilter('always', UserWarning) 

        if not app.ARGS.sigma:
            warnings.warn('SMI is poorly conditioned without prior estimate of sigma')
            sigma = None
        else:
            sigma = image_read(path.from_user(app.ARGS.sigma)).numpy()

        if app.ARGS.compartments:
            compartments = app.ARGS.compartments
            if type(compartments) == str:
                compartments = compartments.split(",")
                compartments = [str(i) for i in compartments]
            else:
                raise MRtrixError(" Compartments must be a comma sepearated string (i.e. IAS,EAS)")
        else: 
            compartments = ['IAS', 'EAS']

        print('...SMI fit...')
        if multi_te_beta:
            smi = SMI(bval=bval_orig, bvec=bvec_orig)
            smi.set_compartments(compartments)
            smi.set_echotime(dwi_metadata['echo_time_per_volume'])
            smi.set_bshape(dwi_metadata['bshape_per_volume'])
            params_smi = smi.fit(dwi_orig, mask=mask, sigma=sigma)
            save_params(params_smi, nii, model='smi', outdir=outdir)
        else:
            smi = SMI(bval=bval, bvec=bvec)
            smi.set_compartments(compartments)
            smi.set_echotime(dwi_metadata['echo_time_per_volume'])
            smi.set_bshape(dwi_metadata['bshape_per_volume'])
            params_smi = smi.fit(dwi, mask=mask, sigma=sigma)
            save_params(params_smi, nii, model='smi', outdir=outdir)

def main():
    import mrtrix3
    mrtrix3.execute() #pylint: disable=no-member

if __name__ == "__main__":
    main()


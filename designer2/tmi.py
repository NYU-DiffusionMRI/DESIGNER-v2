#!/usr/bin/env python3
import os
import logging
import json
from logging import StreamHandler, FileHandler

from lib.designer_input_utils import get_input_info, convert_input_data, create_shell_table, assert_inputs
from lib.designer_fit_wrappers import refit_or_smooth, save_params

# List of keys to exclude from logs
EXCLUDED_KEYS = {
    "msg", "levelname", "levelno", "exc_info", "exc_text", 
    "stack_info", "created", "msecs", "relativeCreated", 
    "thread", "threadName", "processName", "process", "args"
}

# Set up logging in JSON format
class JsonFormatter(logging.Formatter):
    def format(self, record):
        log_record = {
            'level': record.levelname,
            'time': self.formatTime(record, self.datefmt),
            'message': record.getMessage(),
            'name': record.name,
            'pathname': record.pathname,
            'lineno': record.lineno,
        }
        log_record.update({
            key: value for key, value in record.__dict__.items()
            if key not in EXCLUDED_KEYS and key not in log_record
        })
        return json.dumps(log_record)
    
# Simple Formatter for StreamHandler (console output)
class SimpleFormatter(logging.Formatter):
    def format(self, record):
        if (record.levelname == "WARNING") or (record.levelname == "ERROR") or (record.levelname == "CRITICAL"):
            log_message = f"{record.levelname} - {record.getMessage()}"
        else:
            log_message = f"... {record.getMessage()}"
        extra_info = {
            key: value for key, value in record.__dict__.items()
            if key not in EXCLUDED_KEYS and key not in ['levelname', 'message', 'name', 'pathname', 'lineno', 'levelno', 'exc_info', 'exc_text', 'args', 'msg', 'filename', 'module', 'funcName']
        }
        if extra_info:
            log_message += f" | {extra_info}"
        return log_message

def setup_logging(output_dir):
    stream_handler = StreamHandler()
    stream_handler.setFormatter(SimpleFormatter())

    # FileHandler to write logs to a JSON file
    file_handler = FileHandler(f"{output_dir}/execution_log.json")
    file_handler.setFormatter(JsonFormatter())

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    return logger

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
    dki_options.add_argument('-maxb', metavar=('<bmax>'),help='maximum b-value for DKI fitting, default=3.')

    smi_options = cmdline.add_argument_group('tensor options for the TMI script')
    smi_options.add_argument('-SMI', action='store_true',help='Perform estimation of SMI (standard model of Diffusion in White Matter). Please use in conjunction with the -bshape, -echo_time, -sigma, and -compartments options.')  
    smi_options.add_argument('-compartments', metavar=('<compartments>'),help='SMI compartments (IAS, EAS, and FW), default=IAS,EAS')  
    smi_options.add_argument('-sigma', metavar=('<noisemap>'),help='path to noise map for SMI parameter estimation')
    smi_options.add_argument('-lmax', metavar=('<lmax>'),help='lmax for SMI polynomial regression. must be 0,2,4, or 6.')

    #wmti_options = cmdline.add_argument_group('tensor options for the TMI script')
    #wmti_options.add_argument('-WMTI', action='store_true', help='Include WMTI parameters in output folder (awf,ias_params,eas_params)')

logger = None

def execute(): #pylint: disable=unused-variable
    from mrtrix3 import app, path, run, MRtrixError #pylint: disable=no-name-in-module, import-outside-toplevel
    import lib.tensor as tensor
    import numpy as np
    from ants import image_read
    import pandas as pd

    outdir = path.from_user(app.ARGS.output, True)
    if outdir[0]=="'" and outdir[-1]=="'":
        outdir= outdir.replace("'","")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        
    logger = setup_logging(outdir)
    logger.info(f"Output directory created: {outdir}")

    logger.info("Starting execution...")

    app.make_scratch_dir()

    dwi_metadata = get_input_info(
        app.ARGS.input, app.ARGS.fslbval, app.ARGS.fslbvec, app.ARGS.bids)
    logger.info("Input information obtained.", extra={"dwi_metadata": dwi_metadata})

    #assert_inputs(dwi_metadata, None, None)
    convert_input_data(dwi_metadata)
    logger.info("Input data converted successfully.")

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
    logger.info("DWI data converted to NIfTI format.")

    # nii = image_read('dwi.nii')
    # dwi = nii.numpy()
    bvec = np.loadtxt('dwi.bvec')
    bval = np.loadtxt('dwi.bval')
    logger.info("Loaded bvec and bval data.", extra={"bvec_shape": bvec.shape, "bval_shape": bval.shape})

    order = np.floor(np.log(abs(np.max(bval)+1)) / np.log(10))
    if order >= 2:
        bval = bval / 1000
        logger.info("bval converted normalized to ms/um^2.", extra={"bval": list(set(np.round(bval, 2)))})

    grad = np.hstack((bvec.T, bval[None,...].T))
    logger.info("Gradient table created.", extra={"grad_shape": grad.shape})

    if app.ARGS.mask:
        mask = image_read(path.from_user(app.ARGS.mask)).numpy()
        logger.info("Loaded mask from file.", extra={"mask_shape": mask.shape})
    else:
        mask = np.ones(dwi.shape[:-1])
        logger.info("No mask provided. Using default mask with all ones.", extra={"mask_shape": mask.shape})
    mask = mask.astype(bool)

    if app.ARGS.fit_constraints:
        constraints = app.ARGS.fit_constraints
        if type(constraints) == str:
                constraints = constraints.split(",")
                constraints = [int(i) for i in constraints]
        else:
            raise MRtrixError("Constraints must be a 3 element comma separated string (i.e. 0,1,0)")
        logger.info("Fit constraints provided.", extra={"constraints": constraints})
    else:
        constraints = [0,0,0]
        logger.info("No fit constraints provided. Using default constraints.", extra={"constraints": constraints})

    if (len(set(dwi_metadata['bshape_per_volume'])) > 1) or (len(set(dwi_metadata['echo_time_per_volume'])) > 1):
        multi_te_beta = True
        dwi_orig = dwi.copy()
        bval_orig = bval.copy()
        bvec_orig = bvec.copy()
        logger.info("Multi-TE or Multi-beta detected.", extra={"multi_te_beta": multi_te_beta})
    else:
        multi_te_beta = False


    if len(set(dwi_metadata['bshape_per_volume'])) > 1:
        if (app.ARGS.DTI) or (app.ARGS.DKI) or (app.ARGS.WDKI):
            logger.info("Variable b-shape detected, filtering for LTE volumes.")
            lte_idx = (dwi_metadata['bshape_per_volume'] == 1)
            if np.sum(lte_idx) < 6:
                logger.error("Not enough LTE measurements for DTI/DKI.", extra={"lte_count": np.sum(lte_idx)})
                raise MRtrixError("Not enough LTE measurements for DTI/DKI")
            dwi = dwi_orig[:,:,:,lte_idx]
            bval = bval_orig[lte_idx]
            bvec = bvec_orig[:,lte_idx]
            logger.info("Filtered DWI, bval, and bvec for LTE volumes.", extra={"dwi_shape": dwi.shape})


    if len(set(dwi_metadata['echo_time_per_volume'])) == 1:
        logger.info("Single echo time detected.")

        if app.ARGS.DTI:
            from lib.mpunits import vectorize
            logger.info("Starting DTI fit...")
            dtishell = (bval <= 0.1) | ((bval > .5) & (bval <= 1.5))
            dwi_dti = dwi[:,:,:,dtishell]
            bval_dti = bval[dtishell]
            bvec_dti = bvec[:,dtishell]
            grad_dti = np.hstack((bvec_dti.T, bval_dti[None,...].T))
            dti = tensor.TensorFitting(grad_dti, int(app.ARGS.n_cores))
            dt_dti, s0_dti, b_dti = dti.dti_fit(dwi_dti, mask)
            logger.info("DTI fit completed.", extra={"dt_dti_shape": dt_dti.shape})

            dt_ = {}
            dt_dti_ = vectorize(dt_dti, mask)
            dt_['dt'] = dt_dti_
            save_params(dt_, nii, model='dti', outdir=outdir)
            logger.info("DT saved.")

        if app.ARGS.DKI or app.ARGS.WDKI:
            from lib.mpunits import vectorize

            logger.info("Starting DKI fit...")

            if app.ARGS.maxb:
                maxb = float(app.ARGS.maxb)
                logger.info("Maximum b-value set.", extra={"maxb": maxb})
            else:
                maxb = 3.01

            dwi_dki = dwi[:,:,:,bval < maxb]
            bvec_dki = bvec[:,bval < maxb]
            bval_dki = bval[bval < maxb]

            bvec_dki_temp=bvec_dki[:,bval_dki > 0.1]
            ndir=np.shape(np.unique(bvec_dki_temp,axis=1))[1]

            if len(set(np.round(bval_dki, 2))) <= 2:
                logger.warning(f"Fewer than 2 nonzero shells found for DKI fitting. Skipping DKI/WDKI.", extra={"bval_dki": bval_dki.tolist()})
                app.ARGS.DKI=False
                app.ARGS.WDKI=False
                app.ARGS.akc_outliers=False
            elif ndir < 21:
                logger.warning("Fewer than 21 unique directions detected for DKI fitting. Skipping DKI/WDKI.")
                app.ARGS.DKI=False
                app.ARGS.WDKI=False
                app.ARGS.akc_outliers=False
            else:
                if np.max(bval_dki) < 2:
                    logger.warning("Max b-value is <2 for DKI fitting.")

                grad_dki = np.hstack((bvec_dki.T, bval_dki[None,...].T))
                dki = tensor.TensorFitting(grad_dki, int(app.ARGS.n_cores))
                dt_dki, s0_dki, b_dki = dki.dki_fit(dwi_dki, mask, constraints=constraints)
                logger.info("DKI fit completed.", extra={"dt_dki_shape": dt_dki.shape})

                dt_ = {}
                dt_dki_ = vectorize(dt_dki, mask)
                dt_['dt'] = dt_dki_
                save_params(dt_, nii, model='dki', outdir=outdir)
                logger.info("DKT saved.")

        if app.ARGS.polyreg:
            if app.ARGS.WDKI:
                dt_poly_dki = dki.train_rotated_bayes_fit(dwi_dki, dt_dki, s0_dki, b_dki, mask)
                logger.info("Polyreg WDKI fit completed.", extra={"dt_poly_dki_shape": dt_poly_dki.shape})
            if app.ARGS.DTI:
                dt_poly_dti = dti.train_rotated_bayes_fit(dwi_dti, dt_dti, s0_dti, b_dti, mask, True)
                logger.info("Polyreg DTI fit completed.", extra={"dt_poly_dti_shape": dt_poly_dti.shape})

        if app.ARGS.akc_outliers:
            from lib.mpunits import vectorize
            import scipy.io as sio

            logger.info("Starting AKC outlier detection...")
            dwd = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            mat = sio.loadmat(os.path.join(dwd,'constant','dirs10000.mat'))
            dir = mat['dir']

            if not (app.ARGS.DKI or app.ARGS.WDKI):
                logger.error("AKC Outlier detection must be accompanied by DKI option")
                raise MRtrixError("AKC Outlier detection must be accompanied by DKI option")
            else:
                akc_mask = dki.outlierdetection(dt_dki, mask, dir)
                
            akc_mask = vectorize(akc_mask, mask).astype(bool)
            logger.info("Outlier detection completed.", extra={"num_outliers": np.sum(akc_mask)})

            dwi_new = refit_or_smooth(akc_mask, dwi_dki, n_cores=int(app.ARGS.n_cores))
            dt_new,_,_ = dki.dki_fit(dwi_new, akc_mask)

            x,y,z = np.where(akc_mask == 1)
            DT = vectorize(dt_dki, mask)
            DT[x,y,z,:] = dt_new.T
            dt_dki = vectorize(DT, mask)
            akc_mask = dki.outlierdetection(dt_dki, mask, dir)
            akc_mask = vectorize(akc_mask, mask).astype(bool)
            logger.info("AKC outlier post-processing completed.", extra={"num_outliers": np.sum(akc_mask)})
        else:
            akc_mask = np.zeros_like(mask)

        if app.ARGS.fit_smoothing:
            logger.info("Starting tensor based smoothing...")
            if app.ARGS.DTI:
                dwi_new = refit_or_smooth(akc_mask, dwi_dti, mask=mask, smoothlevel=int(app.ARGS.fit_smoothing))
                dt_dti,_,_ = dti.dti_fit(dwi_new, mask)
                logger.info("DTI fit after smoothing completed.", extra={"dt_dti_shape": dt_dti.shape})
            if (app.ARGS.DKI or app.ARGS.WDKI):
                dwi_new = refit_or_smooth(akc_mask, dwi_dki, mask=mask, smoothlevel=int(app.ARGS.fit_smoothing))
                dt_dki,_,_ = dki.dki_fit(dwi_new, mask)
                logger.info("DKI fit after smoothing completed.", extra={"dt_dki_shape": dt_dki.shape})

        if app.ARGS.DTI or app.ARGS.DKI or app.ARGS.polyreg or app.ARGS.WDKI:
            rdwi = tensor.vectorize(dwi, mask)
            rdwi[~np.isfinite(rdwi)] = np.finfo(float).eps
            rdwi[rdwi<=0] = np.finfo(float).eps
            B = np.round(grad[:,-1]*1000)
            uB = np.unique(B)
            trace = np.zeros((rdwi.shape[1], uB.shape[0]))
            for ib in range(0, uB.shape[0]):
                t = np.where(B == uB[ib])[0]
                indices = np.where(t < rdwi.shape[0])
                t = t[indices]
                masked_data = np.ma.masked_where(rdwi[t,:] <= 0.1, rdwi[t,:])
                trace[:,ib] = np.exp(np.ma.mean(np.ma.log(masked_data), axis=0))
            trace = tensor.vectorize(trace.T, mask)
            params_trace = {'trace': trace}
            save_params(params_trace, nii, model='allshells', outdir=outdir)
            logger.info("Trace parameters saved for all shells.", extra={"trace_shape": trace.shape})

        if app.ARGS.DTI:
            logger.info("Extracting and saving DTI maps...")
            params_dti = dti.extract_parameters(dt_dti, b_dti, mask, extract_dti=True, extract_dki=False, fit_w=False)
            save_params(params_dti, nii, model='dti', outdir=outdir)
            logger.info("DTI maps saved.")

        if app.ARGS.DKI:
            logger.info("Extracting and saving DKI maps...")
            params_dki = dki.extract_parameters(dt_dki, b_dki, mask, extract_dti=True, extract_dki=True, fit_w=False)
            save_params(params_dki, nii, model='dki', outdir=outdir)
            logger.info("DKI maps saved.")

        if app.ARGS.polyreg:
            if app.ARGS.WDKI:
                logger.info("Extracting and saving polyreg WDKI maps...")
                params_dki_poly = dki.extract_parameters(dt_poly_dki, b_dki, mask, extract_dti=True, extract_dki=True, fit_w=True)
                save_params(params_dki_poly, nii, model='wdki_poly', outdir=outdir)
                logger.info("Polyreg WDKI maps saved.")
            if app.ARGS.DTI:
                logger.info("Extracting and saving polyreg DTI maps...")
                params_dti_poly = dti.extract_parameters(dt_poly_dti, b_dti, mask, extract_dti=True, extract_dki=False, fit_w=False)
                save_params(params_dti_poly, nii, model='dti_poly', outdir=outdir)
                logger.info("Polyreg DTI maps saved.")

        if app.ARGS.WDKI:
            logger.info("Extracting and saving WDKI maps...")
            params_dwi = dki.extract_parameters(dt_dki, b_dki, mask, extract_dti=True, extract_dki=True, fit_w=True)
            save_params(params_dwi, nii, model='wdki', outdir=outdir)
            logger.info("WDKI maps saved.")
    else:

        for te in set(dwi_metadata['echo_time_per_volume']):
            logger.info(f"Processing echo time {te}...")

            te_idx = (dwi_metadata['echo_time_per_volume'] == te) & (dwi_metadata['bshape_per_volume'] == 1)
            dwi = dwi_orig[:,:,:,te_idx]
            bval = bval_orig[te_idx]
            bvec = bvec_orig[:,te_idx]

            if app.ARGS.DTI:

                if np.sum(te_idx) < 6:
                    logger.warning(f"Fewer than 6 measurements found for TE={te}. Skipping DTI fit.")
                    continue
                
                logger.info(f"Starting DTI fit for TE={te}...")
                dtishell = (bval <= 0.1) | ((bval > .5) & (bval <= 1.5))
                dwi_dti = dwi[:,:,:,dtishell]
                bval_dti = bval[dtishell]
                bvec_dti = bvec[:,dtishell]
                grad_dti = np.hstack((bvec_dti.T, bval_dti[None,...].T))
                dti = tensor.TensorFitting(grad_dti, int(app.ARGS.n_cores))
                dt_dti, s0_dti, b_dti = dti.dti_fit(dwi_dti, mask)
                logger.info(f"DTI fit completed for TE={te}.", extra={"dt_dti_shape": dt_dti.shape})

            if app.ARGS.DKI or app.ARGS.WDKI:

                if np.sum(te_idx) < 21:
                    logger.warning(f"Fewer than 21 measurements found for TE={te}. Skipping DKI fit.")
                    continue

                logger.info(f"Starting DKI fit for TE={te} with constraints = {constraints}...")

                if app.ARGS.maxb:
                    maxb = float(app.ARGS.maxb)
                    logger.info(f"Maximum b-value set to {maxb} for TE={te}.")
                else:
                    maxb = 3.01

                dwi_dki = dwi[:,:,:,bval < maxb]
                bvec_dki = bvec[:,bval < maxb]
                bval_dki = bval[bval < maxb]

                bvec_dki_temp=bvec_dki[:,bval_dki > 0.1]
                ndir=np.shape(np.unique(bvec_dki_temp,axis=1))[1]

                if len(set(np.round(bval_dki, 2))) <= 2:
                    logger.warning(f"Fewer than 2 nonzero shells found for DKI fitting at TE={te}. Skipping DKI/WDKI.", extra={"bval_dki": bval_dki.tolist()})
                    app.ARGS.DKI=False
                    app.ARGS.WDKI=False
                    app.ARGS.akc_outliers=False
                elif ndir < 21:
                    logger.warning("Fewer than 21 unique directions detected for DKI fitting. Skipping DKI/WDKI.")
                    app.ARGS.DKI=False
                    app.ARGS.WDKI=False
                    app.ARGS.akc_outliers=False
                else:
                    if np.max(bval_dki) < 2:
                        logger.warning("Max b-value is <2 for DKI fitting.")

                    grad = np.hstack((bvec_dki.T, bval_dki[None,...].T))
                    dki = tensor.TensorFitting(grad, int(app.ARGS.n_cores))
                    dt_dki, s0_dki, b_dki = dki.dki_fit(dwi_dki, mask, constraints=constraints)
                    logger.info(f"DKI fit completed for TE={te}.", extra={"dt_dki_shape": dt_dki.shape})

                if app.ARGS.polyreg:
                    if app.ARGS.WDKI:
                        dt_poly_dki = dki.train_rotated_bayes_fit(dwi_dki, dt_dki, s0_dki, b_dki, mask)
                        logger.info(f"Polyreg WDKI fit completed for TE={te}.", extra={"dt_poly_dki_shape": dt_poly_dki.shape})
                    if app.ARGS.DTI:
                        dt_poly_dti = dti.train_rotated_bayes_fit(dwi_dti, dt_dti, s0_dti, b_dti, mask, True)
                        logger.info(f"Polyreg DTI fit completed for TE={te}.", extra={"dt_poly_dti_shape": dt_poly_dti.shape})

            if app.ARGS.akc_outliers:
                from lib.mpunits import vectorize
                import scipy.io as sio

                logger.info(f"Starting AKC outlier detection for TE={te}...")
                dwd = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
                mat = sio.loadmat(os.path.join(dwd,'constant','dirs10000.mat'))
                dir = mat['dir']

                if not (app.ARGS.DKI or app.ARGS.WDKI):
                    logger.error(f"AKC Outlier detection for TE={te} must be accompanied by DKI option.")
                    raise MRtrixError("AKC Outlier detection must be accompanied by DKI option")
                else:
                    akc_mask = dki.outlierdetection(dt_dki, mask, dir)
                    
                akc_mask = vectorize(akc_mask, mask).astype(bool)
                logger.info(f"Outlier detection completed for TE={te}.", extra={"num_outliers": np.sum(akc_mask)})

                dwi_new = refit_or_smooth(akc_mask, dwi_dki, n_cores=int(app.ARGS.n_cores))
                dt_new,_,_ = dki.dki_fit(dwi_new, akc_mask)

                x,y,z = np.where(akc_mask == 1)
                DT = vectorize(dt_dki, mask)
                DT[x,y,z,:] = dt_new.T
                dt_dki = vectorize(DT, mask)
                akc_mask = dki.outlierdetection(dt_dki, mask, dir)
                akc_mask = vectorize(akc_mask, mask).astype(bool)
                logger.info(f"AKC outlier post-processing completed for TE={te}.", extra={"num_outliers": np.sum(akc_mask)})
            else:
                akc_mask = np.zeros_like(mask)

            if app.ARGS.fit_smoothing:
                logger.info(f"Starting tensor based smoothing for TE={te}...")
                if app.ARGS.DTI:
                    dwi_new = refit_or_smooth(akc_mask, dwi_dti, mask=mask, smoothlevel=int(app.ARGS.fit_smoothing))
                    dt_dti,_,_ = dti.dti_fit(dwi_new, mask)
                    logger.info(f"DTI fit after smoothing completed for TE={te}.", extra={"dt_dti_shape": dt_dti.shape})
                if (app.ARGS.DKI or app.ARGS.WDKI):
                    dwi_new = refit_or_smooth(akc_mask, dwi_dki, mask=mask, smoothlevel=int(app.ARGS.fit_smoothing))
                    dt_dki,_,_ = dki.dki_fit(dwi_new, mask)
                    logger.info(f"DKI fit after smoothing completed for TE={te}.", extra={"dt_dki_shape": dt_dki.shape})

            if app.ARGS.DTI or app.ARGS.DKI or app.ARGS.polyreg or app.ARGS.WDKI:
                rdwi = tensor.vectorize(dwi, mask)
                B = np.round(grad[:,-1]*1000)
                uB = np.unique(B)
                trace = np.zeros((rdwi.shape[1], uB.shape[0]))
                rdwi[~np.isfinite(rdwi)] = np.finfo(float).eps
                rdwi[rdwi<=0] = np.finfo(float).eps
                for ib in range(0, uB.shape[0]):
                    t = np.where(B == uB[ib])[0]
                    indices = np.where(t < rdwi.shape[0])
                    t = t[indices]
                    masked_data = np.ma.masked_where(rdwi[t,:] <= 0.1, rdwi[t,:])
                    trace[:,ib] = np.exp(np.ma.mean(np.ma.log(masked_data), axis=0))
                trace = tensor.vectorize(trace.T, mask)
                params_trace = {'trace': trace}
                save_params(params_trace, nii, model='te'+str(te)+'_shells', outdir=outdir)
                logger.info(f"Trace parameters saved for TE={te}.")

            if app.ARGS.DTI:
                logger.info(f"Extracting and saving DTI maps for TE={te}...")
                params_dti = dti.extract_parameters(dt_dti, b_dti, mask, extract_dti=True, extract_dki=False, fit_w=False)
                save_params(params_dti, nii, model='dti_te'+str(te), outdir=outdir)
                logger.info(f"DTI maps saved for TE={te}.")

            if app.ARGS.DKI:
                logger.info(f"Extracting and saving DKI maps for TE={te}...")
                params_dki = dki.extract_parameters(dt_dki, b_dki, mask, extract_dti=True, extract_dki=True, fit_w=False)
                save_params(params_dki, nii, model='dki_te'+str(te), outdir=outdir)
                logger.info(f"DKI maps saved for TE={te}.")

            if app.ARGS.polyreg:
                if app.ARGS.WDKI:
                    logger.info(f"Extracting and saving polyreg WDKI maps for TE={te}...")
                    params_dki_poly = dki.extract_parameters(dt_poly_dki, b_dki, mask, extract_dti=True, extract_dki=True, fit_w=True)
                    save_params(params_dki_poly, nii, model='wdki_poly_te'+str(te), outdir=outdir)
                    logger.info(f"Polyreg WDKI maps saved for TE={te}.")
                if app.ARGS.DTI:
                    logger.info(f"Extracting and saving polyreg DTI maps for TE={te}...")
                    params_dti_poly = dti.extract_parameters(dt_poly_dti, b_dti, mask, extract_dti=True, extract_dki=False, fit_w=False)
                    save_params(params_dti_poly, nii, model='dti_poly_te'+str(te), outdir=outdir)
                    logger.info(f"Polyreg DTI maps saved for TE={te}.")

            if app.ARGS.WDKI:
                logger.info(f"Extracting and saving WDKI maps for TE={te}...")
                params_dwi = dki.extract_parameters(dt_dki, b_dki, mask, extract_dti=True, extract_dki=True, fit_w=True)
                save_params(params_dwi, nii, model='wdki_te'+str(te), outdir=outdir)
                logger.info(f"WDKI maps saved for TE={te}.")

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

        logger.info("Starting SMI fitting process...")
        bvec_temp=bvec[:,bval > 0.1]
        ndir=np.shape(np.unique(bvec_temp,axis=1))[1]

        if len(set(np.round(bval, 2))) <= 2:
            logger.warning(f"Fewer than 2 nonzero shells found for SMI fitting. Skipping SMI.")
            app.ARGS.SMI=False
        elif ndir < 21:
            logger.warning("Fewer than 21 unique directions detected for SMI fitting. Skipping SMI.")
            app.ARGS.SMI=False
        else:
            if np.max(bval) < 2:
                logger.warning("Max b-value is <2 for SMI fitting.")

            if not app.ARGS.sigma:
                logger.warning("No sigma map provided. SMI may be poorly conditioned.")
                sigma = None
            else:
                sigma = image_read(path.from_user(app.ARGS.sigma)).numpy()
                logger.info("Sigma map loaded.", extra={"sigma_shape": sigma.shape})

            if app.ARGS.compartments:
                compartments = app.ARGS.compartments
                if type(compartments) == str:
                    compartments = compartments.split(",")
                    compartments = [str(i) for i in compartments]
                else:
                    raise MRtrixError(" Compartments must be a comma sepearated string (i.e. IAS,EAS)")
                logger.info("SMI compartments specified.", extra={"compartments": compartments})
            else: 
                compartments = ['IAS', 'EAS']
                logger.info("Using default SMI compartments.", extra={"compartments": compartments})

            if app.ARGS.lmax:
                if int(app.ARGS.lmax) not in {0, 2, 4, 6}:
                    raise ValueError("lmax value must be 0, 2, 4, or 6.")
                else:
                    lmax = int(app.ARGS.lmax)
                logger.info("lmax for SMI specified.", extra={"lmax": lmax})
            else:
                lmax = None
                logger.info("No lmax specified for SMI. Using default.")

            logger.info("Initializing SMI fitting...")
            echo_times = dwi_metadata['echo_time_per_volume']
            if (np.min(echo_times) < 1.0) and (np.min(echo_times) > 0):
                logger.info("Echo times in s, converting to ms.")
                echo_times *= 1000

            if multi_te_beta:
                smi = SMI(bval=bval_orig, bvec=bvec_orig, rotinv_lmax=lmax,
                        compartments=compartments, echo_time=echo_times,
                        beta=dwi_metadata['bshape_per_volume'])
                
                logger.info("SMI model initialized for multi-TE/beta data.")
                
                params_smi = smi.fit(dwi_orig, mask=mask, sigma=sigma)
                logger.info("SMI fitting completed for multi-TE/beta data.", extra={"params_smi_shape": {key: value.shape for key, value in params_smi.items()}})
                
                save_params(params_smi, nii, model='smi', outdir=outdir)
                logger.info("SMI parameters saved for multi-TE/beta data.", extra={"outdir": outdir})
            else:
                smi = SMI(bval=bval, bvec=bvec, rotinv_lmax=lmax,
                        compartments=compartments, echo_time=echo_times,
                        beta=dwi_metadata['bshape_per_volume'])
                
                logger.info("SMI model initialized for single-TE/beta data.")

                params_smi = smi.fit(dwi, mask=mask, sigma=sigma)
                logger.info("SMI fitting completed for single-TE/beta data.", extra={"params_smi_shape": {key: value.shape for key, value in params_smi.items()}})

                save_params(params_smi, nii, model='smi', outdir=outdir)
                logger.info("SMI parameters saved for single-TE/beta data.", extra={"outdir": outdir})
    
    logger.info("Execution completed successfully.")

def main():
    import mrtrix3
    mrtrix3.execute() #pylint: disable=no-member

if __name__ == "__main__":
    main()


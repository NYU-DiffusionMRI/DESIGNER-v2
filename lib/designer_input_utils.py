"""
Utilities for checking compatibility of inputs to designer

Includes:
splitext_: a function to split .gz compressed images

get_input_info: a function to get datatype, gradient info, and image metadata for each input series

assert_inputs: asserts that any manual inputs match those found in the bids json, if avilable

convert_input_data: Merge input series and gradient info, setting up the working pipeline
"""

def splitext_(path):
    import os
    for ext in ['.tar.gz', '.tar.bz2','.nii.gz']:
        if path.endswith(ext):
            return path[:-len(ext)], path[-len(ext):]
    return os.path.splitext(path)

def get_input_info(input, fslbval, fslbvec, bids):
    import os
    from mrtrix3 import image, MRtrixError

    UserCpath = input.rsplit(',')
    DWIlist = [os.path.realpath(i) for i in UserCpath]
    
    isdicom = False
    for i in DWIlist:
        if not os.path.exists(i):
            print('cannot find input ' + i)
            quit()
        if os.path.isdir(i):
            format = image.Header(i).format()
            if format == 'DICOM':
                isdicom = True
            else:
                raise MRtrixError('input is a directory but does not contain DICOMs, quitting')

    DWIflist = [splitext_(i) for i in DWIlist]
    DWInlist = [i[0] for i in DWIflist]
    DWIext = [i[1] for i in DWIflist]

    if not fslbval:
        bvallist = [i + '.bval' for i in DWInlist]
    else:
        UserBvalpath = fslbval.rsplit(',')
        bvallist = [os.path.realpath(i) for i in UserBvalpath]
    if not fslbvec:
        bveclist = [i + '.bvec' for i in DWInlist]
    else:
        UserBvecpath = fslbvec.rsplit(',')
        bveclist = [os.path.realpath(i) for i in UserBvecpath]

    if not bids:
        bidslist = [i + '.json' for i in DWInlist]
    else:
        UserBidspath = bids.resplit(',')
        bidslist = [os.path.realpath(i) for i in UserBidspath]

    dwi_metadata = {
        'isdicom': isdicom, 
        'dwi_ext': DWIext,
        'bvecs': bveclist,
        'bvals': bvallist,
        'dwi_list': DWInlist,
        'bidslist': bidslist}

    return dwi_metadata

def convert_pe_dir_to_ijk(args_pe_dir):
    from mrtrix3 import MRtrixError

    avail_pe_dirs = [None,'1','-1','2','-2','3','-3','i','i-','j','j-','k','k-','LR','RL','AP','PA','IS','SI']

    if not args_pe_dir in avail_pe_dirs:
        raise MRtrixError('unsupported PE direction - choose one of (1,i,LR,-1,i-, etc...')

    if args_pe_dir == '1':
        pe_dir_app = 'i'
    elif args_pe_dir == '2':
        pe_dir_app = 'j'
    elif args_pe_dir == '3':
        pe_dir_app = 'k'
    elif args_pe_dir == '-1':
        pe_dir_app = 'i-'
    elif args_pe_dir == '-2':
        pe_dir_app = 'j-'
    elif args_pe_dir == '-3':
        pe_dir_app = 'k-'
    elif args_pe_dir == 'RL':
        pe_dir_app = 'i-'
    elif args_pe_dir == 'LR':
        pe_dir_app = 'i'
    elif args_pe_dir == 'AP':
        pe_dir_app = 'j-'
    elif args_pe_dir == 'PA':
        pe_dir_app = 'j'
    elif args_pe_dir == 'IS':
        pe_dir_app = 'k'
    elif args_pe_dir == 'SI':
        pe_dir_app = 'k-'
    elif args_pe_dir == 'i':
        pe_dir_app = 'i'
    elif args_pe_dir == 'j':
        pe_dir_app = 'j'
    elif args_pe_dir == 'k':
        pe_dir_app = 'k'
    elif args_pe_dir == 'i-':
        pe_dir_app = 'i-'
    elif args_pe_dir == 'j-':
        pe_dir_app = 'j-'
    elif args_pe_dir == 'k-':
        pe_dir_app = 'k-'
    else:
        pe_dir_app = None

    return pe_dir_app

def assert_inputs(dwi_metadata, args_pe_dir, args_pf):
    import json
    from mrtrix3 import app, MRtrixError

    bidslist = dwi_metadata['bidslist']
    try:
        bids = [json.load(open(i)) for i in bidslist]
    
        TE_bids = []
        pe_dir = []
        pf = []
        for i in bids:
            TE_bids.append(i['EchoTime'])
            pe_dir.append(i['PhaseEncodingDirection'])
            try:
                pf.append(i['PartialFourier'])
            except:
                pf.append(1)

        if not all(x == pf[0] for x in pf):
            raise MRtrixError('input series have different partial fourier factors,', +
                            'series should be processed separately')
        pf_bids = pf[0]

        if not all(x == pe_dir[0] for x in pe_dir):
            raise MRtrixError('input series have different phase encoding directions, series should be processed separately')
        pe_dir_bids = pe_dir[0]
        pe_dir_app = convert_pe_dir_to_ijk(args_pe_dir)

        if app.ARGS.echo_time:
            TE_app = [float(i) for i in app.ARGS.echo_time.rsplit(',')]
        else:
            TE_app = [0] * len(bidslist)
    except:
        print('no bids files identified')
        pe_dir_app = convert_pe_dir_to_ijk(args_pe_dir)
        pe_dir_bids = None
        pf_bids = None
        TE_bids = None
        if app.ARGS.echo_time:
            TE_app = [float(i) for i in app.ARGS.echo_time.rsplit(',')]
        else:
            TE_app = [0] * len(dwi_metadata['dwi_list'])

    if app.ARGS.bshape:
        bshape = [float(i) for i in app.ARGS.bshape.rsplit(',')]
    else:
        bshape = [1] * len(dwi_metadata['dwi_list'])

    if (TE_app) and (TE_bids) and (TE_bids != TE_app):
        raise MRtrixError('User defined echo times do not match those found in bids .json, please check input data for consistancy')
    elif TE_app:
        TE = TE_app
    elif TE_bids:
        TE = TE_bids
    else:
        TE = [0] * len(dwi_metadata['dwi_list'])

    # if no partial fourier information is found, assume full sampling
    if not args_pf and not pf_bids:
        args_pf = 1
        pf_bids = 1
    
    if (pe_dir_bids and pe_dir_app) and (pe_dir_app != pe_dir_bids):
        raise MRtrixError('input phase encoding direction and phase encoding direction from bids file do not match')
    elif pe_dir_bids:
        pe_dir = pe_dir_bids
    else:
        pe_dir = pe_dir_app

    if (pf_bids and args_pf) and (float(args_pf) != pf_bids):
        raise MRtrixError('input partial fourier fractor and bids PF factor do not match')
    elif pf_bids:
        pf = pf_bids
    else:
        pf = args_pf
    
    def convert_to_float(frac_str):
        try:
            return float(frac_str)
        except ValueError:
            num, denom = frac_str.split('/')
            try:
                leading, num = num.split(' ')
                whole = float(leading)
            except ValueError:
                whole = 0
            frac = float(num) / float(denom)
            return whole - frac if whole < 0 else whole + frac

    dwi_metadata['pe_dir'] = pe_dir
    dwi_metadata['pf'] = convert_to_float(pf)
    dwi_metadata['TE'] = TE
    dwi_metadata['bshape'] = bshape

def convert_input_data(dwi_metadata):
    """
    convert input data to .mif mrtrix format
    concatenates all inputs along 4th dimension
    
    """
    from mrtrix3 import run, image, MRtrixError, app
    import numpy as np

    miflist = []
    idxlist = []
    telist = []
    bshapelist = []
    dwi_ind_size = [[0,0,0,0]]

    dwi_n_list = dwi_metadata['dwi_list']
    isdicom = dwi_metadata['isdicom']
    bveclist = dwi_metadata['bvecs']
    bvallist = dwi_metadata['bvals']
    dwi_ext = dwi_metadata['dwi_ext']
    te_per_series = dwi_metadata['TE']
    bshape_per_series = dwi_metadata['bshape']

    if len(dwi_n_list) == 1:
        if not isdicom:
            cmd = ('mrconvert -fslgrad %s %s %s%s %s/dwi.mif' % 
                (bveclist[0], bvallist[0], ''.join(dwi_n_list), ''.join(dwi_ext), app.SCRATCH_DIR))
            run.command(cmd)
        else:
            cmd = ('mrconvert %s %s/dwi.mif' %
            (''.join(dwi_n_list), app.SCRATCH_DIR))
            run.command(cmd)
        dwi_header = image.Header('%s/dwi.mif' % (app.SCRATCH_DIR))
        dwi_ind_size.append([ int(s) for s in dwi_header.size() ])

    else:
        for idx,i in enumerate(dwi_n_list):
            if not isdicom:
                cmd = ('mrconvert -fslgrad %s %s %s%s %s/dwi%s.mif' % 
                (bveclist[idx], bvallist[idx], i, dwi_ext[idx], app.SCRATCH_DIR, str(idx)))
                run.command(cmd)
            else:
                cmd = ('mrconvert %s %s/dwi%s.mif' %
                (i, app.SCRATCH_DIR, str(idx)))
                run.command(cmd)
            dwi_header = image.Header('%s/dwi%s.mif' % (app.SCRATCH_DIR, str(idx)))
            dwi_ind_size.append([ int(s) for s in dwi_header.size() ])
            miflist.append('%s/dwi%s.mif' % (app.SCRATCH_DIR, str(idx)))

        DWImif = ' '.join(miflist)
        cmd = ('mrcat -axis 3 %s %s/dwi.mif' % (DWImif, app.SCRATCH_DIR))
        run.command(cmd)

    # get diffusion header info - check to make sure all values are valid for processing
    dwi_header = image.Header('%s/dwi.mif' % (app.SCRATCH_DIR))
    dwi_size = [ int(s) for s in dwi_header.size() ]
    grad = dwi_header.keyval()['dw_scheme']
    grad = [ line for line in grad ]
    grad = [ [ float(f) for f in line ] for line in grad ]
    
    grad = np.array(grad)
    grad[:,-1] = grad[:,-1] / 1000
    
    #stride = dwi_header.strides()
    num_volumes = 1
    if len(dwi_size) == 4:
        num_volumes = dwi_size[3]

    nvols = [i[3] for i in dwi_ind_size]
    for idx,i in enumerate(dwi_n_list):
        if len(dwi_n_list) == 1:
            tmpidxlist = range(0,num_volumes)
        else:
            tmpidxlist = range(sum(nvols[:idx+1]),sum(nvols[:idx+1])+nvols[idx+1])
        idxlist.append(','.join(str(i) for i in tmpidxlist))
        telist.append([te_per_series[idx]] * nvols[idx+1])
        bshapelist.append([bshape_per_series[idx]] * nvols[idx+1])

    telist = [item for sublist in telist for item in sublist]
    bshapelist = [item for sublist in bshapelist for item in sublist]

    dwi_metadata['idxlist'] = idxlist
    dwi_metadata['echo_time_per_volume'] = np.array(telist)
    dwi_metadata['bshape_per_volume'] = np.array(bshapelist)
    dwi_metadata['grad'] = grad

    if grad is None:
        raise MRtrixError('No diffusion gradient table found')
    if not len(grad) == num_volumes:
        raise MRtrixError('Number of lines in gradient table (%s) does not match input image (%s volumes); check your input data' % 
        (str(len(grad)), str(num_volumes)))

    run.command('mrconvert %s/dwi.mif %s/working.mif' % (app.SCRATCH_DIR, app.SCRATCH_DIR), show=False)

def create_shell_table(dwi_metadata):
    from lib.smi import SMI

    bvals = dwi_metadata['grad'][:,-1]
    bvecs = dwi_metadata['grad'][:,:-1]

    smi = SMI(bval=bvals, bvec=bvecs)

    bshape = dwi_metadata['bshape_per_volume']
    echo_time = dwi_metadata['echo_time_per_volume']
    smi.set_bshape(bshape)
    smi.set_echotime(echo_time)

    return smi.group_dwi_in_shells_b_beta_te()

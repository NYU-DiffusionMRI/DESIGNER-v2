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

    return isdicom, DWIext, bveclist, bvallist, DWInlist, bidslist

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
        pe_dir_app == 'i'
    elif args_pe_dir == 'j':
        pe_dir_app == 'j'
    elif args_pe_dir == 'k':
        pe_dir_app == 'k'
    elif args_pe_dir == 'i-':
        pe_dir_app == 'i-'
    elif args_pe_dir == 'j-':
        pe_dir_app == 'j-'
    elif args_pe_dir == 'k-':
        pe_dir_app == 'k-'
    else:
        pe_dir_app = None

    return pe_dir_app

def assert_inputs(bidslist, args_pe_dir, args_pf):
    import json
    from mrtrix3 import app, MRtrixError

    try:
        bids = [json.load(open(i)) for i in bidslist]
    
        TE = []
        pe_dir = []
        pf = []
        for i in bids:
            TE.append(i['EchoTime'])
            pf.append(i['PartialFourier'])
            pe_dir.append(i['PhaseEncodingDirection'])

        if not all(x == TE[0] for x in TE):
            raise MRtrixError('input series have different echo times. This is not yet supported')
        else:
            TE = TE=[0]

        if not all(x == pf[0] for x in pf):
            raise MRtrixError('input series have different partial fourier factors,', +
                            'series should be processed separately')
        pf_bids = pf[0]

        if not all(x == pe_dir[0] for x in pe_dir):
            raise MRtrixError('input series have different phase encoding directions, series should be processed separately')
        pe_dir_bids = pe_dir[0]
    except:
        print('no bids files identified')
        pe_dir_app = convert_pe_dir_to_ijk(args_pe_dir)
        pe_dir_bids = None
        pf_bids = None
    
    if pe_dir_bids and pe_dir_app != pe_dir_bids:
        raise MRtrixError('input phase encoding direction and phase encoding direction from bids file do not match')
    elif pe_dir_bids:
        pe_dir = pe_dir_bids
    else:
        pe_dir = pe_dir_app

    if pf_bids and float(args_pf) != pf_bids:
        raise MRtrixError('input partial fourier fractor and bids PF factor do not match')
    elif pf_bids:
        pf = pf_bids
    else:
        pf = args_pf

    return pe_dir, pf

def convert_input_data(isdicom, DWIext, bveclist, bvallist, DWInlist):
    from mrtrix3 import run, image, MRtrixError, app

    miflist = []
    idxlist = []
    dwi_ind_size = [[0,0,0,0]]

    if len(DWInlist) == 1:
        if not isdicom:
            cmd = ('mrconvert -fslgrad %s %s %s%s %s/dwi.mif' % 
                (bveclist[0], bvallist[0], ''.join(DWInlist), ''.join(DWIext), app.SCRATCH_DIR))
            run.command(cmd)
        else:
            cmd = ('mrconvert %s %s/dwi.mif' %
            (''.join(DWInlist), app.SCRATCH_DIR))
            run.command(cmd)
    else:
        for idx,i in enumerate(DWInlist):
            if not isdicom:
                cmd = ('mrconvert -fslgrad %s %s %s%s %s/dwi%s.mif' % 
                (bveclist[idx], bvallist[idx], i, DWIext[idx], app.SCRATCH_DIR, str(idx)))
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
    
    #stride = dwi_header.strides()
    num_volumes = 1
    if len(dwi_size) == 4:
        num_volumes = dwi_size[3]
    bval = [int(i[3]) for i in grad]

    nvols = [i[3] for i in dwi_ind_size]
    for idx,i in enumerate(DWInlist):
        if len(DWInlist) == 1:
            tmpidxlist = range(0,num_volumes)
        else:
            tmpidxlist = range(sum(nvols[:idx+1]),sum(nvols[:idx+1])+nvols[idx+1])
        idxlist.append(','.join(str(i) for i in tmpidxlist))

    if not grad:
        raise MRtrixError('No diffusion gradient table found')
    if not len(grad) == num_volumes:
        raise MRtrixError('Number of lines in gradient table (%s) does not match input image (%s volumes); check your input data' % 
        (str(len(grad)), str(num_volumes)))

    run.command('mrconvert %s/dwi.mif %s/working.mif' % (app.SCRATCH_DIR, app.SCRATCH_DIR), show=False)
    
    return idxlist
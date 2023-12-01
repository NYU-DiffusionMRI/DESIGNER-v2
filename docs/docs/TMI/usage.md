---
layout: default
title: usage
has_children: false
nav_order: 5
parent: TMI
---

# Options and Usage for TMI
{: .no_toc }

This page describes each optional argument to TMI. For specific use-case examples of TMI calls, please see the [examples]({{ site.baseurl }}{% link docs/TMI/examples.md %}) page.

Main usage:
`tmi <input> <output>`

### `input`
Input to `tmi` can be any diffusion MRI image file that is compatible with mrtrix3 as an input. Ideally the input to `tmi` is the output from the `deisgner` preprocessing pipeline.

### `output`
Name of the folder which will contain parameter images.

By default, if none of the below options are used, TMI will not estimate parameters for the input data. TMI uses `mrconvert` from mrtrix3 as its image conversion workhorse and can therefore handle any input image data type that mrtrix3 can, including .nii, .nii.gz, .mgh, .mgz, .mif, .mif.gz, .img/.hdr, or dicom inputs (when the inputs are folder names that contain dicom data). 

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

## Options for tensor estimation

### `-DTI`
- Perform DTI fitting using weighted linear least squares.
- Include diffusion tensor parameters in output folder.
- Parameters extracted from the kurtosis tensor are labelled as `*_dti.nii` (e.g. ad_dti.nii, fa_dti.nii).
- When the DTI option is used, b-shells greater than 1000 are automatically excluded from tensor estimation.

### `-DKI`
- Perform DKI fitting using weighted linear least squares.
- Include diffusion kurtosis parameters in output folder.
- Parameters extracted from the diffusion tensor are labelled as `*_dki.nii` (e.g. ad_dki.nii, fa_dki.nii, mk_dki.nii).
- When DKI option is used, b-shells greater than 3000 are automatically excluded from tensor estimation.

### `-WDKI`
- Include diffusion kurtosis parameters based on the W definition of the kurtosis tensor in output folder.
- Parameters extracted from the W tensor are labelled as `*_wdki.nii` (e.g. ad_wdki.nii, fa_wdki.nii, mk_wdki.nii).

### `-fit_constraints <0,0,0>`
- For use with `-DKI` or `-WDKI`
- A comma separated list of either 0 or 1. Each index corresponds to a hard bound on the kurtosis tensor. Index 1 corresponds to diffusion tensor > 0, index 2 corresponds to kurtosis tensor > 0, index 3 corresponds to kurtosis tensor < 10. An input of `1,1,1` will perform a fully constrained DKI fit, and `0,0,0` will perform an unconstrained fit.
- Default `0,0,0`

## Options for outlier replacement

### `-akc_outliers`
- Brute force K tensor outlier detection and filtering. The kurtosis tensor is projected onto a sphere with 10,000 directions. Voxels where mean kurtosis is less than -1 or greater than 10 are labelled as outliers.
- Outliers are replaced by the median value of neighboring (26 adjacent) non-outlier voxels in all diffusion weighted images and the DKI fit is done again.

### `-fit_smoothing <percentile>`
- Windowed adaptive nonlocal means filter on the dwi.
- Diffusion weighted images undergo an NLM filter over the most similar neighboring pixels and directions prior to fitting.
- For each voxel we take a cubic patch in 4d and take the mean over the $n$ most similar voxels in *xyz* and *q-space*. n is decided by the input argument to tmi, where with n=10 (default) we are smoothing over the 10% most similar voxels, and so on.
- If used in conjunctions with `-akc_outliers` outliers are detected and ignored during the smoothing process.
- The `<percentile>` Defines the degree of similar voxels in a sliding window patch to average. A percentile of 10 will average the 10% most similar voxels and directions in the sliding window. We recommend a percentile of 10 as a starting point for tuning this option.

## Options for model estimation

### `-SMI`
- Perform estimation of the Standard Model of WHite Matter.
- Use in conjunction with `-compartments`, `-sigma`, `-bshape`, and `-echo_time` options.
- We recommend prior estimation of the noise level sigma using `designer` along with the `-denoise` option.

### `-WMTI`
- Include WMTI parameters in output folder (awf,ias_params,eas_params).
- WMTI estimation is performed using the `dipy` package.

### `-compartments <IAS,EAS,FW>`
- SMI compartments (IAS, EAS, and FW) expressed as a comma separated string.
- default: `IAS,EAS`

### `-sigma <noisemap>`
- Path to noise map for SMI parameter estimation. Not required but recommended.
- We recommend computing sigma prior to running `tmi` using `designer`.

### `-bshape <beta1,beta2,...>`
- Specify the b-shape used in the acquisition (comma separated list the same length as number of inputs).

### `-echo_time <TE1,TE2,...>`
- Specify the echo time used in the acquisition (comma separated list the same length as number of inputs).

## Other options for TMI

### `-mask <mask>`
- Perform parameter estimation within the specified brain mask.

### `-datatype <dtype>`
- Specify the output datatype. Valid options are float32, float32le, float32be, float64, float64le, float64be, int64, uint64, int64le, uint64le, int64be, uint64be, int32, uint32, int32le, uint32le, int32be, uint32be, int16, uint16, int16le, uint16le, int16be, uint16be, cfloat32, cfloat32le, cfloat32be, cfloat64, cfloat64le, cfloat64be, int8, uint8, bit.

### `-fslbvec <bvecs>`
- Specify bvec path if path is different from the path to the dwi(s) or the file has an unusual extension.

### `-fslbval <bvals>`
- Specify bval path if path is different from the path to the dwi(s) or the file has an unusual extension.

### `-bids <bids>`
- Specify BIDS .json path(s) if different from the path to the dwi or the file has an unusual extension.

### `-n_cores <ncpus>`
- Specify the number of cores to use in parallel tasks, by default TMI will use total available cores - 2.
- Parallel tasks include `-denoise`, `-degibbs`, and `-eddy` (when `eddy_openmp` is the `eddy` backend).

### `-nocleanup`
- Do not delete intermediate files during script execution, and do not delete scratch directory at script completion.

### `-scratch </path/to/scratch>`
- Manually specify the path in which to generate the scratch directory.

---

## Standard options for Mrtrix3 executables

### `-info`
- Display information messages.

### `-quiet`
- Do not display information messages or progress status. Alternatively, this can be achieved by setting the MRTRIX_QUIET environment variable to a non-empty string.

### `-debug`
- Display debugging messages.

### `-force`
- Force overwrite of output files.

### `-nthreads`
- Use this number of threads in multi-threaded applications (set to 0 to disable multi-threading).
- Not used in Designer, use `-n_cores` instead.

### `-config key value`
- Temporarily set the value of an MRtrix config file entry.
- Not used in Designer.

### `-help`
- Display information page and exit.







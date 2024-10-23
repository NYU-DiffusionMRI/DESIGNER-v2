---
layout: default
title: usage
has_children: false
nav_order: 2
parent: designer
---

# Options and Usage for Designer
{: .no_toc }

This page describes each optional argument to designer. For specific use-case examples of designer calls, please see the [examples]({{ site.baseurl }}{% link docs/designer/examples.md %}) page.

Main usage:
`designer <input1,input2,...> <output>`

### `input`
- Inputs are a comma separated list of diffusion weighted image files. 
- Designer uses `mrconvert` from mrtrix3 as its image conversion workhorse and can therefore handle any input image data type that mrtrix3 can, including .nii, .nii.gz, .mgh, .mgz, .mif, .mif.gz, .img/.hdr, or dicom inputs (when the inputs are folder names that contain DICOM data). 

### `output`
- Name of the output image in any of the above formats (except DICOM)
- By default, if none of the below options are used, designer will not preprocess the input data. It will simply concatenate the comma separated inputs and output the concatenated dataset. 
- Designer will also output gradient directions and magnitudes in FSL (.bvec, .bval) format according to the naming convention for the output file.

{: .ref }
> Ades-Aron, B., Veraart, J., Kochunov, P., McGuire, S., Sherman, P., Kellner, E., ... & Fieremans, E. (2018). Evaluation of the accuracy and precision of the diffusion parameter EStImation with Gibbs and NoisE removal pipeline. Neuroimage, 183, 532-543.
>
> Chen, J., Ades-Aron, B., Lee, H.H., Mehrin, S., Pang, M., Novikov, D. S., Veraart, J., & Fieremans, E. (2024). Optimization and validation of the DESIGNER preprocessing pipeline for clinical diffusion MRI in white matter aging. Imaging Neuroscience, 2 1–17

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

## Options relating to denoising

### `-denoise`
- Performs MPPCA magnitude or complex (see -phase option) denoising using default parameters on input data
- Note that designer always performs MPPCA denoising as the first step in the processing pipeline (if the `-denoise` option is used) in order to best preserve the independent noise statistics of input data.
- By default designer chooses 5x5x5 cube shaped patch.
- When using the `-denoise` option please use the following citation:

{: .ref}
Veraart, J., Novikov, D. S., Christiaens, D., Ades-Aron, B., Sijbers, J., & Fieremans, E. (2016). Denoising of diffusion MRI using random matrix theory. Neuroimage, 142, 394-406.

### `-patch2self`
- Performs patch2self denoising (from dipy) on input data.
- Patch2self is included in Designer as a denoising alternative for datasets where noise has been interpolated, such as on data from GE or Phillips scanners.
- Patch2self option should not be used in conjunction with the `-denoise` option.
- When using this option please include the following citation:

{: .ref }
Fadnavis, S., Batson, J., & Garyfallidis, E. (2020). Patch2Self: Denoising Diffusion MRI with Self-Supervised Learning​. Advances in Neural Information Processing Systems, 33, 16293-16303.

### `-shrinkage <frob>`
- Perform or do not perform singular value shrinkage using the Frobenuous norm described [here]({{ site.baseurl }}{% link docs/designer/background.md %}#singular-value-shrinkage) after MPPCA denoising
- Options are `threshold` or `frob`. default `frob`.
- Used in conjunction with the `-denoise` option.
- When using this option please include the following citation:

{: .ref }
Gavish, M., & Donoho, D. L. (2017). Optimal shrinkage of singular values. IEEE Transactions on Information Theory, 63(4), 2137-2152.

### `-algorithm <cordero-grande>`
- Specification for the MP thresholding algorithm described [here]({{ site.baseurl }}{% link docs/designer/background.md %}#symmetric-pca-thresholding).
- Used in conjunction with the `-denoise` option.
- Options are `veraart`, `cordero-grande`, and `jespersen`; default `cordero-grande`.
- When using the `cordero-grande` algorithm please using the following citation:

{: .ref }
Cordero-Grande, L., Christiaens, D., Hutter, J., Price, A. N., & Hajnal, J. V. (2019). Complex diffusion-weighted image estimation via matrix recovery under general noise models. Neuroimage, 200, 391-404.

-when using the `jespersen` algorithm please use the following citation:

{: .ref }
Olesen, J. L., Ianus, A., Østergaard, L., Shemesh, N., & Jespersen, S. N. (2023). Tensor denoising of multidimensional MRI data. Magnetic Resonance in Medicine, 89(3), 1160-1172.

### `-extent <x,y,z>`
- Manually specify the patch extent for MPPCA.
- Used in conjunction with `denoise` option.
- Inputs should be 3 integers separated by commas, e.g. `5,5,5` for a patch size of $5 \times 5 \times 5$, or `11,11,1` for a two dimensional $11 \times 11$ patch. 
- Users should ensure that their chosen patch size is not larger than the total image size in any dimension.

### `-adaptive_patch`
- Run MPPCA with adaptive patching as described [here]({{ site.baseurl }}{% link docs/designer/background.md %}#adaptive-patching).
- used in conjunction with the `-denoise` option.
- When this option is used an adaptive patch is chosen where N=100
- Note when this argument is used, the patch size will be chosen based on the extent, where the total number of voxels in the adaptive patch is 80% of the total voxels in a patch.
- Without the `-adaptive_patch` option, the patches are chosen as blocks with a step size of floor(extent/2), this increases speed at the expense of accuracy in noise estimation.

### `-phase <phase_image1,phase_image2,...>`
- Include a volume (or volumes) of phase images corresponding to the input(s) to Designer. For example, if two series are input to Designer with a total of 100 directions, there should be two phase series with 100 total volumes.
- Used in conjunction with `-denoise` option.
- This option will perform denoising on the complex-valued DWI data to help reduce the effect of the noise floor as described [here]({{ site.baseurl }}{% link docs/designer/background.md %}#denoising-complex-data).

---

## Options relating to EPI distortion correction, eddy current correction, and motion correction

### `-eddy`
- Perform eddy current and motion correction using FSL eddy. This option requires the use of additional options (`rpe_*`) found below. For more information on the `-rpe_*` options please visit the mrtrix3 [dwifslpreproc](https://mrtrix.readthedocs.io/en/dev/reference/commands/dwifslpreproc.html?highlight=dwifslpreproc) documentation page.
- Users should note that Designer always uses certain eddy arguments including `-repol` slice-wise outlier replacement and `-data_is_shelled` to inform `eddy` that the data is multi-shell diffusion.
- Brain masks are always computed before and after eddy to improve performance.
- When using the `-eddy` option please use the following citations:

{: .ref }
> Andersson, J. L. & Sotiropoulos, S. N. An integrated approach to correction for off-resonance effects and subject movement in diffusion MR imaging. NeuroImage, 2015, 125, 1063-1078
>
> Smith, S. M.; Jenkinson, M.; Woolrich, M. W.; Beckmann, C. F.; Behrens, T. E.; Johansen-Berg, H.; Bannister, P. R.; De Luca, M.; Drobnjak, I.; Flitney, D. E.; Niazy, R. K.; Saunders, J.; Vickers, J.; Zhang, Y.; De Stefano, N.; Brady, J. M. & Matthews, P. M. Advances in functional and structural MR image analysis and implementation as FSL. NeuroImage, 2004, 23, S208-S219
>
> Andersson, J. L. R.; Graham, M. S.; Zsoldos, E. & Sotiropoulos, S. N. Incorporating outlier detection and replacement into a non-parametric framework for movement and distortion correction of diffusion MR images. NeuroImage, 2016, 141, 556-572
>
> If performing recombination of diffusion-weighted volume pairs with opposing phase encoding directions: Skare, S. & Bammer, R. Jacobian weighting of distortion corrected EPI data. Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 5063
>
> If performing EPI susceptibility distortion correction: Andersson, J. L.; Skare, S. & Ashburner, J. How to correct susceptibility distortions in spin-echo echo-planar images: application to diffusion tensor imaging. NeuroImage, 2003, 20, 870-888.

### `-eddy_groups <index1,index2,...>`
- Specify how input series should be grouped when running eddy, for use with variable TE or b-shape data. 
- Inputs should be a comma separated list of integers beginning with 1, i.g. 1,1,1,2,2,2 for 6 series where the first 3 series and second 3 series have different echo times.
- This command informs Designer of how to concatenate data for input to eddy. For data with different echo times or b-shapes, the data with the same echo-time/b-shape will be stacked together and input to eddy. Afterwards, each stacked input to eddy is rigidly aligned as a separate step.
- Use this command in conjunction with `-eddy`, and (`-bshape` or `-echo_time`).

### `-rpe_none`
- Specify that no reversed phase-encoding image data is being provided; eddy will perform eddy current and motion correction only and `topup` will not be performed.
- Use in conjunction with `-eddy`.

### `-rpe_pair <reverse_phase_encoding_b0>`
- Specify the reverse phase encoding b=0 image.
- Use in conjunction with `-eddy`
- Prior to input to `topup`, Designer will perform a rigid registration between the input reverse phase encoded image and the mean b=0 image from the input DWI dataset.

### `-rpe_all <reverse_phase_encoding_dataset>`
- Specify that all DWIs have been acquired with opposing phase-encoding; this information will be used to perform a recombination of image volumes (each pair of volumes with the same b-vector but different phase encoding directions will be combined together into a single volume). 
- The argument to this option is the set of volumes with reverse phase encoding but the same b-vectors as the input image.
- Use in conjunction with `-eddy`.
- *This option has not been tested thoroughly for varying combinations of echo times and b-shapes*.

### `-rpe_header`
- Specify that the phase-encoding information can be found in the image header(s), and that this is the information that the script should use.
- Use in conjunction with `-eddy`.

### `-rpe_te <echo_time (milliseconds)>`
- Specify the echo time of the reverse phase encoded image, if it is not accompanied by a BIDS .json sidecar. If there is a BIDS format sidecar and this option is used, the two will be compared to ensure the correct echo-time was input. 
- For use in conjunction with `-eddy` and `-eddy_groups`.
- This option is only necessary when there are multiple inputs to Designer with differing echo times. It ensures that `topup` is run using forward and reverse phase encoding images with the same TE.

### `-pe_dir <phase_encoding_direction>`
- Specify the phase encoding direction of the main inputs to Designer. All inputs should be acquired with the same phase encoding direction.
- Use in conjunction with `eddy` or `degibbs`.
- By default Designer will search for a BIDS .json sidecar for the phase encoding information of the input data. If this .json file does not exist, users must use this option to inform Designer of the phase encoding direction. 
- Can be a signed axis number (e.g. -0, 1, +2), an axis designator (e.g. RL, PA, IS), or NIfTI axis codes (e.g. i-, j, k).

### `-echo_time <TE1,TE2,... (milliseconds)>`
- Specify the echo time used in the acquisition (comma separated list the same length as number of inputs), i units of ms.
- By default Designer will search for a BIDS .json sidecar for the echo-time information of the input data. If this .json file does not exist, users must use this option to inform Designer of the TE. 
- Use in conjunction with `eddy`.

### `-bshape <beta1,beta2,...>`
- Specify the b-shape used in the acquisition (comma separated list the same length as number of inputs)'.
- B-shape can be 1 for linear encoding data, 0 for spherical encoding, 0.5 for planar encoding, -0.5 for zeppelin encoding.
- Use in conjunction with `eddy`.

### `-pre_align`
- Rigidly align each input series to correct for large motion between series.
- Potentially useful in cases where there are multiple input series with large motion difference between them.

### `-ants_motion_correction`
- Perform rigid motion correction using ANTs (useful for cases where eddy breaks down).
- Runs volume-to-volume registration for each DWI image using a mutual information cost function and b-spline interpolation.

### `-eddy_quad_output <path/to/eddy_quad/output/folder>`
- `eddy_quad` (QC report) will run after `-eddy`
- Path to a not yet existing folder you want to save eddy_quad output to.
- By default, it will save to the scratch (use `-nocleanup` to save scratch directory) directory (path/to/scratch/eddy_processing/dwi_post_eddy.qc)

{: .ref }
> Matteo Bastiani, Michiel Cottaar, Sean P. Fitzgibbon, Sana Suri, Fidel Alfaro-Almagro, Stamatios N.
Sotiropoulos, Saad Joabdi and Jesper L.R. Andersson. (2019). Automated quality control for within and between studies diffusion MRI data using a non-parametric framework for movement and distortion correction. Neurolmage 184:801-812.

---

## Options relating to Gibbs artifact correction

### `-degibbs`
- Perform (RPG) Gibbs artifact correction. 
- Users may include PF factor with `-pf` (e.g. 6/8, 7/8) and PF dimension `-pe_dir` (1, 2 or 3 for i, j  or k respectively) for data acquired with Partial Fourier sampling.
- For fully sampled input data users may leave out `-pf` and `pe_dir` options.
- When using the `-degibbs` option please use the following citations:

{: .ref }
> Kellner, E., Dhital, B., Kiselev, V. G., & Reisert, M. (2016). Gibbs‐ringing artifact removal based on local subvoxel‐shifts. Magnetic resonance in medicine, 76(5), 1574-1581.
>
> Lee, H. H., Novikov, D. S., & Fieremans, E. (2021). Removal of partial Fourier‐induced Gibbs (RPG) ringing artifacts in MRI. Magnetic resonance in medicine, 86(5), 2733-2750.

### `-pf <partial_fourier_factor>`
- Specify the partial fourier factor either as a fraction or decimal value (e.g. 7/8, 6/8).
- By default Designer will search for a BIDS .json sidecar for the PF information of the input data. If this .json file does not exist, users must use this option to inform Designer of the PF level.

---

## Options related to bias, normalization, and inhomogeneity 

### `-rician`
- Perform Rician bias correction using the approximation for the method of moments.
- Background on this option can be found [here]({{ site.baseurl }}{% link docs/designer/background.md %}#rician-bias-correction).
- If used in conjunction with `-denoise`, the noisemap from MPPCA is used to compute the bias free maginutude DWI dataset.

{: .warning }
The `-rician` method of correcting for bias is specific to data with SNR > 2. If the SNR of your data is lower than 2, using this option can have detrimental effects. Please use this option with care.

- If `-denoise` is not used, $\sigma$ is estimated using MPPCA, but the denoising is not applied. 
- When using the `-rician` option please use the following citations:

{: .ref }
> Koay CG, Basser PJ. Analytically exact correction scheme for signal extraction from noisy magnitude MR signals. J Magn Reson 2006; 179: 317– 322.

### `-normalize`
- Normalize the DWI data such that all b=0 images have the same mean intensity.
- Normalization is performed voxelwise by smoothing each b=0 image with a Gaussian kernel with $\sigma=3$, then rescaling each input series based on the ratio between its smoothed b=0 image and the mean smoothed b=0 image.
- Particularly important for multiple input series, where there may be changes in scanner gain.

### `-b1correct`
- Include a bias field correction step in dwi preprocessing using ANTs N4 bias field correction.
- Designer uses the `dwibiascorrect` program in mrtrix. The default parameters are -b [ 100, 3 ] -c [ 1000, 0.0 ] -s 4. This fits a relatively smooth bias field compared to the N4 defaults for structural images.
- When using the `-b1correct` option please use the following citations:

{: .ref }
> Tustison, N. J., Avants, B. B., Cook, P. A., Zheng, Y., Egan, A., Yushkevich, P. A., & Gee, J. C. (2010). N4ITK: improved N3 bias correction. IEEE transactions on medical imaging, 29(6), 1310-1320.

---

## Other options for Designer

### `-mask`
- Masks the final output image. Has no impact on intermediate masks used to optimize the performance of topup and eddy.
- Uses FSL bet with `-f .2`.
- Saves the output mask in the scratch directory as `brain_mask.nii`

### `-datatype <dtype>`
- Specify the output datatype. Valid options are float32, float32le, float32be, float64, float64le, float64be, int64, uint64, int64le, uint64le, int64be, uint64be, int32, uint32, int32le, uint32le, int32be, uint32be, int16, uint16, int16le, uint16le, int16be, uint16be, cfloat32, cfloat32le, cfloat32be, cfloat64, cfloat64le, cfloat64be, int8, uint8, bit.

### `-fslbvec <bvec1,bvec2,...>`
- Specify bvec path if path is different from the path to the dwi(s) or the file has an unusual extension.

### `-fslbval <bval1,bval2,...>`
- Specify bval path if path is different from the path to the dwi(s) or the file has an unusual extension.

### `-bids <bids1,bids2,...>`
- Specify BIDS .json path(s) if different from the path to the dwi or the file has an unusual extension.

### `-n_cores <ncpus>`
- Specify the number of cores to use in parallel tasks, by default Designer will use total available cores - 2.
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







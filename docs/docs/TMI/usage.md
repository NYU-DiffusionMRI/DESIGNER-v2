---
layout: default
title: usage
has_children: false
nav_order: 6
parent: TMI
---

# Options and Usage for TMI
{: .no_toc }

This page describes each optional argument to TMI. For specific use-case examples of TMI calls, please see the [examples]({% link docs/TMI/examples.md %}) page.

See [Rician bias]({% link docs/TMI/Rician-bias.md %}) for the rationale, all estimator alternatives, noise assumptions, and compatibility with AKC and b0 restoration.

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
- Perform DTI fitting using weighted linear least squares by default, or unconstrained Rician-mean NLLS with `-RicianExp`.
- Include diffusion tensor parameters in output folder.
- Parameters extracted from the diffusion tensor are labelled as `*_dti.nii` (e.g. ad_dti.nii, fa_dti.nii).
- When the DTI option is used, b-shells greater than 1000 are automatically excluded from tensor estimation.

### `-DKI`
- Perform DKI fitting using weighted linear least squares by default, or unconstrained Rician-mean NLLS with `-RicianExp`.
- Include diffusion kurtosis parameters in output folder.
- Parameters extracted from the diffusion tensor are labelled as `*_dki.nii` (e.g. ad_dki.nii, fa_dki.nii, mk_dki.nii).
- When DKI option is used, b-shells greater than 3000 are automatically excluded from tensor estimation.

### `-WDKI`
- Output folder will also contain kurtosis parameters that are computed using definitions with the kurtosis tensor, unlike the conventional definitions which compute them using averages of the directional kurtosis. For example, for mean kurtosis mk_wdki is simply the angular average of the dimensionless cumulant tensor $\overline{W({\bf \hat{n}})} =  \tfrac{1}{4 \pi}  \int_{\mathbb{S}^2} d\hat{\mathbf{n}}\, W(\hat{\mathbf{n}})$ while the more standard definition (mk_dki) is the angular average of the directional kurtosis $\overline{K({\bf \hat{n}})} =  \tfrac{1}{4 \pi} \int_{\mathbb{S}^2} d\hat{\mathbf{n}}\, K(\hat{\mathbf{n}})=  \tfrac{1}{4 \pi}  \int_{\mathbb{S}^2} d\hat{\mathbf{n}}\,  \frac{\bar{D}^2 W(\hat{\mathbf{n}})}{ {D(\hat{\mathbf{n}})}^2}$.
- Parameters extracted from the W tensor are labelled as `*_wdki.nii` (e.g. ad_wdki.nii, fa_wdki.nii, mk_wdki.nii).
- Note that the estimated tensors (D and W) are identical but how we combine the tensor elements to make kurtosis scalar maps differs (DTI maps with _wdki are identical to those with _dki since there is no conflict with DTI-derived metrics).

### `-gnc <field>`

Use a voxelwise gradient field on the DWI grid for DTI/DKI fitting. Pass the original graddev NIfTI. Set `-gnc_format` to specify its source. See the [TORTOISE and gradunwarp examples]({% link docs/TMI/DTI-DKI.md %}#gradient-nonlinearity-correction).

### `-gnc_format <convention>`

Choose `tortoise` for a TORTOISE graddev, `gradunwarp` for a full Jacobian in FSL axes, or `scanner` for a coil tensor already expressed in scanner axes. An explicit choice overrides the header. If omitted, TMI uses the DESIGNER backend name in the NIfTI header’s `descrip` field when recognized. For an empty or unrecognized description, a negative-determinant grid defaults to `gradunwarp`, while a positive-determinant grid requires an explicit choice. The determinant is taken from the 3 × 3 spatial part of the NIfTI affine. MIF fields default to `scanner`.

### `-fit_constraints <0,0,0>`
- For use with `-DKI` or `-WDKI`
- A comma separated list of either 0 or 1. Each index corresponds to a hard bound on the kurtosis tensor. Index 1 corresponds to diffusion tensor > 0, index 2 corresponds to kurtosis tensor > 0, index 3 corresponds to kurtosis tensor < 10. An input of `1,1,1` will perform a fully constrained DKI fit, and `0,0,0` will perform an unconstrained fit.
- Default `0,0,0`
- Supported with `-RicianClipMax` or `-RicianClipAbs`; not yet supported with `-RicianExp` (not applied; a warning is issued).

## Options for outlier replacement

### `-akc_outliers`
- Brute force K tensor outlier detection and filtering. The kurtosis tensor is projected onto a sphere with 10,000 directions. Voxels where mean kurtosis is less than -1 or greater than 10 are labelled as outliers.
- Outliers are replaced by the median value of neighboring (26 adjacent) non-outlier voxels in all diffusion weighted images and the DKI fit is done again.
- Supported with `-RicianClipMax` or `-RicianClipAbs`; not yet supported with `-RicianExp` (not applied; a warning is issued).

### `-fit_smoothing <percentile>`
- Windowed adaptive nonlocal means filter on the dwi.
- Diffusion weighted images undergo an NLM filter over the most similar neighboring pixels and directions prior to fitting.
- For each voxel we take a cubic patch in 4d and take the mean over the $n$ most similar voxels in *xyz* and *q-space*. n is decided by the input argument to tmi, where with n=10 we are smoothing over the 10% most similar voxels, and so on.
- If used in conjunctions with `-akc_outliers` outliers are detected and ignored during the smoothing process.
- The `<percentile>` Defines the degree of similar voxels in a sliding window patch to average. A percentile of 10 will average the 10% most similar voxels and directions in the sliding window. We recommend a percentile of 10 as a starting point for tuning this option.

## Options for model estimation

### `-SH`

- Estimate signal rotational invariants without SMI model fitting or training priors. Use either `-SH` or `-SMI`.
- **`-SH` alone uses ordinary linear least squares (LLS)** with automatic per-shell fitting orders and no signal nonnegativity constraint. No noise map is required.
- **`-SH -RicianClipMax -sigma <noisemap>` fits max-zero corrected magnitudes with ordinary LLS.**
- **`-SH -RicianClipAbs -sigma <noisemap>` fits absolute-value corrected magnitudes with ordinary LLS.**
- **`-SH -RicianExp -sigma <noisemap>` fits the expected Rician magnitude signal without constraints by default. Add `-nonNeg` for fitted-signal nonnegativity.** Supplying `-sigma` alone does not change the default LLS estimator.
- Save separate `S0.nii`, `S2.nii`, ..., `S<maximum_order>.nii` files in original signal units. Each file is 4D, with one volume per acquisition shell, including a singleton volume dimension for a single shell.
- `SH_shells.json` maps each zero-based volume index to its b-value (s/mm²), b-shape, echo time (ms), number of measurements, and fitted SH order. Shells are ordered by b-value, b-shape, then echo time.
- `-lmaxSH 6` fits through order 6 and saves S0, S2, S4, and S6. Without `-lmaxSH`, use the existing automatic per-shell orders and save every fitted even order. Higher-order volumes are zero for shells fitted at a lower order; b=0 and isotropic shells use scalar estimation with the selected correction and zero higher orders.
- A comma-separated list sets each shell's fitting order independently for all estimators, for example `-lmaxSH 0,2,4,4,6` for b=0,1000,2000,7000,25000 at one b-shape and TE. Include every shell, sorted by b-value, b-shape, then echo time; `SH_shells.json` records the mapping. See the [per-shell example]({% link docs/TMI/SMI.md %}#choose-a-fitting-order-for-each-shell).
- A single diffusion shell is supported, including acquisitions without b=0. A requested fitting order must have enough independent directions in every directional shell; otherwise reduce `-lmaxSH`.
- All estimators retain the existing SH normalization and rank-one processing across shells; no microstructure regression is performed.

### `-SMI`

- Perform estimation of the Standard Model of White Matter.
- SH fitting uses LLS by default. Choose `-RicianClipMax -sigma <noisemap>` or `-RicianClipAbs -sigma <noisemap>` for corrected LLS. Use `-RicianExp -sigma <noisemap>` for expected Rician magnitude fitting, with optional `-nonNeg`, exactly as in standalone `-SH`.
- `-lmaxSH` sets standalone SH fitting orders and is supported only with `-SH`. With `-SMI`, a single `-lmax` value selects the invariant order for model regression while SH fitting orders remain automatic per shell.
- Use in conjunction with `-compartments`, `-sigma`, `-bshape`, and `-echo_time` options.
- We recommend prior estimation of the noise level sigma using `designer` along with the `-denoise` option.

### `-RicianExp`

- Fit the expected Rician magnitude to original observations for DTI/DKI/WDKI/SH/SMI. Requires `-sigma`.
- Recommended for SH/SMI data including $b\geq3000\,\mathrm{s/mm^2}$. Use `-RicianClipAbs` when all shells are below this threshold, and for DTI/DKI/WDKI. The selected estimator applies to every shell in the call; there is no automatic switching by b-value.
- Slower than `-RicianClipAbs` because it uses nonlinear optimization. For typical SH/SMI datasets entirely below $3000\,\mathrm{s/mm^2}$, the additional benefit generally does not justify the extra runtime.
- Initialization: ordinary LLS on absolute-value corrected measurements for tensors, max-zero corrected measurements for SH/SMI. Negative observations after denoising are preserved.
- **Tensors:** unconstrained signal-domain NLLS; `-fit_constraints` and `-akc_outliers` are not yet supported (not applied; a warning is issued). `-db0restore`, `-fit_smoothing`, `-polyreg` and external `-dt` are not yet supported.
- **SH/SMI:** unconstrained by default; add `-nonNeg` for signal nonnegativity. b0/isotropic shells invert the Rician mean using original sigma; rank-one processing is unchanged.
- Overrides either clipping flag for all selected models, without changing these initialization choices.
- See [Rician bias]({% link docs/TMI/Rician-bias.md %}) for noise assumptions and [tensor details]({% link docs/TMI/DTI-DKI.md %}#rician-correction-options) for diagnostics.

### `-RicianClipMax`

Correct measurements using

$$\widetilde M=\sqrt{\max(0,M^2-\sigma^2)}.$$

Negative squared differences become zero.

### `-RicianClipAbs`

Recommended Rician correction for DTI/DKI/WDKI, followed by WLLS.
Also recommended for SH/SMI when all shells are below $3000\,\mathrm{s/mm^2}$,
followed by ordinary SH LLS.

Correct measurements using

$$\widetilde M=\sqrt{\lvert M^2-\sigma^2\rvert}.$$

Negative squared differences are reflected upward. This avoids most corrected
zeros, but remains zero at $M=\sigma$ and is nonmonotonic below it.

For both clipping flags:

- Requires `-sigma` and at least one of `-DTI`, `-DKI`, `-WDKI`, `-SH` or `-SMI`. Choose one clipping flag; selecting both is an error.
- **Tensors:** WLLS on corrected measurements. Correction is applied once, including b0; constraints, AKC, smoothing, polynomial regression, GNC and b0 restoration are supported. Corrected zeros use the existing positive log floor.
- **SH/SMI:** ordinary corrected LLS, without signal constraints. b0/isotropic shells average individually corrected measurements; rank-one processing and SMI regression are unchanged.
- Supply original data with a retained noise floor and an aligned original sigma map. Both formulas are approximate corrections.
- A combined run uses the selected formula for every model; use separate calls for different choices. `-RicianExp` takes precedence.

### `-nonNeg`

- Requires `-RicianExp` with `-SH` or `-SMI`. Invalid with clipping alone or tensor-only fitting.
- Enforces a nonnegative SH signal at acquired directions and on a 65,536-direction sphere grid, before rank-one processing. It does not constrain tensor fits.
- Optional: a positive physical signal can have a negative truncated SH projection, so the constraint can introduce truncation bias.

### `-RicianSH_maxiter N`

Requires `-RicianExp -nonNeg` with `-SH` or `-SMI`. Caps constrained refinement
and saves [SH diagnostics]({% link docs/TMI/SMI.md %}#iteration-limits-and-diagnostics).
This cap does not apply to tensor fits.

### `-WMTI`
- Include WMTI parameters in output folder (awf,ias_params,eas_params).
- WMTI estimation is performed using the `dipy` package.

### `-compartments <IAS,EAS,FW>`
- SMI compartments (IAS, EAS, and FW) expressed as a comma separated string.
- default: `IAS,EAS`

### `-select_prior tag.mat` (choose from 4 predetermined training data files)
- SMI training data options: `SMI_Gaussian_noFWPrior.mat`, `SMI_Gaussian_wFWPrior.mat`, `SMI_Uniform_noFWPrior.mat`, or `SMI_Uniform_wFWPrior.mat`.
- default: `Gaussian` (that matches the choice in `-compartments`)

### `-load_prior filename.mat` (load a user defined training data file)
- SMI training data should be [ Ntraining x 9], with each column being a model parameter: [f, Da, Depar, Deperp, fw, T2a, T2e, p2, p4]

### `-sigma <noisemap>`
- Path to noise map. Required with `-RicianExp`, `-RicianClipMax` or `-RicianClipAbs`, including DTI/DKI/WDKI correction; optional but recommended for SMI. Standalone `-SH` LLS does not need a noise map.
- We recommend computing sigma prior to running `tmi` using `designer`.

### `-lmaxSH <orders>`

- Use with `-SH` to set the maximum SH fitting and output order for any estimator. A single nonnegative even integer applies to every directional shell: `-lmaxSH 6` saves S0, S2, S4, and S6; `-lmaxSH 0` saves only S0. If omitted, fitting orders are selected automatically per shell.
- A comma-separated list specifies each shell's order independently, sorted by b-value, b-shape, then echo time, including b=0/isotropic shells. For example, `-lmaxSH 0,2,4,4,6` assigns orders 0, 2, 4, 4, and 6 to b=0,1000,2000,7000,25000 when these are the input shells at one b-shape and TE. Use nonnegative even integers, no spaces or trailing comma, and exactly one entry per shell.
- Higher-order volumes for shells fitted at a lower order are zero placeholders, not measurements. `SH_shells.json` records each shell's fitted order.
- This option is not accepted with `-SMI`. When adapting an earlier standalone SH command, replace `-lmax` with `-lmaxSH`; `-SH -lmax` produces a migration message.

### `-lmax <lmax>`

- With `-SMI`: a single invariant order used for model regression (0, 2, 4, or 6). Per-shell lists are not accepted; SH fitting orders remain automatic per shell, for all three SH estimators.
- For standalone SH fitting orders, use `-SH -lmaxSH` instead.

### `-bshape <beta1,beta2,...>`
- Specify the b-shape used in the acquisition (comma separated list the same length as number of inputs).

### `-echo_time <TE1,TE2,...>`
- Specify the echo time used in the acquisition (comma separated list the same length as number of inputs).

### `-maxb <b-value>`
- Specify the maxmimum b-value to use during DKI estimation. 
- Units should be in ms/µm^2.
- default maximum b-value for DKI estimation is 3.0 ms/µm^2

### `-maxb_dti <b-value>`
- Specify the maxmimum b-value to use during DTI estimation. 
- Units should be in ms/µm^2.
- default maximum b-value for DTI estimation is 1.5 ms/µm^2

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







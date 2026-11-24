---
layout: default
title: examples
has_children: false
nav_order: 7
parent: TMI
---

# Example Usage for TMI
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---


## Basic Example

A basic call to `tmi`, performing single shell DTI estimation and multi-shell DKI estimation will look like:
```
tmi -DTI -DKI dwi.mif parameters
```
In this case `tmi` assumes that gradient information has been stored in the `.mif` header and this information will be used as inputs to fitting functions. All output parameter maps are saved as nifti files and stored in the directory `parameters`.

In general, we recommend running parameter estimation within a tissue mask, to decrease the likelihood of fit degeneracies and unnecessary memory overhead. A mask can be included in the `tmi` call using the `-mask` option.
```
tmi -DTI -DKI -mask /path/to/brain_mask.nii dwi.mif parameters
```

---

## Extract the W tensor 
```
tmi -DKI -WDKI dwi.mif parameters
```
In this case the W form of the kurtosis tensor is returned in addition to the K form. Differing sets of tensor parameters will be saved in the `parameters` directory with the suffixes `*_dki.nii` and `*_wdki.nii`.

---

## Outlier replacement and filtering
This call to `tmi` will perform  outlier replacement along with nonlocal means smoothing (user may adjust how much smoothing by adjusting `-fit_smoothing <percentile>`) as a filtering step prior to fitting the kurtosis tensor.
```
tmi -WDKI -DKI -akc_outliers -fit_smoothing 25 dwi.mif parameters
```

{: .note }
It is important to note that when performing outlier correction that constrained fitting is *not* used. Constrained fitting will bias the output parameters and potentially hide voxels that should be labelled as outliers.

---

## constrained fitting
Constraints can be used to bound the kurtosis tensor estimation, an example of running a kurtosis fit with a positivity constraint on kurtosis looks like:
```
tmi -WDKI -DKI -fit_constraints 0,1,0 dwi.mif parameters
```
- The order of values input to the `-fit_constraints` option is: diffusion tensor > 0, kurtosis tensor > 0, and kurtosis tensor < 10. These settings are supported by WLLS, including with either clipping flag; they are not yet supported with `-RicianExp`.

--- 

## SMI
Inputting data to SMI can be slightly more complex, particularly if the data comes from a multi echo-time or multi bshape acquisition. An example is shown here:
```
dwi1=dwi_te92_LTE.nii
dwi2=dwi_te92_ZTE.nii
dwi3=dwi_te92_STE.nii
dwi4=dwi_te62_LTE.nii
dwi5=dwi_te78_ZTE.nii
dwi6=dwi_te130_LTE.nii

tmi \
-SMI \
-sigma /path/to/designer_processing/sigma.nii \
-compartments EAS,IAS,FW \
-echo_time 92,92,92,62,78,130 \
-bshape 1,-0.5,0,1,-0.5,1 \
-mask /path/to/designer_processing/brain_mask.nii \
-scratch tmi_processing_variable_te_beta -nocleanup \
$dwi1,$dwi2,$dwi3,$dwi4,$dwi5,$dwi6 parameters
```
This command will run SMI on data with multiple echo times and multiple bshapes, it will include all cellular compartments and use the noisemap computed by a prior `designer` call.

If the input to SMI is purely LTE, the call to SMI can be simplified:
```
tmi \
-SMI \
-sigma /path/to/designer_processing/sigma.nii \
-compartments EAS,IAS \
-mask /path/to/designer_processing/brain_mask.nii \
-scratch tmi_processing_lte -nocleanup \
dwi.mif parameters
```

## Docker example

`docker run -it` will take you inside a container. You can add `-v` to attach your data to the container for processing. Use `-v host_path:container_path` to mount your folder (host_path) to a specified path inside the container (container_path). 

The below example mounts the host_path (`/path/to/folder/with/dataset`) to the container_path (`/data`). The host_path contains the input series `dwi.mif`. The following is an example that calls tmi and will extract WDKI and DKI maps with outlier correction and nonlocal means smoothing:
```
docker run -it -v /path/to/folder/with/dataset:/data \
nyudiffusionmri/designer2:<tag> tmi \
-WDKI -DKI -akc_outliers -fit_smoothing 25 /data/dwi.mif /data/parameters
```

---

## Singularity example

`singularity run` will run a Singularity container. `--bind` allows you to attach your data to the container for processing. Use `--bind host_path:container_path` to mount your folder (host_path) to a specified path inside the container (container_path). 

The below example mounts the host_path (`/path/to/folder/with/dataset`) to the container_path (`/mnt`). The host_path contains the input series `dwi.mif`. The following is an example that calls tmi and will extract WDKI and DKI maps with outlier correction and nonlocal means smoothing:
```
singularity run --bind /path/to/folder/with/dataset:/mnt \
designer2_<tag>.sif tmi \
-WDKI -DKI -akc_outliers -fit_smoothing 25 /mnt/dwi.mif /mnt/parameters
```

---


### Recommended Rician corrections

Use the original DWI with a retained Rician noise floor and an aligned original
noise map. Add `-mask brain_mask.nii` to restrict fitting to the brain.
For SH/SMI, use `-RicianExp` when the data include
$b\geq3000\,\mathrm{s/mm^2}$; use `-RicianClipAbs` when all shells are below
this threshold. The choice applies to the whole fit, including mixed-shell data.
Expectation fitting is slower; for typical datasets below this threshold,
its additional benefit generally does not justify the extra runtime.

```bash
# DTI + DKI with absolute-value correction, followed by WLLS
tmi dwi.nii tensor_abs -DTI -DKI -RicianClipAbs -sigma sigma.nii

# DTI + DKI with absolute-value correction and b0 restoration
tmi dwi.nii tensor_abs_restore -DTI -DKI -RicianClipAbs -sigma sigma.nii -db0restore

# SH/SMI with all shells below 3000 s/mm^2 (choose the desired output)
tmi dwi.nii sh_abs -SH -RicianClipAbs -sigma sigma.nii
tmi dwi.nii smi_abs -SMI -RicianClipAbs -sigma sigma.nii

# SH invariants only, using Rician expectation
tmi dwi.nii sh_exp -SH -RicianExp -sigma sigma.nii

# SMI using Rician expectation
tmi dwi.nii smi_exp -SMI -RicianExp -sigma sigma.nii

# SMI + separate SH maps with nonnegative fitted SH signals (two calls)
tmi dwi.nii smi_positive -SMI -RicianExp -sigma sigma.nii -nonNeg
tmi dwi.nii sh_positive -SH -RicianExp -sigma sigma.nii -nonNeg
```

`-SMI` already saves invariants in `rotinvs_smi.nii`. The extra `-SH` call in
the last workflow saves separate `S0.nii`, `S2.nii`, … files. Choose either
`-SMI` or `-SH` per call; the flags cannot be combined.

Tensor clipping uses WLLS. Both clipping choices support b0 restoration and
AKC replacement; add `-akc_outliers` when needed. SH/SMI expectation fits are
unconstrained unless `-nonNeg` is supplied. Use separate calls for tensor and
SH/SMI correction, because `-RicianExp` overrides either clipping flag for
all selected models. See [Rician bias]({% link docs/TMI/Rician-bias.md %}) for
alternative estimators and compatibility.

### Standalone SH fitting orders

For data including $b\geq3000\,\mathrm{s/mm^2}$, fit through SH order 6 with
Rician expectation:

```bash
tmi dwi.nii sh_exp_l6 -SH -RicianExp -lmaxSH 6 \
  -sigma sigma.nii -mask brain_mask.nii
```

This writes `S0.nii`, `S2.nii`, `S4.nii`, and `S6.nii`, with one volume per
shell. `SH_shells.json` records the b-values, b-shapes, echo times,
measurement counts, and fitted orders.

Omit `-lmaxSH` for automatic per-shell orders, or use a comma-separated list
including b0, such as `-lmaxSH 0,2,4`. Every requested order must be supported
by the acquired directions. Use `-lmaxSH 0` for S0 only. For multi-TE or
multi-b-shape data, supply the same acquisition metadata as for SMI.

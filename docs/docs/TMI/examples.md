---
layout: default
title: examples
has_children: false
nav_order: 6
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
- The order of values input to the `-fit_constraints` option is: diffusion tensor > 0, kurtosis tensor > 10, and kurtosis tensor < 10.

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


### Rician-aware SH estimation

For magnitude data with an appropriate noise map and a retained Rician noise floor:

```bash
tmi dwi.nii smi_rician \
  -SMI -RicianSH \
  -sigma sigma.nii -mask brain_mask.nii
```

Omit `-RicianSH` to use the existing linear SH fit:

```bash
tmi dwi.nii smi_linear -SMI -sigma sigma.nii -mask brain_mask.nii
```

Both commands produce the same named SMI parameter and rotational-invariant maps.
The flag changes only SMI's SH estimation, including the multi-TE/b-shape path;
it does not change DTI/DKI fits requested in the same command.

### Rotational invariants without SMI

Fit through SH order 6 using the Rician expected magnitude and nonnegative
signal constraints, without training or fitting the Standard Model:

```bash
tmi dwi.nii sh_rician -SH -RicianSH -lmaxSH 6 \
  -sigma sigma.nii -mask brain_mask.nii
```

This writes `S0.nii`, `S2.nii`, `S4.nii`, and `S6.nii`. Each image contains one
volume per shell; `SH_shells.json` provides the corresponding b-values,
b-shapes, echo times, measurement counts, and fitted orders. No SMI parameter
maps or combined `rotinvs_smi.nii` are generated.

For the linear fit, no noise map is required:

```bash
tmi dwi.nii sh_linear -SH -lmaxSH 6 -mask brain_mask.nii
```

Omit `-lmaxSH` to use automatic per-shell SH orders, or use `-lmaxSH 0` to save
only S0. The requested order must be supported by the directions in each
diffusion shell. A single shell without b=0 is supported. For multi-TE or
multi-b-shape data, use the same acquisition metadata options as for SMI.
Choose either `-SH` or `-SMI` in a command.

### Clipped-corrected SH fitting with signal nonnegativity

```bash
tmi dwi.nii sh_clip -SH -RicianClip -sigma sigma.nii -mask brain_mask.nii
tmi dwi.nii smi_clip -SMI -RicianClip -sigma sigma.nii -mask brain_mask.nii
```

The fitting target is `sqrt(max(0, M²−sigma²))`, with nonnegativity enforced
on the underlying SH signal. Supply the original magnitude DWI.
`-RicianSH` takes precedence if both correction flags are present. Per-shell orders in standalone mode
are set with `-lmaxSH`, as for the other estimators.

### Clipped tensor inputs with an independent SH estimator

```bash
# Clipped correction before the existing DTI and DKI fits
tmi dwi.nii tensor_clip -DTI -DKI -RicianClip -sigma sigma.nii

# Clipped correction for tensors; Rician expectation + nonnegativity for SH
tmi dwi.nii combined -DTI -DKI -SH -RicianSH -sigma sigma.nii
```

`-RicianSH` alone enables clipped correction for tensors and keeps original
magnitudes for the SH expectation fit. Adding `-RicianClip` is unnecessary
and does not correct the data twice. With tensor fitting alone, the two flags
are equivalent. Replace `-SH` with `-SMI` to include Standard Model fitting.

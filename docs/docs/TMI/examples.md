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
tmi -DKI -akc_outliers -fit_smoothing 25 dwi.mif parameters
```

{: .note }
It is important to note that when performing outlier correction that constrained fitting is *not* used. Constrained fitting will bias the output parameters and potentially hide voxels that should be labelled as outliers.

---

## constrained fitting
Constraints can be used to bound the kurtosis tensor estimation, an example of running a kurtosis fit with a positivity constraint on kurtosis looks like:
```
tmi -DKI -fit_constraints 0,1,0 dwi.mif parameters
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

The below example mounts the host_path (`/path/to/folder/with/dataset`) to the container_path (`/data`). The host_path contains the input series `dwi.mif`. The following is an example that calls tmi and will extract DKI maps with outlier correction and nonlocal means smoothing:
```
docker run -it -v /path/to/folder/with/dataset:/data \
nyudiffusionmri/designer2:<tag> tmi \
-DKI -akc_outliers -fit_smoothing 25 /data/dwi.mif /data/parameters
```

---

## Singularity example

`singularity run` will run a Singularity container. `--bind` allows you to attach your data to the container for processing. Use `--bind host_path:container_path` to mount your folder (host_path) to a specified path inside the container (container_path). 

The below example mounts the host_path (`/path/to/folder/with/dataset`) to the container_path (`/mnt`). The host_path contains the input series `dwi.mif`. The following is an example that calls tmi and will extract DKI maps with outlier correction and nonlocal means smoothing:
```
singularity run --bind /path/to/folder/with/dataset:/mnt \
designer2_<tag>.sif tmi \
-DKI -akc_outliers -fit_smoothing 25 /mnt/dwi.mif /mnt/parameters
```

---

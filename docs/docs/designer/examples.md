---
layout: default
title: examples
has_children: false
nav_order: 4
parent: designer
---

# Example Usage for Designer
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Basic Example

A basic call to designer, applying MPPCA denoising to a single input series called `dwi.nii` might look like:
```
designer -denoise dwi.nii dwi_designer.nii
```
While this call does not require gradient information in order to run, some options such as `-eddy` do. Therefore we assume that there are gradient files called `dwi.bvec` and `dwi.bval` located in the same directory as `dwi.nii`.

Users may also include gradient directions into the header of the image file by converting to `.mif` format prior to running designer:
```
mrconvert -fslgrad dwi.bvec dwi.bval dwi.nii dwi.mif
designer -denoise dwi.mif dwi_designer.mif
```

Using `.mif` format is helpful overall as gradient information is automatically rotated along with image volumes during registration steps.

---

## Recommended: Denoising, RPG degibbsing, eddy current and motion correction (default [DESIGNER pipeline](https://doi.org/10.1162/imag_a_00125))

Given a partial Fourier (6/8) sampled magnitude dataset with two shells: dwi_b1000.nii and dwi_b2000.nii (along with accompanying .bvec and .bval files), acquired in the AP direction, along with a reverse phase encoding b=0 image called "rpe_b0.nii", the call to designer to perform complete preprocessing would look like:
```
designer -denoise -shrinkage frob -adaptive_patch -rician \
-degibbs -pf 0.75 -pe_dir j- \
-eddy -rpe_pair rpe_b0.nii \
-normalize -mask \
-scratch designer_processing -nocleanup \
dwi_b1000.nii,dwi_b2000.nii dwi_designer.nii
```

---

## Denoising Magnitude and Phase data

A more complex call to designer, applying MPPCA denoising on complex valued data (including phase images called `phase.nii`) with custom optional arguments to perform denoising using the `jespersen` PCA cutoff and singular value shrinkage looks like:
```
designer -denoise -shrinkage frob -algorithm jespersen -phase phase.nii dwi.nii dwi_designer.nii
```

---

## Processing data with variable b-tensor shapes (multiple b-shapes)

A more complex call to designer is required to process data with variable b-shapes. Given 2 input series with severe motion artifacts, one acquired with single diffusion encoding (SDE) and another one with double diffusion encoding (DDE), an example of such a call might look like:
```
dwi1=SDE_b1000_b2000_TE106.nii.gz
dwi2=DDE_b2000_TE106.nii.gz

designer \
-denoise -degibbs \
-pre_align -ants_motion_correction \
-eddy -rpe_pair rpe_b0.nii -rpe_te 106 -pe_dir AP \
-normalize -mask \
-scratch designer_variable_beta -nocleanup \
-eddy_fakeb 1,1.25 -bshape 1,-0.5 -echo_time 106,106 \
$dwi1,$dwi2 dwi_designer.mif
```
In this case, eddy is run for the full data, making sure different b, and b-shape combinations are not mixed. 


## Processing data with variable echo time and b-shape

A more complex call to designer is required to process data with variable echo time and b-shape. Given 6 input series with severe motion artifacts and with varying echo times and b-shapes, an example of such a call might look like:
```
dwi1=dwi_te92_LTE.nii
dwi2=dwi_te92_ZTE.nii
dwi3=dwi_te92_STE.nii
dwi4=dwi_te62_LTE.nii
dwi5=dwi_te78_ZTE.nii
dwi6=dwi_te130_LTE.nii

designer \
-denoise -degibbs \
-pre_align -ants_motion_correction \
-eddy -rpe_pair rpe_b0.nii -rpe_te 62 -pe_dir AP \
-normalize -mask \
-scratch designer_variable_te_beta -nocleanup \
-eddy_fakeb 1,1.4,0.3,0.75,1.7,1.25 -bshape 1,0.6,0,1,0.6,1 -echo_time 92,92,92,62,80,130 \
$dwi1,$dwi2,$dwi3,$dwi4,$dwi5,$dwi6 dwi_designer.mif
```
In this case, eddy is run for the full data, making sure different b, b-shape, and echo-time combinations are not mixed. 

---

## Docker example

`docker run -it` will take you inside a container. You can add `-v` to attach your data to the container for processing. Use `-v host_path:container_path` to mount your folder (host_path) to a specified path inside the container (container_path). 

The below example mounts the host_path (`/path/to/folder/with/dataset`) to the container_path (`/data`). The host_path contains the input series `dwi.nii` with its bvec (`dwi.bvec`),  bval (`dwi.bvec`), and json (`dwi.json` containing PF factor, PE dir, and TE) along with reverse phase-encoding b=0 image (`rpe_b0.nii`). The following is an example that calls designer and will process the data using the recommended designer pipeline:
```
pa=/data/rpe_b0.nii
dwi=/data/dwi.nii
dwi_out=/data/dwi_designer.nii

docker run -it -v /path/to/folder/with/dataset:/data \
nyudiffusionmri/designer2:<tag> designer \
-denoise -shrinkage frob -adaptive_patch -rician \
-degibbs \
-eddy -rpe_pair $pa \
-normalize -mask \
-scratch /data/processing -nocleanup \
$dwi $dwi_out
```

---

## Singularity example

`singularity run` will run a Singularity container. `--bind` allows you to attach your data to the container for processing. Use `--bind host_path:container_path` to mount your folder (host_path) to a specified path inside the container (container_path). 

The below example mounts the host_path (`/path/to/folder/with/dataset`) to the container_path (`/mnt`). The host_path contains the input series `dwi.nii` with its bvec (`dwi.bvec`),  bval (`dwi.bvec`), and json (`dwi.json` containing PF factor, PE dir, and TE) along with reverse phase-encoding b=0 image (`rpe_b0.nii`). The following is an example that calls designer and will process the data using the recommended designer pipeline:
```
pa=/mnt/rpe_b0.nii
dwi=/mnt/dwi.nii
dwi_out=/mnt/dwi_designer.nii

singularity run --bind /path/to/folder/with/dataset:/mnt \
designer2_<tag>.sif designer \
-denoise -shrinkage frob -adaptive_patch -rician \
-degibbs \
-eddy -rpe_pair $pa \
-normalize -mask \
-scratch /mnt/processing -nocleanup \
$dwi $dwi_out
```

If you're running this on a machine that uses GPU and CUDA, you can include `-nv` to use eddy_cuda. The following is an example that will use eddy_cuda when running on a GPU node on a HPC server:
```
module load singularity/3.9.8
module load cuda/11.8

pa=/mnt/rpe_b0.nii
dwi=/mnt/dwi.nii
dwi_out=/mnt/dwi_designer.nii

singularity run --nv --bind /path/to/folder/with/dataset:/mnt \
designer2_<tag>.sif designer \
-denoise -shrinkage frob -adaptive_patch -rician \
-degibbs \
-eddy -rpe_pair $pa \
-normalize -mask \
-scratch /mnt/processing -nocleanup \
$dwi $dwi_out
```
This is compatible with CUDA versions 9.1 to 11.8. Other versions have not been tested.

---

## General recommended usage

In general, we recommend running designer using the `-nocleanup` option and the `-mask` option turned on. This way, certain files that are useful for subsequent processing are retained. These include the brain mask generated after eddy current correction and the noisemap generated during denoising. For users intending to run TMI to estimate parametric maps after running Designer, these images are useful inputs.

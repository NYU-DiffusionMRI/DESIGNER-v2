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

## Denoising Magnitude and Phase data

A more complex call to designer, applying MPPCA denoising on complex valued data (including phase images called `phase.nii`) with custom optional arguments to perform denoising using the `jespersen` PCA cutoff and singular value shrinkage looks like:
```
designer -denoise -shrinkage frob -algorithm jespersen -phase phase.nii dwi.nii dwi_designer.nii
```

---

## Denoising, RPG degibbsing, eddy current and motion correction

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
-echo_time 92,92,92,62,78,130 \
-bshape 1,-0.5,0,1,-0.5,1 \
-eddy_groups 1,2,3,4,5,6 \
-normalize -mask \
-scratch designer_processing_variable_te_beta -nocleanup \
$dwi1,$dwi2,$dwi3,$dwi4,$dwi5,$dwi6 dwi_designer.mif
```
In this case, eddy is run 6 separate times, for each unique echo-time b-shape combination. 

---

## General recommended usage

In general, we recommend running designer using the `-nocleanup` option and the `-mask` option turned on. This way, certain files that are useful for subsequent processing are retained. These include the brain mask generated after eddy current correction and the noisemap generated during denoising. For users intending to run TMI to estimate parametric maps after running Designer, these images are useful inputs.

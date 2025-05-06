# designer-v2

## Introduction
Designer is a python tool for diffusion MRI preprocessing. It includes:

* denoising using MPPCA (or optionally using patch2self through dipy)
* RPG Gibbs artifact correction
* Rician bias correction
* EPI distortion correction
* Eddy current and motion correction
* b0 normalization for multiple input series
* b1 inhomogeneity correction

For complete documentation on designer & tmi installation and usage please visit the [documentation](https://nyu-diffusionmri.github.io/DESIGNER-v2).

## Installation

designer-v2 is made available as a [python package](https://pypi.org/project/designer2/), and a [docker image](https://hub.docker.com/r/nyudiffusionmri/designer2/) bundling all dependencies. See the [Installation](https://nyu-diffusionmri.github.io/DESIGNER-v2/docs/designer/installation/) documentation for detailed stps


## Example usage for "meso" data

An example script can be found in the examples folder. It is copied here as well. Preprocessing and fitting are now split into two separate functions: designer for preprocessing and tmi for fitting. 

```
#!/bin/bash

datapath=/mnt/labspace/Projects/MESO_v2.0/ALLSUBJS_2.0/M9734
meso1=M9734_074YM_DIFF_meso.nii
meso2=M9734_074YM__DIFF_meso_research.nii
pa=M9734_074YM_DIFF_meso_PA.nii

cd $datapath

designer \
-denoise -algorithm jespersen -extent 7,7,7 \
-degibbs \
-mask \
-scratch designer2_processing_test -nocleanup \
$meso1,$meso2 designer2_test.mif

tmi \
-DTI -DKI -WDKI -SMI \
-mask designer2_processing_test/brain_mask.nii \
-sigma designer2_processing_test/sigma.nii \
-nocleanup \
designer2_test.mif designer2_params_test
```


# References
1. Ades-Aron, B., Veraart, J., Kochunov, P., McGuire, S., Sherman, P., Kellner, E., … & Fieremans, E. (2018). Evaluation of the accuracy and precision of the diffusion parameter EStImation with Gibbs and NoisE removal pipeline. *Neuroimage*, 183, 532-543.

2. Chen, J., Ades-Aron, B., Lee, H.H., Mehrin, S., Pang, M., Novikov, D. S., Veraart, J., & Fieremans, E. (2024). Optimization and validation of the DESIGNER preprocessing pipeline for clinical diffusion MRI in white matter aging. *Imaging Neuroscience*, 2, 1–17.

3. Tournier, J.-D., Smith, R. E., Raffelt, D., Tabbara, R., Dhollander, T., Pietsch, M., Christiaens, D., Jeurissen, B., Yeh, C.-H., & Connelly, A. (2019). MRtrix3: A fast, flexible and open software framework for medical image processing and visualisation. *NeuroImage*, 202, 116137.


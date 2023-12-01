# designer-v2

## For complete documentation on designer & tmi installation and usage please visit the [documentation](https://nyu-diffusionmri.github.io/DESIGNER-v2).

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

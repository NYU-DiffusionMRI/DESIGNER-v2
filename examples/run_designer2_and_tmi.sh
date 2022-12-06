#!/bin/bash
export PYTHONPATH=/data/apps/mrtrix3/bin

datapath=/mnt/labspace/Projects/MESO_v2.0/ALLSUBJS_2.0/M9734
meso1=M9734_074YM_DIFF_meso.nii
meso2=M9734_074YM__DIFF_meso_research.nii
pa=M9734_074YM_DIFF_meso_PA.nii

cd $datapath

# python /cbi05data/data1/Hamster/Ben/designer2/bin/designer \
# -denoise -algorithm jespersen -extent 7,7,7 \
# -degibbs \
# -mask \
# -scratch designer2_processing_test -nocleanup \
# $meso1,$meso2 designer2_test.mif

python /cbi05data/data1/Hamster/Ben/designer_v2_dev/bin/tmi \
-DTIparams -DKIparams -WDKI -SMIparams \
-mask designer2_processing_test/brain_mask.nii \
-sigma designer2_processing_test/sigma.nii \
-nocleanup \
designer2_test.mif designer2_params_test

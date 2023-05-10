#!/bin/bash

#root=/Volumes/Research/fieree01lab
root=/mnt
datapath=$root/labspace/Santiago/PyDesigner/ValidationData/subj_1_22oct21/LTE_ZTE_variableTE
dwi1=20211022_194435sc1DIFF2isoTE92LTEs031a001.nii.gz
dwi2=20211022_194435sc1FWF2isoTE92PTEs005a001.nii.gz
dwi3=20211022_194435sc1FWF2isoTE92STEs007a001.nii.gz
dwi4=20211022_194435sc1FWF2isoTE62LTEs009a001.nii.gz
dwi5=20211022_194435sc1FWF2isoTE78PTEs013a001.nii.gz
dwi6=20211022_194435sc1FWF2isoTE130LTEs011a001.nii.gz
pa=20211022_194435sc1FWF2isoTE62LTEPAs015a1001.nii.gz

cd $datapath

# designer \
# -eddy -rpe_pair $pa -rpe_te 62 -pe_dir AP \
# -echo_time 92,92,92,62,78,130 \
# -bshape 1,-0.5,0,1,-0.5,1 \
# -eddy_groups 1,2,3,4,5,6 \
# -scratch designer2_processing_test -nocleanup \
# $dwi1,$dwi2,$dwi3,$dwi4,$dwi5,$dwi6 designer2_test.mif

designer \
-pre_align -qc -ants_motion_correction \
-eddy -rpe_pair $pa -rpe_te 62 -pe_dir AP \
-echo_time 92,92,92,62,78,130 \
-bshape 1,-0.5,0,1,-0.5,1 \
-eddy_groups 1,2,3,4,5,6 \
-scratch designer2_processing_test -nocleanup \
$dwi1,$dwi2,$dwi3,$dwi4,$dwi5,$dwi6 designer2_test.mif


# tmi \
# -DTI -DKI -WDKI -SMI \
# -mask designer2_processing_test/brain_mask.nii \
# -sigma designer2_processing_test/sigma.nii \
# -nocleanup \
# designer2_test.mif designer2_params_test

#!/bin/bash


datapath=/mnt/labspace/Projects/MESO_v2.0/ALLSUBJS_2.0/M9734
meso1=M9734_074YM_DIFF_meso.nii
meso2=M9734_074YM__DIFF_meso_research.nii
pa=M9734_074YM_DIFF_meso_PA.nii

cd $datapath
#run tests on whether inputs and their parameters are recognized.

# 1 input series
python designer \
-scratch designer2_processing_test -nocleanup \
$meso1 designer1_test.mif

# 2 input series
python designer \
-scratch designer2_processing_test -nocleanup \
$meso1,$meso2 designer1_test.mif

# a dicom series as input
python designer \
-scratch designer2_processing_test -nocleanup \
$meso1_dicoms designer1_test.mif

# A test with a purposefully wrong input
python designer \
-scratch designer2_processing_test -nocleanup \
BLAH designer1_test.mif

# test different output types.
# nifti output
python designer \
-scratch designer2_processing_test -nocleanup \
$meso1 designer1_test.nii

# zipped output
python designer \
-scratch designer2_processing_test -nocleanup \
$meso1 designer1_test.nii.gz

# test the denoising
# defaults
python designer \
-denoise \
-scratch designer2_processing_test -nocleanup \
$meso1,$meso2 designer1_test.mif




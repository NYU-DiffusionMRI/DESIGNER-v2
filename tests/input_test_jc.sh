
meso=/Users/jc8329/Desktop/test/M0033/M0033_043YF_DIFF_meso.nii
pa=/Users/jc8329/Desktop/test/M0033/M0033_043YF_DIFF_meso_PA.nii
research=/Users/jc8329/Desktop/test/M0033/M0033_043YF_DIFF_meso_research.nii

meso_=/Users/jc8329/Desktop/test/M0033/temp/M0033_043YF_DIFF_meso.nii
pa_=/Users/jc8329/Desktop/test/M0033/temp/M0033_043YF_DIFF_meso_PA.nii
research_=/Users/jc8329/Desktop/test/M0033/temp/M0033_043YF_DIFF_meso_research.nii

out=/Users/jc8329/Desktop/test/M0033
meso_json_wrong=/Users/jc8329/Desktop/test/M0033/M0033_043YF_DIFF_meso_.json
meso_json=/Users/jc8329/Desktop/test/M0033/M0033_043YF_DIFF_meso.json
r_json_wrong=/Users/jc8329/Desktop/test/M0033/M0033_043YF_DIFF_meso_research_.json
r_json=/Users/jc8329/Desktop/test/M0033/M0033_043YF_DIFF_meso_research.json

#test

    designer \
    -eddy -rpe_pair $pa \
    -mask -nocleanup \
    -scratch $out/processing_json \
    $meso,$research $out/dwi_designer_json.nii

    designer \
    -eddy -rpe_pair $pa_ \
    -pf 6/8 -pe_dir j- \
    -mask -nocleanup \
    -scratch $out/temp/processing_nojson \
    $meso_,$research_ $out/dwi_designer_nojson.nii

# ==================================================================================================
#single input file
#no user defined args
# designer $meso $out

#multiple series input
#no user defined args
# designer $meso,$research $out
# designer -bids $meso_json_wrong,$r_json $meso,$research $out

#single input file
#w user defined args
# designer -pf 6/8 -pe_dir j- $meso $out #right
# designer -pf 7/8 -pe_dir 1 $meso $out #wrong

#multiple series input
#w user defined args
# designer -pf 6/8 -pe_dir j- $meso,$research $out #right
# designer -pf 7/8 -pe_dir 1 $meso,$research $out #wrong


# ==================================================================================================
dwi=/Users/jc8329/Desktop/test/dwi/sub-03286_ses-NOT1ACH001_dir-AP_dwi.nii.gz
out=/Users/jc8329/Desktop/test/dwi

#single input file
#no user defined args
# designer -degibbs -pf 1 $dwi $out

#single input file
#w user defined args
# designer -pf 1 -pe_dir j $dwi $out #right
# designer -pf 7/8 -pe_dir j- $dwi $out #wrong

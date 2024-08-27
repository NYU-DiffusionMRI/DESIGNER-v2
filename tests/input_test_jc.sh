
meso=/Users/jc8329/Desktop/test/M0033/M0033_043YF_DIFF_meso.nii
research=/Users/jc8329/Desktop/test/M0033/M0033_043YF_DIFF_meso_research.nii
out=/Users/jc8329/Desktop/test/M0033/dwi.nii
meso_json_wrong=/Users/jc8329/Desktop/test/M0033/M0033_043YF_DIFF_meso_.json
meso_json=/Users/jc8329/Desktop/test/M0033/M0033_043YF_DIFF_meso.json
r_json_wrong=/Users/jc8329/Desktop/test/M0033/M0033_043YF_DIFF_meso_research_.json
r_json=/Users/jc8329/Desktop/test/M0033/M0033_043YF_DIFF_meso_research.json

#test

#nifti
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

dwi=/Users/jc8329/Desktop/test/dwi/sub-03286_ses-NOT1ACH001_dir-AP_dwi.nii.gz
out=/Users/jc8329/Desktop/test/dwi
#nifti
#single input file
#no user defined args
designer -degibbs -pf 1 $dwi $out

#single input file
#w user defined args
# designer -pf 1 -pe_dir j $dwi $out #right
# designer -pf 7/8 -pe_dir j- $dwi $out #wrong

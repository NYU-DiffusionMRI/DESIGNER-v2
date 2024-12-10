
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
# root=/Users/jc8329/Desktop/test/M0033
# meso=/data/M0033_043YF_DIFF_meso.nii
# pa=/data/M0033_043YF_DIFF_meso_PA.nii
# research=/data/M0033_043YF_DIFF_meso_research.nii
#     docker run --rm -it -v ${root}:/data test_petable/designer2:main designer \
#     -eddy -rpe_pair $pa \
#     -mask -nocleanup \
#     -scratch /data/processing_json \
#     $meso,$research /data/dwi_designer_json.nii

# meso=/data/temp/M0033_043YF_DIFF_meso.nii
# pa=/data/temp/M0033_043YF_DIFF_meso_PA.nii
# research=/data/temp/M0033_043YF_DIFF_meso_research.nii
#     docker run --rm -it -v ${root}:/data test_petable/designer2:main designer \
#     -eddy -rpe_pair $pa_ \
#     -pf 6/8 -pe_dir j- \
#     -mask -nocleanup \
#     -scratch /data/temp/processing_nojson \
#     $meso_,$research_ /data/dwi_designer_nojson.nii

folder=M0033
docker run --rm --platform linux/amd64 -it -v /Users/jc8329/Desktop/test/$folder:/data test_fit/designer2:main \
    tmi -DKI -DTI -mask /data/brain_mask.nii \
    /data/dwi_designer.mif /data/params_test

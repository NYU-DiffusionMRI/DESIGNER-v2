
input=/Users/benaron/Documents/datasets/1.5T/NIFTI/DIFF_2p0_12dir_B1000_22.nii
output=/Users/benaron/Documents/datasets/1.5T/output_test

# Test 1
# single input file
# no user defined args

designer $input $output
# we expect: TE 0.094, PE_DIR j-, PF 0.75


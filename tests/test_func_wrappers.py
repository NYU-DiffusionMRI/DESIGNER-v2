import pytest
from designer.func_wrappers import run_mppca

def test_run_mppca():
    # Replace with actual test
    args_extent = "5,5,5"
    args_phase = None
    args_shrinkage = "threshold"
    args_algorithm = "cordero-grande"
    dwi_metadata = {"stride_3dim": "-1,+2,+3,+4"}
    
    # Call the function and check results
    run_mppca(args_extent, args_phase, args_shrinkage, args_algorithm, dwi_metadata)
    # Add assertions to verify the expected behavior
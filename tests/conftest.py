"""
Pytest configuration and fixtures for DESIGNER-v2 integration tests.
"""
import pytest
import tempfile
import shutil
import os
from pathlib import Path
import numpy as np
import nibabel as nib


@pytest.fixture(scope="session")
def test_data_dir():
    """Create a temporary directory with test data for the session."""
    temp_dir = tempfile.mkdtemp(prefix="designer_test_")
    yield Path(temp_dir)
    # Cleanup after tests
    shutil.rmtree(temp_dir, ignore_errors=True)


@pytest.fixture
def sample_dwi_data(test_data_dir):
    """Create sample DWI test data for integration tests."""
    # Create a minimal DWI dataset with realistic dimensions but synthetic data
    dwi_shape = (32, 32, 16, 30)  # Small volume with 30 directions
    dwi_data = np.random.randn(*dwi_shape).astype(np.float32)
    
    # Create realistic b-values and b-vectors
    b_values = np.concatenate([
        np.zeros(3),  # 3 b=0 volumes
        np.full(27, 1000)  # 27 b=1000 volumes
    ])
    
    # Create random but normalized b-vectors
    b_vectors = np.zeros((30, 3))
    b_vectors[3:] = np.random.randn(27, 3)
    # Normalize non-zero b-vectors
    for i in range(3, 30):
        norm = np.linalg.norm(b_vectors[i])
        if norm > 0:
            b_vectors[i] /= norm
    
    # Create NIfTI image
    affine = np.eye(4)
    affine[0, 0] = affine[1, 1] = affine[2, 2] = 2.0  # 2mm isotropic
    dwi_img = nib.Nifti1Image(dwi_data, affine)
    
    # Save files
    dwi_file = test_data_dir / "test_dwi.nii.gz"
    bval_file = test_data_dir / "test_dwi.bval"
    bvec_file = test_data_dir / "test_dwi.bvec"
    
    nib.save(dwi_img, str(dwi_file))
    
    # Save b-values and b-vectors in FSL format
    np.savetxt(bval_file, b_values.reshape(1, -1), fmt='%d')
    np.savetxt(bvec_file, b_vectors.T, fmt='%.6f')
    
    return {
        'dwi': str(dwi_file),
        'bval': str(bval_file),
        'bvec': str(bvec_file),
        'shape': dwi_shape,
        'b_values': b_values,
        'b_vectors': b_vectors
    }


@pytest.fixture
def sample_phase_encoding_data(test_data_dir):
    """Create sample phase encoding data for eddy correction tests."""
    # Create a small b=0 image with reverse phase encoding
    pa_shape = (32, 32, 16)
    pa_data = np.random.randn(*pa_shape).astype(np.float32)
    
    affine = np.eye(4)
    affine[0, 0] = affine[1, 1] = affine[2, 2] = 2.0
    pa_img = nib.Nifti1Image(pa_data, affine)
    
    pa_file = test_data_dir / "test_pa.nii.gz"
    nib.save(pa_img, str(pa_file))
    
    return str(pa_file)


@pytest.fixture
def output_dir(test_data_dir):
    """Create a temporary output directory for each test."""
    output_path = test_data_dir / "output"
    output_path.mkdir(exist_ok=True)
    return str(output_path)


@pytest.fixture
def mock_environment(monkeypatch):
    """Mock environment variables and paths for testing."""
    # Mock FSL environment
    monkeypatch.setenv("FSLDIR", "/usr/local/fsl")
    monkeypatch.setenv("FSLOUTPUTTYPE", "NIFTI_GZ")
    
    # Mock MRtrix3 environment  
    monkeypatch.setenv("PATH", "/usr/local/mrtrix3/bin:" + os.environ.get("PATH", ""))
    
    return True


@pytest.fixture
def designer_args_basic():
    """Basic arguments for designer without preprocessing steps."""
    return {
        'denoise': False,
        'degibbs': False,
        'rician': False,
        'eddy': False,
        'normalize': False,
        'mask': False,
        'b1correct': False,
        'n_cores': 1,
        'datatype': None
    }


@pytest.fixture
def designer_args_full():
    """Full preprocessing arguments for designer."""
    return {
        'denoise': True,
        'degibbs': True,
        'rician': True,
        'eddy': False,  # Skip eddy for basic tests (requires FSL)
        'normalize': True,
        'mask': True,
        'b1correct': False,
        'n_cores': 1,
        'extent': '5,5,5',
        'shrinkage': 'threshold',
        'algorithm': 'cordero-grande',
        'pf': '6/8',
        'pe_dir': '1',
        'datatype': 'float32'
    }

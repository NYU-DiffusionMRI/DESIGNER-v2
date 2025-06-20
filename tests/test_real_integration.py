"""
Real integration tests for DESIGNER-v2 that run with actual dependencies.

These tests translate the integration.sh script into Python tests that run the actual
DESIGNER pipeline with real MRtrix3, FSL, and other dependencies. They use synthetic
test data but exercise the complete processing pipeline.

Requirements:
- MRtrix3 must be installed and available in PATH
- FSL may be required for some tests (eddy correction)
- ANTs may be required for some tests (motion correction)

To run these tests:
    pytest tests/test_real_integration.py -v

To skip if dependencies are missing:
    pytest tests/test_real_integration.py -v -k "not real_deps"
"""
import pytest
import os
import sys
import tempfile
import subprocess
import shutil
from pathlib import Path
import numpy as np
import nibabel as nib
import time
import json
import logging
from typing import Dict, List, Optional, Tuple

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Add the project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent))


def check_dependency(cmd: str) -> bool:
    """Check if a command-line dependency is available."""
    try:
        result = subprocess.run([cmd, '--help'], 
                              capture_output=True, 
                              timeout=10)
        return result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return False


def check_mrtrix3() -> bool:
    """Check if MRtrix3 is available."""
    return check_dependency('mrinfo')


def check_fsl() -> bool:
    """Check if FSL is available."""
    return check_dependency('fslinfo')


def check_ants() -> bool:
    """Check if ANTs is available."""
    return check_dependency('antsRegistration')


class RealTestDataGenerator:
    """Generate realistic test data for integration testing."""
    
    @staticmethod
    def create_dwi_dataset(output_dir: Path, 
                          name: str,
                          shape: Tuple[int, int, int, int] = (64, 64, 32, 30),
                          bvals: Optional[List[int]] = None,
                          noise_level: float = 0.1) -> Dict[str, str]:
        """
        Create a realistic DWI dataset with proper NIfTI structure.
        
        Args:
            output_dir: Directory to save files
            name: Base name for files
            shape: Image dimensions (x, y, z, volumes)
            bvals: B-values (default: [0]*3 + [1000]*27)
            noise_level: Gaussian noise level
            
        Returns:
            Dictionary with file paths
        """
        if bvals is None:
            bvals = [0] * 3 + [1000] * (shape[3] - 3)
        
        # Create realistic DWI signal
        # Start with a simple tensor model
        S0 = 1000  # Reference signal
        D = 1e-3   # Diffusivity (typical brain tissue)
        
        # Generate DWI data with exponential decay
        dwi_data = np.zeros(shape, dtype=np.float32)
        
        for vol in range(shape[3]):
            if bvals[vol] == 0:
                # b=0 image
                signal = S0
            else:
                # Exponential decay with some spatial variation
                x, y, z = np.meshgrid(range(shape[0]), range(shape[1]), range(shape[2]), indexing='ij')
                # Add some spatial structure
                tissue_factor = 0.8 + 0.4 * np.exp(-((x-shape[0]//2)**2 + (y-shape[1]//2)**2) / (shape[0]*shape[1]/8))
                signal = S0 * tissue_factor * np.exp(-bvals[vol] * D)
            
            # Add Rician noise
            noise_real = np.random.normal(0, noise_level * signal, shape[:3])
            noise_imag = np.random.normal(0, noise_level * signal, shape[:3])
            signal_noisy = np.sqrt((signal + noise_real)**2 + noise_imag**2)
            
            dwi_data[:, :, :, vol] = signal_noisy
        
        # Create NIfTI image with proper header
        affine = np.diag([2.0, 2.0, 2.0, 1.0])  # 2mm isotropic
        dwi_img = nib.Nifti1Image(dwi_data, affine)
        
        # Set realistic header information
        dwi_img.header.set_xyzt_units('mm', 'sec')
        
        # Save files
        dwi_file = output_dir / f"{name}.nii.gz"
        bval_file = output_dir / f"{name}.bval"
        bvec_file = output_dir / f"{name}.bvec"
        
        nib.save(dwi_img, str(dwi_file))
        
        # Create b-vectors (random unit vectors for non-zero b-values)
        bvecs = np.zeros((len(bvals), 3))
        for i, bval in enumerate(bvals):
            if bval > 0:
                # Random direction
                vec = np.random.randn(3)
                bvecs[i] = vec / np.linalg.norm(vec)
        
        # Save in FSL format
        np.savetxt(bval_file, bvals, fmt='%d')
        np.savetxt(bvec_file, bvecs.T, fmt='%.6f')
        
        return {
            'dwi': str(dwi_file),
            'bval': str(bval_file),
            'bvec': str(bvec_file),
            'shape': shape,
            'bvals': bvals,
            'bvecs': bvecs
        }
    
    @staticmethod
    def create_phase_encoding_pair(output_dir: Path, 
                                  dwi_shape: Tuple[int, int, int]) -> str:
        """Create a reverse phase encoding b=0 image."""
        # Create b=0 volume with different distortion pattern
        pa_data = np.random.randn(*dwi_shape).astype(np.float32) * 100 + 1000
        
        affine = np.diag([2.0, 2.0, 2.0, 1.0])
        pa_img = nib.Nifti1Image(pa_data, affine)
        
        pa_file = output_dir / "pa_b0.nii.gz"
        nib.save(pa_img, str(pa_file))
        
        return str(pa_file)


class TestDesignerRealIntegration:
    """Real integration tests based on integration.sh script."""
    
    @pytest.fixture(scope="class")
    def test_workspace(self):
        """Create a temporary workspace for all tests."""
        workspace = tempfile.mkdtemp(prefix="designer_real_test_")
        yield Path(workspace)
        # Cleanup
        shutil.rmtree(workspace, ignore_errors=True)
    
    @pytest.fixture
    def test_data(self, test_workspace):
        """Generate test data for each test."""
        data_dir = test_workspace / "data"
        data_dir.mkdir(exist_ok=True)
        
        generator = RealTestDataGenerator()
        
        # Create primary DWI dataset (equivalent to M9734_074YM_DIFF_meso.nii)
        meso1 = generator.create_dwi_dataset(
            data_dir, 
            "meso1",
            shape=(32, 32, 16, 30),  # Smaller for faster testing
            noise_level=0.05
        )
        
        # Create secondary DWI dataset (equivalent to M9734_074YM__DIFF_meso_research.nii)
        meso2 = generator.create_dwi_dataset(
            data_dir,
            "meso2", 
            shape=(32, 32, 16, 30),
            bvals=[0] * 3 + [2000] * 27,  # Different b-value
            noise_level=0.05
        )
        
        # Create reverse phase encoding image
        pa_file = generator.create_phase_encoding_pair(
            data_dir,
            (32, 32, 16)
        )
        
        return {
            'meso1': meso1,
            'meso2': meso2,
            'pa_file': pa_file,
            'data_dir': data_dir
        }
    
    def run_designer_command(self, args: List[str], timeout: int = 300) -> subprocess.CompletedProcess:
        """
        Run designer command and return result.
        
        Args:
            args: Command line arguments for designer
            timeout: Timeout in seconds
            
        Returns:
            CompletedProcess result
        """
        # Import here to avoid import errors if designer2 isn't available
        designer_script = Path(__file__).parent.parent / "designer2" / "designer.py"
        cmd = [sys.executable, str(designer_script)] + args
        logger.info(f"Running: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
                cwd=Path(__file__).parent.parent
            )
            return result
        except subprocess.TimeoutExpired:
            pytest.fail(f"Command timed out after {timeout} seconds")
    
    def assert_output_files_exist(self, output_path: str, expected_extensions: List[str]):
        """Assert that expected output files exist."""
        base_path = Path(output_path)
        
        for ext in expected_extensions:
            if ext.startswith('.'):
                file_path = base_path.with_suffix(ext)
            else:
                file_path = Path(f"{output_path}.{ext}")
            
            assert file_path.exists(), f"Expected output file not found: {file_path}"
            assert file_path.stat().st_size > 0, f"Output file is empty: {file_path}"
    
    def assert_nifti_properties(self, nifti_path: str, expected_shape: Optional[Tuple] = None):
        """Assert properties of output NIfTI file."""
        assert os.path.exists(nifti_path), f"NIfTI file not found: {nifti_path}"
        
        try:
            img = nib.load(nifti_path)
            
            # Check that image loads successfully
            assert img is not None, "Failed to load NIfTI image"
            
            # Check shape if provided
            if expected_shape:
                assert img.shape == expected_shape, f"Shape mismatch: got {img.shape}, expected {expected_shape}"
            
            # Check data is not all zeros or NaN
            data = img.get_fdata()
            assert not np.all(data == 0), "Output image is all zeros"
            assert not np.any(np.isnan(data)), "Output image contains NaN values"
            
            logger.info(f"NIfTI validation passed for {nifti_path}: shape={img.shape}, dtype={data.dtype}")
            
        except Exception as e:
            pytest.fail(f"NIfTI validation failed for {nifti_path}: {e}")
    
    def assert_bval_bvec_files(self, base_path: str, expected_volumes: int):
        """Assert that bval and bvec files are correctly formatted."""
        bval_path = f"{base_path}.bval"
        bvec_path = f"{base_path}.bvec"
        
        assert os.path.exists(bval_path), f"bval file not found: {bval_path}"
        assert os.path.exists(bvec_path), f"bvec file not found: {bvec_path}"
        
        # Check bval file
        bvals = np.loadtxt(bval_path)
        if bvals.ndim == 0:
            bvals = np.array([bvals])
        assert len(bvals) == expected_volumes, f"bval count mismatch: got {len(bvals)}, expected {expected_volumes}"
        
        # Check bvec file
        bvecs = np.loadtxt(bvec_path)
        if bvecs.ndim == 1:
            bvecs = bvecs.reshape(3, -1)
        assert bvecs.shape[1] == expected_volumes, f"bvec count mismatch: got {bvecs.shape[1]}, expected {expected_volumes}"
        
        logger.info(f"bval/bvec validation passed: {expected_volumes} volumes")

    # Test 1: Single input series (basic processing)
    @pytest.mark.real_deps
    @pytest.mark.skipif(not check_mrtrix3(), reason="MRtrix3 not available")
    def test_single_input_basic_processing(self, test_data, test_workspace):
        """Test 1: Single input series - basic processing without preprocessing steps."""
        output_path = test_workspace / "single_input_test"
        
        args = [
            test_data['meso1']['dwi'],
            str(output_path),
            '-scratch', str(test_workspace / 'scratch1'),
            '-nocleanup'
        ]
        
        result = self.run_designer_command(args)
        
        # Log output for debugging
        logger.info(f"STDOUT: {result.stdout}")
        if result.stderr:
            logger.warning(f"STDERR: {result.stderr}")
        
        # Assert command succeeded
        assert result.returncode == 0, f"Command failed with return code {result.returncode}: {result.stderr}"
        
        # Assert output files exist
        self.assert_output_files_exist(str(output_path), ['.nii.gz', '.bval', '.bvec'])
        
        # Assert NIfTI properties
        self.assert_nifti_properties(f"{output_path}.nii.gz", expected_shape=(32, 32, 16, 30))
        
        # Assert bval/bvec files
        self.assert_bval_bvec_files(str(output_path), 30)
        
        logger.info("✓ Single input basic processing test passed")

    # Test 2: Multiple input series
    @pytest.mark.real_deps
    @pytest.mark.skipif(not check_mrtrix3(), reason="MRtrix3 not available")
    def test_multiple_input_series(self, test_data, test_workspace):
        """Test 2: Multiple input series concatenation."""
        output_path = test_workspace / "multi_input_test"
        
        # Combine input files with comma separation
        input_files = f"{test_data['meso1']['dwi']},{test_data['meso2']['dwi']}"
        
        args = [
            input_files,
            str(output_path),
            '-scratch', str(test_workspace / 'scratch2'),
            '-nocleanup'
        ]
        
        result = self.run_designer_command(args)
        
        # Log output for debugging
        logger.info(f"STDOUT: {result.stdout}")
        if result.stderr:
            logger.warning(f"STDERR: {result.stderr}")
        
        # Assert command succeeded
        assert result.returncode == 0, f"Command failed with return code {result.returncode}: {result.stderr}"
        
        # Assert output files exist
        self.assert_output_files_exist(str(output_path), ['.nii.gz', '.bval', '.bvec'])
        
        # Should have concatenated volumes (30 + 30 = 60)
        self.assert_nifti_properties(f"{output_path}.nii.gz", expected_shape=(32, 32, 16, 60))
        
        # Assert bval/bvec files for combined data
        self.assert_bval_bvec_files(str(output_path), 60)
        
        logger.info("✓ Multiple input series test passed")

    # Test 3: Processing with denoising
    @pytest.mark.real_deps
    @pytest.mark.skipif(not check_mrtrix3(), reason="MRtrix3 not available")
    def test_denoising_processing(self, test_data, test_workspace):
        """Test 3: Single input with denoising enabled."""
        output_path = test_workspace / "denoise_test"
        
        args = [
            test_data['meso1']['dwi'],
            str(output_path),
            '-denoise',
            '-extent', '3,3,3',  # Smaller extent for faster processing
            '-scratch', str(test_workspace / 'scratch3'),
            '-nocleanup'
        ]
        
        result = self.run_designer_command(args, timeout=600)  # Longer timeout for denoising
        
        # Log output for debugging
        logger.info(f"STDOUT: {result.stdout}")
        if result.stderr:
            logger.warning(f"STDERR: {result.stderr}")
        
        # Assert command succeeded
        assert result.returncode == 0, f"Command failed with return code {result.returncode}: {result.stderr}"
        
        # Assert output files exist
        self.assert_output_files_exist(str(output_path), ['.nii.gz', '.bval', '.bvec'])
        
        # Assert NIfTI properties
        self.assert_nifti_properties(f"{output_path}.nii.gz", expected_shape=(32, 32, 16, 30))
        
        # Assert bval/bvec files
        self.assert_bval_bvec_files(str(output_path), 30)
        
        logger.info("✓ Denoising processing test passed")

    # Test 4: Error handling - invalid input
    def test_invalid_input_handling(self, test_workspace):
        """Test 4: Error handling with invalid input file."""
        output_path = test_workspace / "error_test"
        
        args = [
            "/nonexistent/file.nii",  # Invalid file
            str(output_path),
            '-scratch', str(test_workspace / 'scratch4'),
            '-nocleanup'
        ]
        
        result = self.run_designer_command(args)
        
        # Log output for debugging
        logger.info(f"STDOUT: {result.stdout}")
        logger.info(f"STDERR: {result.stderr}")
        
        # Assert command failed as expected
        assert result.returncode != 0, "Command should have failed with invalid input"
        
        logger.info("✓ Invalid input handling test passed")

    # Test 5: Custom bval/bvec files
    @pytest.mark.real_deps
    @pytest.mark.skipif(not check_mrtrix3(), reason="MRtrix3 not available")
    def test_custom_bval_bvec_files(self, test_data, test_workspace):
        """Test 5: Using custom bval/bvec files."""
        output_path = test_workspace / "custom_bval_test"
        
        args = [
            test_data['meso1']['dwi'],
            str(output_path),
            '-fslbval', test_data['meso1']['bval'],
            '-fslbvec', test_data['meso1']['bvec'],
            '-scratch', str(test_workspace / 'scratch5'),
            '-nocleanup'
        ]
        
        result = self.run_designer_command(args)
        
        # Log output for debugging
        logger.info(f"STDOUT: {result.stdout}")
        if result.stderr:
            logger.warning(f"STDERR: {result.stderr}")
        
        # Assert command succeeded
        assert result.returncode == 0, f"Command failed with return code {result.returncode}: {result.stderr}"
        
        # Assert output files exist
        self.assert_output_files_exist(str(output_path), ['.nii.gz', '.bval', '.bvec'])
        
        logger.info("✓ Custom bval/bvec files test passed")

    # Test 6: Performance benchmark
    @pytest.mark.performance
    @pytest.mark.real_deps
    @pytest.mark.skipif(not check_mrtrix3(), reason="MRtrix3 not available")
    def test_processing_performance(self, test_data, test_workspace):
        """Test 6: Basic performance benchmark."""
        output_path = test_workspace / "performance_test"
        
        start_time = time.time()
        
        args = [
            test_data['meso1']['dwi'],
            str(output_path),
            '-scratch', str(test_workspace / 'scratch6'),
            '-nocleanup'
        ]
        
        result = self.run_designer_command(args)
        
        elapsed_time = time.time() - start_time
        
        # Log output for debugging
        logger.info(f"STDOUT: {result.stdout}")
        if result.stderr:
            logger.warning(f"STDERR: {result.stderr}")
        
        # Assert command succeeded
        assert result.returncode == 0, f"Command failed with return code {result.returncode}: {result.stderr}"
        
        # Performance assertion - basic processing should complete reasonably quickly
        # Allow more time for real processing
        assert elapsed_time < 300, f"Basic processing took too long: {elapsed_time:.2f}s"
        
        logger.info(f"✓ Performance test passed: {elapsed_time:.2f}s")

    # Test 7: Output validation
    @pytest.mark.real_deps
    @pytest.mark.skipif(not check_mrtrix3(), reason="MRtrix3 not available")
    def test_output_data_validation(self, test_data, test_workspace):
        """Test 7: Detailed validation of output data quality."""
        output_path = test_workspace / "validation_test"
        
        args = [
            test_data['meso1']['dwi'],
            str(output_path),
            '-scratch', str(test_workspace / 'scratch7'),
            '-nocleanup'
        ]
        
        result = self.run_designer_command(args)
        
        # Log output for debugging
        logger.info(f"STDOUT: {result.stdout}")
        if result.stderr:
            logger.warning(f"STDERR: {result.stderr}")
        
        # Assert command succeeded
        assert result.returncode == 0, f"Command failed with return code {result.returncode}: {result.stderr}"
        
        # Load and validate output data
        output_img = nib.load(f"{output_path}.nii.gz")
        output_data = output_img.get_fdata()
        
        # Load input data for comparison
        input_img = nib.load(test_data['meso1']['dwi'])
        input_data = input_img.get_fdata()
        
        # Validate output properties
        assert output_data.shape == input_data.shape, "Shape should be preserved"
        assert output_data.dtype in [np.float32, np.float64], "Output should be floating point"
        
        # Validate signal range is reasonable
        assert np.all(output_data >= 0), "Output should be non-negative"
        assert np.max(output_data) > 0, "Output should not be all zeros"
        
        # Validate b=0 volumes have higher signal than DWI volumes
        bvals = np.loadtxt(f"{output_path}.bval")
        b0_indices = np.where(bvals == 0)[0]
        dwi_indices = np.where(bvals > 0)[0]
        
        if len(b0_indices) > 0 and len(dwi_indices) > 0:
            b0_mean = np.mean(output_data[:, :, :, b0_indices])
            dwi_mean = np.mean(output_data[:, :, :, dwi_indices])
            assert b0_mean > dwi_mean, "b=0 signal should be higher than DWI signal"
        
        logger.info("✓ Output data validation test passed")


# Utility function to run all real integration tests
def run_real_integration_tests():
    """Run all real integration tests."""
    pytest.main([
        __file__,
        "-v",
        "--tb=short",
        "-m", "real_deps",
        "--cov=designer2",
        "--cov=lib",
        "--cov-report=term-missing",
        "--cov-report=html:htmlcov_real",
        "--cov-report=xml:coverage_real.xml"
    ])


if __name__ == "__main__":
    run_real_integration_tests()

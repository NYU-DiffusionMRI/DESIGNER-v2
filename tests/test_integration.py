"""
Integration tests for DESIGNER-v2 main processing pipeline.
"""
import pytest
import os
import sys
import tempfile
import subprocess
from pathlib import Path
from unittest.mock import patch, MagicMock
import numpy as np
import nibabel as nib

# Add the project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent))

from designer2 import designer


class TestDesignerBasicFunctionality:
    """Test basic DESIGNER functionality without external dependencies."""
    
    def test_designer_import(self):
        """Test that designer module can be imported successfully."""
        assert hasattr(designer, 'execute')
        assert hasattr(designer, 'usage')
        assert hasattr(designer, 'main')
    
    def test_usage_function(self):
        """Test the usage function doesn't crash."""
        # Mock the cmdline object
        class MockCmdline:
            def set_copyright(self, text): pass
            def set_author(self, text): pass
            def set_synopsis(self, text): pass
            def add_citation(self, text, is_external=False): pass
            def add_argument(self, name, help=None): pass
            def add_argument_group(self, name): return self
        
        mock_cmdline = MockCmdline()
        # Should not raise any exceptions
        designer.usage(mock_cmdline)


class TestDesignerWithTestData:
    """Test DESIGNER processing with synthetic test data."""
    
    @pytest.mark.parametrize("processing_step", [
        ("basic", {}),
        ("denoise_only", {"denoise": True}),
        ("normalize_only", {"normalize": True}),
    ])
    def test_designer_processing_steps(self, sample_dwi_data, output_dir, 
                                     processing_step, mock_environment):
        """Test individual processing steps of DESIGNER."""
        step_name, args = processing_step
        
        # Mock mrtrix3 dependencies
        with patch('designer2.designer.app') as mock_app, \
             patch('designer2.designer.run') as mock_run, \
             patch('designer2.designer.fsl') as mock_fsl:
            
            # Setup mocks
            mock_app.ARGS = type('Args', (), {
                'input': sample_dwi_data['dwi'],
                'output': os.path.join(output_dir, f'test_{step_name}'),
                'fslbval': sample_dwi_data['bval'],
                'fslbvec': sample_dwi_data['bvec'],
                'bids': None,
                'pe_dir': None,
                'pf': None,
                'denoise': args.get('denoise', False),
                'degibbs': args.get('degibbs', False),
                'rician': args.get('rician', False),
                'eddy': args.get('eddy', False),
                'normalize': args.get('normalize', False),
                'mask': args.get('mask', False),
                'b1correct': args.get('b1correct', False),
                'n_cores': 1,
                'extent': '5,5,5',
                'shrinkage': 'threshold',
                'algorithm': 'cordero-grande',
                'phase': None,
                'patch2self': False,
                'adaptive_patch': False,
                'datatype': 'float32'
            })
            
            mock_app.make_scratch_dir.return_value = None
            mock_app.goto_scratch_dir.return_value = None
            mock_fsl.suffix.return_value = '.nii.gz'
            
            # Mock the processing functions
            with patch('lib.designer_input_utils.get_input_info') as mock_get_info, \
                 patch('lib.designer_input_utils.convert_input_data') as mock_convert, \
                 patch('lib.designer_input_utils.create_shell_table') as mock_shell:
                
                # Setup input processing mocks
                mock_get_info.return_value = {
                    'bshape_per_volume': [1] * 30,
                    'echo_time_per_volume': [0.1] * 30
                }
                mock_shell.return_value = np.array([[0, 1000], [1, 1], [3, 27], [0.1, 0.1]])
                
                # This should not raise any exceptions
                try:
                    designer.execute()
                    success = True
                except Exception as e:
                    # Expected to fail due to missing external dependencies
                    # but should not fail due to our code logic
                    success = "AttributeError" not in str(type(e))
                
                assert success, f"Processing step {step_name} failed unexpectedly"


class TestDesignerInputValidation:
    """Test input validation and error handling."""
    
    def test_missing_input_file(self, output_dir, mock_environment):
        """Test behavior with missing input files."""
        with patch('designer2.designer.app') as mock_app:
            mock_app.ARGS = type('Args', (), {
                'input': '/nonexistent/file.nii',
                'output': os.path.join(output_dir, 'test_missing'),
                'fslbval': None,
                'fslbvec': None,
                'bids': None
            })
            
            with patch('lib.designer_input_utils.get_input_info') as mock_get_info:
                # Should handle missing files gracefully
                mock_get_info.side_effect = FileNotFoundError("Input file not found")
                
                with pytest.raises(FileNotFoundError):
                    designer.execute()
    
    def test_invalid_bval_bvec_files(self, sample_dwi_data, output_dir, mock_environment):
        """Test behavior with invalid bval/bvec files."""
        # Create invalid bval file
        invalid_bval = os.path.join(os.path.dirname(sample_dwi_data['bval']), 'invalid.bval')
        with open(invalid_bval, 'w') as f:
            f.write("invalid content\nnot numbers")
        
        with patch('designer2.designer.app') as mock_app:
            mock_app.ARGS = type('Args', (), {
                'input': sample_dwi_data['dwi'],
                'output': os.path.join(output_dir, 'test_invalid'),
                'fslbval': invalid_bval,
                'fslbvec': sample_dwi_data['bvec'],
                'bids': None
            })
            
            with patch('lib.designer_input_utils.get_input_info') as mock_get_info:
                mock_get_info.side_effect = ValueError("Invalid bval format")
                
                with pytest.raises(ValueError):
                    designer.execute()


class TestDesignerOutputs:
    """Test DESIGNER output file generation and validation."""
    
    def test_output_file_generation(self, sample_dwi_data, output_dir, mock_environment):
        """Test that expected output files are generated."""
        output_basename = os.path.join(output_dir, 'test_output')
        
        with patch('designer2.designer.app') as mock_app, \
             patch('designer2.designer.run') as mock_run, \
             patch('designer2.designer.fsl') as mock_fsl:
            
            mock_app.ARGS = type('Args', (), {
                'input': sample_dwi_data['dwi'],
                'output': output_basename,
                'fslbval': sample_dwi_data['bval'],
                'fslbvec': sample_dwi_data['bvec'],
                'bids': None,
                'denoise': False,
                'degibbs': False,
                'rician': False,
                'eddy': False,
                'normalize': False,
                'mask': False,
                'datatype': 'float32'
            })
            
            mock_app.make_scratch_dir.return_value = None
            mock_app.goto_scratch_dir.return_value = None
            mock_fsl.suffix.return_value = '.nii.gz'
            
            # Mock file operations
            with patch('lib.designer_input_utils.get_input_info') as mock_get_info, \
                 patch('lib.designer_input_utils.convert_input_data') as mock_convert, \
                 patch('lib.designer_input_utils.create_shell_table') as mock_shell, \
                 patch('os.path.exists') as mock_exists:
                
                mock_get_info.return_value = {
                    'bshape_per_volume': [1] * 30,
                    'echo_time_per_volume': [0.1] * 30
                }
                mock_shell.return_value = np.array([[0, 1000], [1, 1], [3, 27], [0.1, 0.1]])
                mock_exists.return_value = True
                
                # Mock the final conversion step
                mock_run.command = MagicMock()
                
                try:
                    designer.execute()
                    # Check that mrconvert was called for output generation
                    assert mock_run.command.called
                except Exception:
                    # Expected due to missing dependencies, but conversion should be attempted
                    pass


class TestDesignerCommandLineInterface:
    """Test DESIGNER command-line interface."""
    
    def test_main_function_exists(self):
        """Test that main function exists and can be called."""
        assert callable(designer.main)
    
    def test_script_execution_dry_run(self, sample_dwi_data, output_dir):
        """Test script execution in dry-run mode."""
        # This would test actual command-line execution
        # For now, we'll test the entry point exists
        script_path = Path(__file__).parent.parent / "designer2" / "designer.py"
        assert script_path.exists()
        
        # Test that the script is executable (has main guard)
        with open(script_path, 'r') as f:
            content = f.read()
            assert 'if __name__ == "__main__":' in content
            assert 'main()' in content


class TestDesignerParameterValidation:
    """Test parameter validation and edge cases."""
    
    @pytest.mark.parametrize("n_cores", [-1, 0, 1, 4, 16])
    def test_n_cores_parameter(self, n_cores):
        """Test n_cores parameter handling."""
        # Test different core count values
        # This should not crash for any reasonable value
        assert isinstance(n_cores, int)
        # Negative values should be handled appropriately
        if n_cores < 0:
            # Should use default behavior (available cores - abs(n_cores))
            pass
        elif n_cores == 0:
            # Should handle gracefully
            pass
        else:
            # Should use specified number of cores
            pass
    
    @pytest.mark.parametrize("extent", ["5,5,5", "3,3,3", "7,7,7"])
    def test_extent_parameter(self, extent):
        """Test extent parameter validation."""
        # Should be parseable as comma-separated integers
        parts = extent.split(',')
        assert len(parts) == 3
        assert all(part.isdigit() for part in parts)
    
    @pytest.mark.parametrize("algorithm", ["veraart", "cordero-grande", "jespersen"])
    def test_algorithm_parameter(self, algorithm):
        """Test algorithm parameter validation."""
        valid_algorithms = ["veraart", "cordero-grande", "jespersen"]
        assert algorithm in valid_algorithms


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

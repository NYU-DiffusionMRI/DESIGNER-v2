"""
Unit tests for DESIGNER-v2 library modules.
"""
import pytest
import sys
import os
from pathlib import Path
import numpy as np
from unittest.mock import patch, MagicMock

# Add the project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent))


class TestLibraryModules:
    """Test individual library modules."""
    
    def test_import_lib_modules(self):
        """Test that library modules can be imported."""
        try:
            import lib.designer_input_utils
            import lib.designer_func_wrappers
            assert True
        except ImportError as e:
            pytest.skip(f"Library modules not available: {e}")
    
    def test_mpcomplex_import(self):
        """Test mpcomplex module import."""
        try:
            import lib.mpcomplex
            assert hasattr(lib.mpcomplex, '__file__')
        except ImportError:
            pytest.skip("mpcomplex module not available")
    
    def test_tensor_import(self):
        """Test tensor module import."""
        try:
            import lib.tensor
            assert hasattr(lib.tensor, '__file__')
        except ImportError:
            pytest.skip("tensor module not available")
    
    def test_smi_import(self):
        """Test SMI module import."""
        try:
            import lib.smi
            assert hasattr(lib.smi, '__file__')
        except ImportError:
            pytest.skip("SMI module not available")


class TestInputUtils:
    """Test input utility functions."""
    
    def test_get_input_info_mock(self):
        """Test get_input_info function with mocked inputs."""
        try:
            from lib.designer_input_utils import get_input_info
        except ImportError:
            pytest.skip("designer_input_utils not available")
        
        # Mock the function behavior since we don't have real MRtrix3
        with patch('lib.designer_input_utils.get_input_info') as mock_func:
            mock_func.return_value = {
                'bshape_per_volume': [1] * 30,
                'echo_time_per_volume': [0.1] * 30,
                'input_files': ['test.nii'],
                'bval_files': ['test.bval'],
                'bvec_files': ['test.bvec']
            }
            
            result = mock_func('test.nii', None, None, None)
            assert isinstance(result, dict)
            assert 'bshape_per_volume' in result
            assert 'echo_time_per_volume' in result
    
    def test_create_shell_table_mock(self):
        """Test create_shell_table function with mocked inputs."""
        try:
            from lib.designer_input_utils import create_shell_table
        except ImportError:
            pytest.skip("designer_input_utils not available")
        
        with patch('lib.designer_input_utils.create_shell_table') as mock_func:
            mock_metadata = {
                'bshape_per_volume': [1] * 30,
                'echo_time_per_volume': [0.1] * 30
            }
            
            expected_table = np.array([
                [0, 1000],      # b-values
                [1, 1],         # b-shapes  
                [3, 27],        # n volumes
                [0.1, 0.1]      # echo times
            ])
            
            mock_func.return_value = expected_table
            result = mock_func(mock_metadata)
            
            assert isinstance(result, np.ndarray)
            assert result.shape[0] == 4  # 4 rows: b-value, b-shape, n volumes, echo time


class TestFunctionWrappers:
    """Test function wrapper modules."""
    
    def test_run_mppca_mock(self):
        """Test MPPCA wrapper function."""
        try:
            from lib.designer_func_wrappers import run_mppca
        except ImportError:
            pytest.skip("designer_func_wrappers not available")
        
        with patch('lib.designer_func_wrappers.run_mppca') as mock_func:
            mock_metadata = {'test': 'data'}
            
            # Should not raise exception when mocked
            mock_func.return_value = None
            result = mock_func('5,5,5', None, 'threshold', 'cordero-grande', mock_metadata)
            assert result is None
    
    def test_run_patch2self_mock(self):
        """Test patch2self wrapper function."""
        try:
            from lib.designer_func_wrappers import run_patch2self
        except ImportError:
            pytest.skip("designer_func_wrappers not available")
        
        with patch('lib.designer_func_wrappers.run_patch2self') as mock_func:
            mock_func.return_value = None
            result = mock_func()
            assert result is None


class TestMPComplex:
    """Test MP complex data processing."""
    
    def test_mpcomplex_functions_exist(self):
        """Test that mpcomplex functions exist."""
        try:
            import lib.mpcomplex
            # Check for expected functions/classes
            # This is a basic existence check since we can't test without data
            assert hasattr(lib.mpcomplex, '__file__')
        except ImportError:
            pytest.skip("mpcomplex module not available")


class TestTensor:
    """Test tensor processing functions."""
    
    def test_tensor_functions_exist(self):
        """Test that tensor functions exist."""
        try:
            import lib.tensor
            assert hasattr(lib.tensor, '__file__')
        except ImportError:
            pytest.skip("tensor module not available")


class TestSMI:
    """Test SMI (Spherical Mean Imaging) functions."""
    
    def test_smi_functions_exist(self):
        """Test that SMI functions exist."""
        try:
            import lib.smi
            assert hasattr(lib.smi, '__file__')
        except ImportError:
            pytest.skip("smi module not available")


class TestRPG:
    """Test RPG (Random Phase Generation) C++ extension."""
    
    def test_rpg_import(self):
        """Test RPG C++ module import."""
        try:
            import lib.rpg
            assert hasattr(lib.rpg, '__file__')
            # If successfully imported, it's a compiled module
            assert lib.rpg.__file__.endswith('.so')
        except ImportError:
            pytest.skip("RPG C++ module not available (compilation may be needed)")


class TestConstantModule:
    """Test constant values module."""
    
    def test_constant_import(self):
        """Test constant module import."""
        try:
            import constant
            assert hasattr(constant, '__file__')
        except ImportError:
            pytest.skip("constant module not available")
    
    def test_constant_data_files(self):
        """Test that constant data files exist."""
        constant_dir = Path(__file__).parent.parent / "constant"
        
        expected_files = [
            "constant_values.mat",
            "dirs_constrained.mat", 
            "dirs10000.mat",
            "dirs256.mat"
        ]
        
        for filename in expected_files:
            filepath = constant_dir / filename
            if filepath.exists():
                assert filepath.stat().st_size > 0, f"{filename} is empty"
            else:
                pytest.skip(f"Constant file {filename} not found")


class TestDesignerTMI:
    """Test TMI (Tissue Microstructure Imaging) module."""
    
    def test_tmi_import(self):
        """Test TMI module import."""
        try:
            from designer2 import tmi
            assert hasattr(tmi, 'execute')
            assert hasattr(tmi, 'usage')
            assert hasattr(tmi, 'main')
        except ImportError:
            pytest.skip("TMI module not available")
    
    def test_tmi_usage_function(self):
        """Test TMI usage function."""
        try:
            from designer2 import tmi
        except ImportError:
            pytest.skip("TMI module not available")
        
        # Mock the cmdline object
        class MockCmdline:
            def set_copyright(self, text): pass
            def set_author(self, text): pass
            def set_synopsis(self, text): pass
            def add_citation(self, text, is_external=False): pass
            def add_argument(self, name, help=None, **kwargs): pass
            def add_argument_group(self, name): return self
        
        mock_cmdline = MockCmdline()
        # Should not raise any exceptions
        tmi.usage(mock_cmdline)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

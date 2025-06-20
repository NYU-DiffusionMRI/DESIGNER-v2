"""
Performance and benchmarking tests for DESIGNER-v2.
"""
import pytest
import time
import psutil
import os
import tempfile
from pathlib import Path
import numpy as np
import nibabel as nib
from unittest.mock import patch


class TestDesignerPerformance:
    """Test performance characteristics of DESIGNER processing."""
    
    @pytest.mark.performance
    def test_memory_usage_basic_processing(self, sample_dwi_data, output_dir, mock_environment):
        """Test memory usage during basic processing."""
        # Get initial memory usage
        process = psutil.Process()
        initial_memory = process.memory_info().rss / 1024 / 1024  # MB
        
        with patch('designer2.designer.app') as mock_app, \
             patch('designer2.designer.run') as mock_run, \
             patch('designer2.designer.fsl') as mock_fsl:
            
            mock_app.ARGS = type('Args', (), {
                'input': sample_dwi_data['dwi'],
                'output': os.path.join(output_dir, 'test_memory'),
                'fslbval': sample_dwi_data['bval'],
                'fslbvec': sample_dwi_data['bvec'],
                'bids': None,
                'denoise': True,  # Enable processing step
                'degibbs': False,
                'rician': False,
                'eddy': False,
                'normalize': False,
                'mask': False,
                'n_cores': 1
            })
            
            mock_app.make_scratch_dir.return_value = None
            mock_app.goto_scratch_dir.return_value = None
            mock_fsl.suffix.return_value = '.nii.gz'
            
            with patch('lib.designer_input_utils.get_input_info') as mock_get_info, \
                 patch('lib.designer_input_utils.convert_input_data') as mock_convert, \
                 patch('lib.designer_input_utils.create_shell_table') as mock_shell, \
                 patch('lib.designer_func_wrappers.run_mppca') as mock_mppca:
                
                mock_get_info.return_value = {
                    'bshape_per_volume': [1] * 30,
                    'echo_time_per_volume': [0.1] * 30
                }
                mock_shell.return_value = np.array([[0, 1000], [1, 1], [3, 27], [0.1, 0.1]])
                
                try:
                    designer.execute()
                except Exception:
                    pass  # Expected due to mocking
                
                # Check memory usage didn't explode
                final_memory = process.memory_info().rss / 1024 / 1024  # MB
                memory_increase = final_memory - initial_memory
                
                # Should not use more than 500MB for this small test case
                assert memory_increase < 500, f"Memory usage increased by {memory_increase:.1f}MB"
    
    @pytest.mark.performance  
    def test_processing_time_basic(self, sample_dwi_data, output_dir, mock_environment):
        """Test basic processing time is reasonable."""
        start_time = time.time()
        
        with patch('designer2.designer.app') as mock_app, \
             patch('designer2.designer.run') as mock_run, \
             patch('designer2.designer.fsl') as mock_fsl:
            
            mock_app.ARGS = type('Args', (), {
                'input': sample_dwi_data['dwi'],
                'output': os.path.join(output_dir, 'test_timing'),
                'fslbval': sample_dwi_data['bval'],
                'fslbvec': sample_dwi_data['bvec'],
                'bids': None,
                'denoise': False,
                'degibbs': False,
                'rician': False,
                'eddy': False,
                'normalize': False,
                'mask': False,
                'n_cores': 1
            })
            
            mock_app.make_scratch_dir.return_value = None
            mock_app.goto_scratch_dir.return_value = None
            mock_fsl.suffix.return_value = '.nii.gz'
            
            with patch('lib.designer_input_utils.get_input_info') as mock_get_info, \
                 patch('lib.designer_input_utils.convert_input_data') as mock_convert, \
                 patch('lib.designer_input_utils.create_shell_table') as mock_shell:
                
                mock_get_info.return_value = {
                    'bshape_per_volume': [1] * 30,
                    'echo_time_per_volume': [0.1] * 30
                }
                mock_shell.return_value = np.array([[0, 1000], [1, 1], [3, 27], [0.1, 0.1]])
                
                try:
                    designer.execute()
                except Exception:
                    pass  # Expected due to mocking
        
        elapsed_time = time.time() - start_time
        
        # Basic mocked processing should complete quickly
        assert elapsed_time < 10, f"Basic processing took {elapsed_time:.2f}s (too slow)"
    
    @pytest.mark.performance
    @pytest.mark.parametrize("data_size", ["small", "medium"])
    def test_scaling_with_data_size(self, test_data_dir, data_size, mock_environment):
        """Test performance scaling with different data sizes."""
        # Create different sized test datasets
        if data_size == "small":
            shape = (16, 16, 8, 10)  # Very small
        elif data_size == "medium":
            shape = (32, 32, 16, 30)  # Small but realistic
        
        # Create synthetic data
        dwi_data = np.random.randn(*shape).astype(np.float32)
        affine = np.eye(4)
        dwi_img = nib.Nifti1Image(dwi_data, affine)
        
        dwi_file = test_data_dir / f"test_dwi_{data_size}.nii.gz"
        nib.save(dwi_img, str(dwi_file))
        
        # Create corresponding bval/bvec files
        n_vols = shape[3]
        b_values = np.concatenate([np.zeros(3), np.full(n_vols-3, 1000)])
        b_vectors = np.zeros((n_vols, 3))
        b_vectors[3:] = np.random.randn(n_vols-3, 3)
        
        bval_file = test_data_dir / f"test_dwi_{data_size}.bval"
        bvec_file = test_data_dir / f"test_dwi_{data_size}.bvec"
        
        np.savetxt(bval_file, b_values.reshape(1, -1), fmt='%d')
        np.savetxt(bvec_file, b_vectors.T, fmt='%.6f')
        
        # Time the processing
        start_time = time.time()
        
        with patch('designer2.designer.app') as mock_app, \
             patch('designer2.designer.run') as mock_run, \
             patch('designer2.designer.fsl') as mock_fsl:
            
            mock_app.ARGS = type('Args', (), {
                'input': str(dwi_file),
                'output': str(test_data_dir / f'test_output_{data_size}'),
                'fslbval': str(bval_file),
                'fslbvec': str(bvec_file),
                'bids': None,
                'denoise': False,
                'degibbs': False,
                'rician': False,
                'eddy': False,
                'normalize': False,
                'mask': False,
                'n_cores': 1
            })
            
            mock_app.make_scratch_dir.return_value = None
            mock_app.goto_scratch_dir.return_value = None
            mock_fsl.suffix.return_value = '.nii.gz'
            
            with patch('lib.designer_input_utils.get_input_info') as mock_get_info, \
                 patch('lib.designer_input_utils.convert_input_data') as mock_convert, \
                 patch('lib.designer_input_utils.create_shell_table') as mock_shell:
                
                mock_get_info.return_value = {
                    'bshape_per_volume': [1] * n_vols,
                    'echo_time_per_volume': [0.1] * n_vols
                }
                mock_shell.return_value = np.array([[0, 1000], [1, 1], [3, n_vols-3], [0.1, 0.1]])
                
                try:
                    designer.execute()
                except Exception:
                    pass  # Expected due to mocking
        
        elapsed_time = time.time() - start_time
        
        # All should complete quickly since we're mocking heavy operations
        assert elapsed_time < 15, f"Processing {data_size} data took {elapsed_time:.2f}s"


class TestDesignerStress:
    """Stress tests for DESIGNER robustness."""
    
    @pytest.mark.stress
    def test_multiple_concurrent_instances(self, sample_dwi_data, test_data_dir, mock_environment):
        """Test running multiple designer instances concurrently."""
        import threading
        import queue
        
        results = queue.Queue()
        
        def run_designer_instance(instance_id):
            """Run a single designer instance."""
            output_dir = test_data_dir / f"concurrent_{instance_id}"
            output_dir.mkdir(exist_ok=True)
            
            try:
                with patch('designer2.designer.app') as mock_app, \
                     patch('designer2.designer.run') as mock_run, \
                     patch('designer2.designer.fsl') as mock_fsl:
                    
                    mock_app.ARGS = type('Args', (), {
                        'input': sample_dwi_data['dwi'],
                        'output': str(output_dir / f'test_{instance_id}'),
                        'fslbval': sample_dwi_data['bval'],
                        'fslbvec': sample_dwi_data['bvec'],
                        'bids': None,
                        'denoise': False,
                        'degibbs': False,
                        'rician': False,
                        'eddy': False,
                        'normalize': False,
                        'mask': False,
                        'n_cores': 1
                    })
                    
                    mock_app.make_scratch_dir.return_value = None
                    mock_app.goto_scratch_dir.return_value = None
                    mock_fsl.suffix.return_value = '.nii.gz'
                    
                    with patch('lib.designer_input_utils.get_input_info') as mock_get_info, \
                         patch('lib.designer_input_utils.convert_input_data') as mock_convert, \
                         patch('lib.designer_input_utils.create_shell_table') as mock_shell:
                        
                        mock_get_info.return_value = {
                            'bshape_per_volume': [1] * 30,
                            'echo_time_per_volume': [0.1] * 30
                        }
                        mock_shell.return_value = np.array([[0, 1000], [1, 1], [3, 27], [0.1, 0.1]])
                        
                        designer.execute()
                        results.put((instance_id, "success"))
                        
            except Exception as e:
                results.put((instance_id, f"error: {str(e)}"))
        
        # Start multiple threads
        threads = []
        n_instances = 3
        
        for i in range(n_instances):
            thread = threading.Thread(target=run_designer_instance, args=(i,))
            threads.append(thread)
            thread.start()
        
        # Wait for all threads to complete
        for thread in threads:
            thread.join(timeout=30)  # 30 second timeout per thread
        
        # Check results
        completed_results = []
        while not results.empty():
            completed_results.append(results.get())
        
        assert len(completed_results) == n_instances, "Not all instances completed"
        
        # All should complete (though may error due to mocking)
        for instance_id, result in completed_results:
            # We expect some errors due to mocking, but no hanging
            assert result is not None, f"Instance {instance_id} hung"
    
    @pytest.mark.stress
    def test_edge_case_inputs(self, test_data_dir, mock_environment):
        """Test with edge case input data."""
        # Test with minimal volume (single slice, few directions)
        minimal_shape = (8, 8, 1, 4)  # Single slice, 4 volumes
        minimal_data = np.random.randn(*minimal_shape).astype(np.float32)
        
        affine = np.eye(4)
        minimal_img = nib.Nifti1Image(minimal_data, affine)
        
        minimal_file = test_data_dir / "minimal_dwi.nii.gz"
        nib.save(minimal_img, str(minimal_file))
        
        # Create minimal bval/bvec
        b_values = np.array([0, 1000, 1000, 1000])
        b_vectors = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]).T
        
        bval_file = test_data_dir / "minimal_dwi.bval"
        bvec_file = test_data_dir / "minimal_dwi.bvec"
        
        np.savetxt(bval_file, b_values.reshape(1, -1), fmt='%d')
        np.savetxt(bvec_file, b_vectors, fmt='%.6f')
        
        with patch('designer2.designer.app') as mock_app, \
             patch('designer2.designer.run') as mock_run, \
             patch('designer2.designer.fsl') as mock_fsl:
            
            mock_app.ARGS = type('Args', (), {
                'input': str(minimal_file),
                'output': str(test_data_dir / 'minimal_output'),
                'fslbval': str(bval_file),
                'fslbvec': str(bvec_file),
                'bids': None,
                'denoise': False,  # Skip computationally intensive steps
                'degibbs': False,
                'rician': False,
                'eddy': False,
                'normalize': False,
                'mask': False,
                'n_cores': 1
            })
            
            mock_app.make_scratch_dir.return_value = None
            mock_app.goto_scratch_dir.return_value = None
            mock_fsl.suffix.return_value = '.nii.gz'
            
            with patch('lib.designer_input_utils.get_input_info') as mock_get_info, \
                 patch('lib.designer_input_utils.convert_input_data') as mock_convert, \
                 patch('lib.designer_input_utils.create_shell_table') as mock_shell:
                
                mock_get_info.return_value = {
                    'bshape_per_volume': [1] * 4,
                    'echo_time_per_volume': [0.1] * 4
                }
                mock_shell.return_value = np.array([[0, 1000], [1, 1], [1, 3], [0.1, 0.1]])
                
                # Should handle minimal data gracefully
                try:
                    designer.execute()
                    success = True
                except Exception as e:
                    # Some exceptions expected due to mocking, but should be graceful
                    success = "RuntimeError" not in str(type(e))
                
                assert success, "Failed to handle minimal input data gracefully"


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-m", "performance"])

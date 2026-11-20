"""Ensure comparison metrics use the same physical voxels after MRtrix reorientation."""
import nibabel as nib
import numpy as np
import pytest

from tests.scripts.compare_rician_sh_maps import mask_on_grid


def test_mask_axis_flip_preserves_physical_selection():
    data=np.zeros((3,4,5),np.float32)
    data[0,1,2]=1
    mask=nib.Nifti1Image(data,np.diag([2.,2.,2.,1.]))
    flipped=mask.as_reoriented(np.array([[0,-1],[1,1],[2,1]]))
    reference=nib.Nifti1Image(np.zeros((3,4,5,12)),mask.affine)
    aligned=mask_on_grid(flipped,reference)
    np.testing.assert_array_equal(aligned.get_fdata(),data)
    np.testing.assert_allclose(aligned.affine,mask.affine)


def test_mask_axis_permutation_preserves_physical_selection():
    data=np.arange(60,dtype=np.float32).reshape(3,4,5)
    mask=nib.Nifti1Image(data,np.diag([2.,3.,4.,1.]))
    permuted=mask.as_reoriented(np.array([[1,-1],[2,1],[0,-1]]))
    aligned=mask_on_grid(permuted,mask)
    np.testing.assert_array_equal(aligned.get_fdata(),data)
    np.testing.assert_allclose(aligned.affine,mask.affine)


@pytest.mark.parametrize('mismatch',['translation','shape'])
def test_different_grid_is_rejected(mismatch):
    affine=np.eye(4)
    shape=(3,4,5)
    reference=nib.Nifti1Image(np.ones(shape),affine)
    if mismatch=='translation':
        affine=affine.copy();affine[0,3]=10
    else:
        shape=(3,4,6)
    with pytest.raises(ValueError,match='grid differs'):
        mask_on_grid(nib.Nifti1Image(np.ones(shape),affine),reference)

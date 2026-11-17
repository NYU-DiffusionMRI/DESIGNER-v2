"""Read gradient-coil tensors for use with scanner-frame diffusion gradients."""

import nibabel as nib
import numpy as np

GRADDEV_DESCRIPTIONS = {
    "tortoise": b"DESIGNER graddev (TORTOISE voxel axes)",
    "gradunwarp": b"DESIGNER graddev (gradunwarp FSL axes)",
}


def sanitize_graddev_matrix(mat):
    """Replace any 3x3 gradient-deviation block containing a NaN with the
    identity matrix as a whole block (not entrywise). Works around a
    confirmed TORTOISE bug: CreateGradientNonlinearityBMatrix emits NaN for
    the entire 3x3 block exactly at the isocenter voxel, where the true
    value is Identity. Accepts (...,9) or (...,3,3); returns the same shape.
    """
    flat = mat.reshape(*mat.shape[:-1], 3, 3) if mat.shape[-1] == 9 else mat
    bad = np.isnan(flat).any(axis=(-2, -1))
    if not bad.any():
        return mat
    flat = np.where(bad[..., None, None], np.eye(3, dtype=flat.dtype), flat)
    return flat.reshape(mat.shape)


def load_gnc_field(filename, dwi, gnc_format=None):
    """Return L on the DWI grid, with L @ g mapping nominal to actual gradients."""
    from mrtrix3 import MRtrixError, app, run

    if not filename.lower().endswith((".nii", ".nii.gz")):
        if gnc_format in ("gradunwarp", "tortoise"):
            raise MRtrixError("For a native graddev field, supply the original NIfTI with its voxel axes.")
        run.command(["mrconvert", filename, "L_gnc.nii", "-force"], show=False)
        filename = "L_gnc.nii"
        gnc_format = "scanner"

    field = nib.load(filename)
    if field.shape[3:] not in ((9,), (3, 3)):
        raise MRtrixError(f"-gnc field must be (X,Y,Z,9) or (X,Y,Z,3,3); got {field.shape}")
    if gnc_format is None:
        desc = field.header["descrip"].item()
        gnc_format = "gradunwarp" if desc == GRADDEV_DESCRIPTIONS["gradunwarp"] else "tortoise"

    # MRtrix may permute or reverse the DWI's spatial axes during import.
    # Arrange the scalar component volumes on that grid without interpolation.
    # Their component basis is still described by the original field affine.
    dwi_affine = dwi.transform.copy()
    dwi_affine[:3, :3] *= dwi.vox[:3]
    voxel_to_voxel = np.linalg.solve(dwi_affine, field.affine)
    aligned = field.as_reoriented(nib.orientations.io_orientation(voxel_to_voxel))
    if (aligned.shape[:3] != dwi.shape[:3]
            or not np.allclose(aligned.affine, dwi_affine, rtol=0, atol=1e-4)):
        raise MRtrixError("-gnc field must be on the same spatial grid as the input DWI.")

    blocks = np.asanyarray(aligned.dataobj).reshape(*aligned.shape[:3], 3, 3)
    app.console(f"Reading -gnc field in {gnc_format} convention")
    if gnc_format == "scanner":
        return blocks

    axes = field.affine[:3, :3]
    R = axes / np.linalg.norm(axes, axis=0)
    if not np.allclose(R.T @ R, np.eye(3), atol=1e-4):
        raise MRtrixError("-gnc native field must have orthogonal voxel axes.")
    if gnc_format == "gradunwarp" and np.linalg.det(R) > 0:
        R[:, 0] *= -1  # FSL uses a radiological basis.
    # J has one row per nominal coil and one column per spatial derivative.
    # Its transpose maps a nominal gradient to the actual gradient in native axes.
    J = sanitize_graddev_matrix(blocks)
    return R @ J.swapaxes(-1, -2) @ R.T

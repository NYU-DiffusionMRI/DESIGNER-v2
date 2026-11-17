"""Exercise gradunwarp/FSL axes through TMI with an analytic coil field."""

import shutil
import subprocess

import nibabel as nib
import numpy as np
import pytest

from lib.gnc.backends import derive_gnc_graddev


@pytest.mark.skipif(
    any(shutil.which(tool) is None for tool in ("gradient_unwarp.py", "calc_grad_perc_dev", "fslmerge")),
    reason="requires gradunwarp and FSL",
)
def test_gradunwarp_positive_determinant_fiber_fit(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("LOKY_MAX_CPU_COUNT", "1")
    # A single Siemens A(3,0) term gives dz = k*(z^3 - 1.5*z*(x^2+y^2)).
    # Its analytic derivatives supply a reference independent of the backend.
    coeff = tmp_path / "analytic.grad"
    coeff.write_text(" 0.2 m = R0\n  1 A( 3, 0)                 0.5                      z\n")
    c, s = np.cos(np.deg2rad(20.)), np.sin(np.deg2rad(20.))
    R = np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])
    affine = nib.affines.from_matvec(R @ np.diag([1., 1.5, 2.]), [35., 40., 40.])
    shape = (5, 5, 5)
    nib.save(nib.Nifti1Image(np.zeros(shape, np.float32), affine), "ref.nii")
    tagged = derive_gnc_graddev("ref.nii", coeff, "siemens", "gnc", "gradunwarp")
    native = nib.load(tagged)
    native.header["descrip"] = b""
    nib.save(native, "preproc_graddev.nii")

    indices = np.moveaxis(np.indices(shape), 0, -1)
    x, y, z = np.moveaxis(nib.affines.apply_affine(affine, indices), -1, 0)
    k = 0.5 / 200**2
    L = np.broadcast_to(np.eye(3), shape + (3, 3)).copy()
    L[..., 0, 2] = -3*k*x*z
    L[..., 1, 2] = -3*k*y*z
    L[..., 2, 2] += k*(3*z*z - 1.5*(x*x + y*y))
    dirs = np.random.default_rng(12).normal(size=(48, 3))
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    scanner_gradients = np.vstack([np.zeros((4, 3)), dirs])
    bvals = np.r_[np.zeros(4), np.full(48, 1000.)]
    # FSL reverses x on this positive-determinant grid.
    bvecs = (scanner_gradients @ R) * [-1., 1., 1.]
    np.savetxt("dwi.bval", bvals[None])
    np.savetxt("dwi.bvec", bvecs.T)
    actual = np.einsum("...ij,nj->...ni", L, scanner_gradients)
    fiber = np.array([1., 2., 3.]) / np.sqrt(14.)
    D = 0.35*np.eye(3) + 1.2*np.outer(fiber, fiber)
    attenuation = np.einsum("...ni,ij,...nj->...n", actual, D, actual)
    signal = 100*np.exp(-bvals/1000 * attenuation)
    nib.save(nib.Nifti1Image(signal, affine), "dwi.nii")

    command = ["tmi", "-DTI", "-n_cores", "1"]
    ambiguous = subprocess.run(command + ["-gnc", "preproc_graddev.nii", "dwi.nii", "ambiguous"],
                               capture_output=True, text=True)
    assert ambiguous.returncode != 0
    assert "Supply -gnc_format tortoise" in ambiguous.stderr + ambiguous.stdout

    for output, options in [
        ("tagged", ["-gnc", tagged]),
        ("explicit", ["-gnc", "preproc_graddev.nii", "-gnc_format", "gradunwarp"]),
        ("wrong_frame", ["-gnc", "preproc_graddev.nii", "-gnc_format", "tortoise"]),
    ]:
        subprocess.run(command + options + ["dwi.nii", output], check=True)
        fe = nib.load(f"{output}/fe_dti.nii").get_fdata()
        principal = fe / np.linalg.norm(fe, axis=-1, keepdims=True)
        angles = np.degrees(np.arccos(np.clip(principal @ fiber, 0., 1.)))
        if output == "wrong_frame":
            assert np.median(angles) > 0.5
        else:
            # FSL differentiates the cubic field by finite differences.
            assert np.max(angles) < 0.01
            md = nib.load(f"{output}/md_dti.nii").get_fdata()
            np.testing.assert_allclose(md, np.trace(D)/3, rtol=2e-4)

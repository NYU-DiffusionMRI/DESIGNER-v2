"""Regression for voxelwise b-matrices in AKC replacement and smoothing refits."""
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest


@pytest.mark.parametrize("echo_times", [(60,), (60, 90)])
def test_gnc_survives_tensor_refits(tmp_path, echo_times, monkeypatch):
    monkeypatch.setenv("LOKY_MAX_CPU_COUNT", "1")
    # Use the existing two-shell gradient table with a known diffusion signal.
    series = Path(__file__).parent / "data" / "D5" / "dwi1_ds"
    bvals = np.loadtxt(series.with_suffix(".bval"))
    bvecs = np.loadtxt(series.with_suffix(".bvec"))
    shape = (5, 5, 5)
    rng = np.random.default_rng(42)
    s0 = rng.uniform(90, 110, size=shape)
    diffusivity = rng.uniform(0.6, 0.9, size=shape)
    bD = diffusivity[..., None] * bvals / 1000
    signals = []
    for te in echo_times:
        kurtosis = np.full(shape, 0.7)
        # Exercise two outliers in the single-TE case, and one then zero at the two TEs.
        if te == 60:
            kurtosis[2, 2, 2] = -2
            if len(echo_times) == 1:
                kurtosis[2, 2, 3] = -2
        signals.append(s0[..., None] * np.exp(-bD + kurtosis[..., None] * bD**2 / 6))
    dwi = np.concatenate(signals, axis=-1)
    nib.save(nib.Nifti1Image(dwi, np.eye(4)), tmp_path / "dwi.nii")
    np.savetxt(tmp_path / "dwi.bval", [np.tile(bvals, len(echo_times))])
    np.savetxt(tmp_path / "dwi.bvec", np.tile(bvecs, (1, len(echo_times))))
    np.savetxt(tmp_path / "dwi.echotime", np.repeat(echo_times, bvals.size))

    # L = scale * I multiplies b by scale**2 at each voxel.
    scale = np.linspace(1.0, 1.5, num=np.prod(shape)).reshape(shape)
    field = scale[..., None] * np.eye(3).ravel()
    nib.save(nib.Nifti1Image(field, np.eye(4)), tmp_path / "L.nii")

    command = ["tmi", "-DTI", "-DKI", "-akc_outliers", "-fit_smoothing", "10",
               "-fit_constraints", "0,0,0", "-n_cores", "1"]
    for output, options in [("nominal", []), ("corrected", ["-gnc", "L.nii"])]:
        subprocess.run(command + options + ["dwi.nii", output], cwd=tmp_path, check=True)

    suffixes = [""] if len(echo_times) == 1 else [f"_te{te:.1f}" for te in echo_times]
    for suffix in suffixes:
        for model in ("dti", "dki"):
            filename = f"md_{model}{suffix}.nii"
            nominal = nib.load(tmp_path / "nominal" / filename).get_fdata()
            corrected = nib.load(tmp_path / "corrected" / filename).get_fdata()
            assert np.all(nominal > 0)
            np.testing.assert_allclose(corrected, nominal / scale**2, rtol=1e-6)

    if len(echo_times) == 1:
        # The AKC intermediate contains the six diffusion tensor elements first.
        nominal = nib.load(tmp_path / "nominal/dt_dki_akc.nii").get_fdata()[..., :6]
        corrected = nib.load(tmp_path / "corrected/dt_dki_akc.nii").get_fdata()[..., :6]
        np.testing.assert_allclose(corrected, nominal / scale[..., None]**2, atol=1e-7)

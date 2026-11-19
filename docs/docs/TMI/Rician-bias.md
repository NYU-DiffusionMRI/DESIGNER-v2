---
layout: default
title: Rician bias
has_children: false
parent: TMI
nav_order: 5
---

# Rician bias: correction and expectation fitting
{: .no_toc }

TMI offers three optional ways to account for a retained Rician magnitude
noise floor. The user chooses the estimator explicitly; it is not selected
automatically from the model or SNR.

## Contents
{: .no_toc .text-delta }

1. TOC
{:toc}

## Why magnitude signals are biased

If independent real and imaginary Gaussian noise components have standard
deviation `sigma`, the measured magnitude is

$$M=\sqrt{(s+n_1)^2+n_2^2},\qquad n_1,n_2\sim\mathcal N(0,\sigma^2),$$

where `s` is the underlying nonnegative signal. The mean magnitude
`mu_R(s,sigma) = E[M | s,sigma]` exceeds `s`, particularly at low directional
SNR. Averaging directions reduces random variation but does not remove this
mean offset. A high b0 SNR does not imply high SNR in every diffusion-weighted
direction.

Magnitude denoising can reduce the residual variance while retaining the
original Rician mean bias. This distinction matters for choosing the noise map
and for interpreting the expectation estimator below.

## Available estimators

| Option | Correction or observation model | DTI/DKI/WDKI | SH/SMI |
|---|---|---|---|
| No correction flag | Original input | Existing WLLS tensor pipeline | Ordinary SH LLS, without signal nonnegativity constraints |
| `-RicianClipMax` | `sqrt(max(0, M²−sigma²))` | Correct input, then retain the existing tensor pipeline | Corrected LLS with signal nonnegativity |
| `-RicianClipAbs` | `sqrt(abs(M²−sigma²))` | Correct input, then retain the existing tensor pipeline | Corrected LLS with signal nonnegativity |
| `-RicianExp` | Expected Rician magnitude fitted to original observations | Unconstrained signal-domain NLLS | Rician-mean NLLS with signal nonnegativity |

All three flags require `-sigma <noisemap>` and at least one of `-DTI`, `-DKI`,
`-WDKI`, `-SH` or `-SMI`. A sigma map by itself does not enable correction.
`-SH` saves rotational invariants without SMI model fitting; `-SMI` continues
from those invariants to Standard Model regression.

Choose **one** clipping flag. Specifying both `-RicianClipMax` and
`-RicianClipAbs` is an error. The former `-RicianClip` flag now reports a
migration message asking the user to choose explicitly. `-RicianExp` takes
precedence over either clipping flag for all selected models; specifying an
additional clipping flag is unnecessary and does not apply correction twice.

## Max-zero correction

`-RicianClipMax` subtracts **one sigma squared**, sets a negative squared
difference to zero, and takes the square root:

$$\widetilde M_{\rm max}=\sqrt{\max(0,M^2-\sigma^2)}.$$

This is an approximate signal correction. It is not an exact inverse of the
Rician mean and does not guarantee unbiased fitted parameters.

The distinction between fitting domains is important:

- **DTI/DKI:** the existing tensor fitter takes a logarithm. Corrected zeros
  use its existing small positive floor, which can create extreme log values
  and large parameter outliers.
- **SH/SMI:** SH LLS fits in signal space and does not take a logarithm. Zero
  targets are permitted. The fitted angular signal is constrained to be
  nonnegative; individual SH coefficients need not be nonnegative.

## Absolute-value correction

`-RicianClipAbs` uses

$$\widetilde M_{\rm abs}=\sqrt{|M^2-\sigma^2|}.$$

It is identical to max-zero when `M >= sigma`. Below sigma, it reflects the
negative squared difference upward. This avoids most corrected zeros before
tensor logarithms, but still gives zero at `M=sigma`, so the existing positive
log floor remains necessary.

This reflection changes the statistical correction: as `M` decreases from
sigma to zero, the corrected signal increases from zero to sigma. Thus the
formula is not monotonic below sigma and is not an exact Rician-bias inverse.
Improved log-domain stability does not by itself establish greater accuracy.
SH fitting has no logarithm, but users can still select this formula to compare
its statistical behavior with max-zero.

Both clipping formulas subtract one sigma squared. Neither is the second-moment
correction obtained from `E[M²] = s² + 2 sigma²`.

## Rician expectation fitting

`-RicianExp` models the expected magnitude instead of correcting each
measurement before the fit:

$$\widehat\theta=\arg\min_\theta\sum_i
\left[M_i-\mu_R(s_i(\theta),\sigma)\right]^2.$$

Every nonlinear search, including retries, starts from ordinary LLS on
**max-zero corrected measurements**. This initialization does not change if
`-RicianClipAbs` is also supplied. The objective always uses the original
observations. Occasional negative observations after denoising are preserved,
not folded upward with an absolute value.

For tensors, `s_i(theta) = exp((X theta)_i)` is the DTI/DKI signal model. The
fit has **no diffusion-eigenvalue, kurtosis, or monotonicity constraints** and
ignores `-fit_constraints` and `-akc_outliers`, with a log message. Positivity
of the exponential model is intrinsic to that signal parameterization.

For SH/SMI, `s=Yc` is constrained to be nonnegative at acquired directions and
on a deterministic 65,536-direction sphere grid. This constrains the fitted
signal, not the noisy observations. A finite grid does not guarantee
nonnegativity everywhere on the continuous sphere.

This is least squares with a mean-signal observation model, **not Rician
maximum likelihood**. It can describe Gaussian residuals around a retained
Rician mean after denoising, provided the supplied sigma describes that mean.
Correct mean modeling does not guarantee unbiased finite-sample tensor
parameters, rotational invariants or downstream SMI estimates. At sigma=0,
tensor fitting becomes signal-domain NLLS and SH fitting becomes constrained
LLS.

## Noise map and preprocessing

Supply a spatially aligned 3D map of the **original complex-component Gaussian
standard deviation**, in the same intensity units as the DWI. Values inside
the fitting mask must be finite and nonnegative. If denoising reduces residual
SD but leaves the original magnitude bias, use the original sigma responsible
for that bias, not the smaller residual SD. Match any intensity rescaling.

Use these options for data whose magnitude noise floor is still present.
Do not apply them on top of `designer -rician` or another correction that has
already removed that floor. Complex or phase-unwound denoising and other
preprocessing may change the mean/noise model and require different treatment.

## Compatibility and output conventions

| Setting or behavior | ClipMax / ClipAbs | RicianExp |
|---|---|---|
| Tensor fitting constraints | Retained | Ignored; unconstrained fit |
| `-akc_outliers` for tensors | Supported | Ignored |
| `-db0restore` for tensors | Supported, including together with AKC | Not yet supported; explicit error |
| Tensor fit smoothing, polynomial regression, external tensors | Existing behavior retained | Not yet supported; explicit error |
| GNC, shell selection, separate echo-time/LTE subsets | Supported | Supported |
| SH/SMI fitted-signal nonnegativity | Enforced | Enforced |
| SH/SMI b0 and isotropic-shell means | Raw arithmetic means retained | Raw arithmetic means retained |

For tensor fits, either clipping formula is applied once, including to b0
measurements, before the existing fitting/refinement pipeline. Restored b0
samples are not corrected a second time when saving trace images. The trace
output with `-RicianExp` remains measurement-derived, not an
expectation-corrected model prediction.

For SH/SMI, corrections apply to directional-shell SH fits. The existing raw
b0/isotropic means, SH normalization, cross-shell rank-one processing of
higher-order coefficients, and downstream SMI regression remain unchanged.
Nonnegativity is enforced before rank-one processing; S0 is not rank-one
processed.

The execution log identifies the selected estimator/formula. Standalone
`SH_shells.json` records the estimator and shell orders. Tensor expectation
fits also save convergence diagnostics. See the [tensor details]({% link docs/TMI/DTI-DKI.md %}#rician-correction-options)
and [SH diagnostics]({% link docs/TMI/SMI.md %}#rician-bias-expectation-correction-with-a-nonnegative-signal).

## Command examples

Both clipping choices work with DKI b0 restoration and AKC replacement:

```bash
tmi dwi.nii dki_max -DKI -RicianClipMax -sigma sigma.nii -mask mask.nii -db0restore -akc_outliers
tmi dwi.nii dki_abs -DKI -RicianClipAbs -sigma sigma.nii -mask mask.nii -db0restore -akc_outliers
```

Either formula can also be used for standalone SH or SMI:

```bash
tmi dwi.nii sh_max -SH -RicianClipMax -sigma sigma.nii -mask mask.nii
tmi dwi.nii sh_abs -SH -RicianClipAbs -sigma sigma.nii -mask mask.nii
tmi dwi.nii smi_abs -SMI -RicianClipAbs -sigma sigma.nii -mask mask.nii
```

A combined run uses the same selected formula for every requested model:

```bash
tmi dwi.nii combined_abs -DTI -DKI -SH -RicianClipAbs -sigma sigma.nii
tmi dwi.nii combined_exp -DTI -DKI -SH -RicianExp -sigma sigma.nii
```

To choose different formulas for tensors and SH/SMI, use separate calls:

```bash
tmi dwi.nii tensor_abs -DTI -DKI -RicianClipAbs -sigma sigma.nii
tmi dwi.nii sh_max -SH -RicianClipMax -sigma sigma.nii
```

Omit correction flags to retain the existing estimators. For independent
standalone SH orders, add `-lmaxSH`, for example `0,2,4` for an input containing
b=0,1000,2000 at one b-shape and echo time.

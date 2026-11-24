---
layout: default
title: Rician bias
has_children: false
parent: TMI
nav_order: 5
---

# Rician bias
{: .no_toc }

Magnitude noise raises the mean signal, particularly in low-SNR directions:

$$
\begin{gathered}
M=\sqrt{(s+n_1)^2+n_2^2},\qquad
n_1,n_2\overset{\mathrm{iid}}{\sim}\mathcal N(0,\sigma^2),\\[6pt]
\mu_R(s,\sigma)=\mathbb E[M\mid s,\sigma]>s\quad(\sigma>0).
\end{gathered}
$$

Here $s$ is the underlying signal and $\sigma$ is the standard deviation of
each complex Gaussian noise component. Averaging directions reduces variance,
but retains this bias; high b0 SNR does not guarantee high directional SNR.

## Recommended correction

- **DTI/DKI/WDKI:** `-RicianClipAbs -sigma sigma.nii`, followed by weighted linear least squares (WLLS).
- **SH/SMI with data at $b\geq3000\,\mathrm{s/mm^2}$:** `-RicianExp -sigma sigma.nii`, fitting the expected Rician magnitude. Add `-nonNeg` only if a nonnegative SH signal is desired.
- **SH/SMI with all shells below $3000\,\mathrm{s/mm^2}$:** `-RicianClipAbs -sigma sigma.nii`, followed by ordinary SH LLS.

One estimator applies to all shells in a call. For mixed-shell data including
$b\geq3000\,\mathrm{s/mm^2}$, use `-RicianExp` for the whole SH/SMI fit;
TMI does not automatically switch estimators at this threshold.

## Choose an estimator

All correction options apply to `-DTI`, `-DKI`, `-WDKI`, `-SH` and `-SMI`.

| Option | DTI/DKI/WDKI | SH/SMI |
|---|---|---|
| No correction flag | WLLS | Ordinary LLS |
| `-RicianClipMax` | Max-zero correction, then WLLS | Ordinary corrected LLS |
| `-RicianClipAbs` | Absolute-value correction, then WLLS | Ordinary corrected LLS |
| `-RicianExp` | Unconstrained Rician-mean NLLS | Unconstrained Rician-mean NLLS |

Choose one clipping flag. `-RicianExp` overrides either clipping flag for all
selected models. Use separate calls to choose different methods for tensors
and SH/SMI. Supplying `-sigma` alone does not enable correction.

## Clipping corrections

$$
\begin{aligned}
\texttt{-RicianClipMax}:\quad
\widetilde M_{\mathrm{max}}&=\sqrt{\max(0,M^2-\sigma^2)},\\[6pt]
\texttt{-RicianClipAbs}:\quad
\widetilde M_{\mathrm{abs}}&=\sqrt{\lvert M^2-\sigma^2\rvert}.
\end{aligned}
$$

Both subtract one $\sigma^2$ and agree when $M\geq\sigma$. Below this threshold,
max-zero returns zero; absolute value reflects the squared difference upward.
These are approximate corrections, not exact inverses of the Rician mean.

For tensors, corrected zeros use a small positive floor before logarithms and
can produce outliers. Absolute value avoids most zeros, but remains zero at
$M=\sigma$ and is nonmonotonic below it. SH fits in signal space, so zero
measurements do not create this logarithm problem.

## Rician expectation fitting

For $\sigma>0$, the expected magnitude is

$$
\mu_R(s,\sigma)=\sigma\sqrt{\frac{\pi}{2}}\,e^{-z}
\left[(1+2z)I_0(z)+2zI_1(z)\right],
\qquad z=\frac{s^2}{4\sigma^2}.
$$

Here $I_0$ and $I_1$ are modified Bessel functions of the first kind, of
orders zero and one. At zero noise, $\mu_R(s,0)=s$ for $s\geq0$.

`-RicianExp` fits the original observations using the expected magnitude:

$$
\widehat\theta=\arg\min_\theta\sum_i
\left[M_i-\mu_R\!\left(s_i(\theta),\sigma\right)\right]^2.
$$

Initialization uses ordinary LLS on **absolute-value corrected measurements for
tensors**, and **max-zero corrected measurements for SH/SMI**, including retries.
Negative observations after denoising are preserved.

For tensors, $s_i(\theta)=\exp((X\theta)_i)$; diffusion and kurtosis parameters
are unconstrained. For SH/SMI, $s=Yc$ is unconstrained by default and the
Rician mean depends on $\lvert Yc\rvert$. Only `-RicianExp -nonNeg` enforces

$$Y_{\mathrm{grid}}c\geq0$$

at acquired directions and on a 65,536-direction sphere grid, before cross-shell
rank-one processing. This constrains the signal, not individual coefficients;
a finite grid is not a continuous-sphere guarantee. A positive physical signal
can have a negative truncated SH projection, so this optional constraint can
introduce truncation bias. Nonnegativity can also substantially reduce the
variance of `RicianExp` estimates: this bias–variance tradeoff is left to the
user. `-nonNeg` is not available with either clipping flag.

Expectation fitting models the mean rather than the Rician likelihood. None
of these methods guarantees unbiased finite-sample parameters.

## Noise map

All three options require `-sigma <noisemap>`: an aligned 3D map with finite,
nonnegative values inside the mask, in matching signal units. Use the
**original $\sigma$ responsible for the retained Rician bias**, even if
denoising has reduced residual variance. Do not correct data whose noise floor
has already been removed, for example with `designer -rician`.

## Compatibility

| Setting | ClipMax / ClipAbs | RicianExp |
|---|---|---|
| Tensor constraints and `-akc_outliers` | Supported | Not yet supported |
| `-db0restore`, including with AKC | Supported | Not yet supported |
| `-fit_smoothing`, `-polyreg`, external `-dt` | Supported | Not yet supported |
| GNC, shell selection, multi-TE subsets | Supported | Supported |

With tensor `-RicianExp`, constraints and AKC replacement are not applied;
requesting them produces a warning. Other unsupported combinations produce an error.
Tensor clipping includes b0 and is applied once. SH/SMI also correct b0 and
isotropic (STE) shells:
clipping precedes averaging; expectation fitting inverts $\mu_R$ at the sample
mean using the original $\sigma$, with zero amplitude below the Rayleigh floor.
Higher-order invariants are zero for these shells. Rank-one processing and SMI
regression are unchanged. In combined corrected tensor/SH runs, b0 restoration
applies to tensors; SH/SMI use the original observations.

## Examples

The SH/SMI expectation examples below are recommended when the data include
$b\geq3000\,\mathrm{s/mm^2}$. For data entirely below this threshold, use the
absolute-value correction examples instead.

```bash
# DTI + DKI with absolute-value correction, followed by WLLS
tmi dwi.nii tensor_abs -DTI -DKI -RicianClipAbs -sigma sigma.nii

# DTI + DKI with absolute-value correction and b0 restoration
tmi dwi.nii tensor_abs_restore -DTI -DKI -RicianClipAbs -sigma sigma.nii -db0restore

# SH/SMI with all shells below 3000 s/mm^2 (choose the desired output)
tmi dwi.nii sh_abs -SH -RicianClipAbs -sigma sigma.nii
tmi dwi.nii smi_abs -SMI -RicianClipAbs -sigma sigma.nii

# SH invariants only, using Rician expectation
tmi dwi.nii sh_exp -SH -RicianExp -sigma sigma.nii

# SMI using Rician expectation
tmi dwi.nii smi_exp -SMI -RicianExp -sigma sigma.nii

# SMI + separate SH maps with nonnegative fitted SH signals (two calls)
tmi dwi.nii smi_positive -SMI -RicianExp -sigma sigma.nii -nonNeg
tmi dwi.nii sh_positive -SH -RicianExp -sigma sigma.nii -nonNeg
```

`-SMI` already saves invariants in `rotinvs_smi.nii`; the extra `-SH` call saves
separate `S0.nii`, `S2.nii`, … files. The two flags cannot be combined in one call.
Add `-mask brain_mask.nii` to restrict fitting to the brain. Add `-lmaxSH`
to set standalone SH orders independently, for example `0,2,4` for
b=0,1000,2000 at one b-shape and echo time. See [tensor details]({% link docs/TMI/DTI-DKI.md %}#rician-correction-options)
and [SH details]({% link docs/TMI/SMI.md %}#rician-bias-expectation-correction)
for solver diagnostics.

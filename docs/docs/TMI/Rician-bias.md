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

## Choose an estimator

All correction options apply to `-DTI`, `-DKI`, `-WDKI`, `-SH` and `-SMI`.

| Option | DTI/DKI/WDKI | SH/SMI |
|---|---|---|
| No correction flag | Existing WLLS pipeline | Ordinary LLS |
| `-RicianClipMax` | Max-zero correction, then existing tensor fit | Corrected LLS with signal nonnegativity |
| `-RicianClipAbs` | Absolute-value correction, then existing tensor fit | Corrected LLS with signal nonnegativity |
| `-RicianExp` | Unconstrained Rician-mean NLLS | Rician-mean NLLS with signal nonnegativity |

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

`-RicianExp` fits the original observations using the expected magnitude:

$$
\widehat\theta=\arg\min_\theta\sum_i
\left[M_i-\mu_R\!\left(s_i(\theta),\sigma\right)\right]^2.
$$

Every search, including retries, starts from **ordinary LLS on max-zero
corrected measurements**. This remains true if `-RicianClipAbs` is also given.
Negative observations after denoising are preserved.

For tensors, $s_i(\theta)=\exp((X\theta)_i)$; diffusion and kurtosis parameters
are unconstrained. For SH/SMI, $s=Yc$ and all corrected estimators enforce

$$Y_{\mathrm{grid}}c\geq0$$

at acquired directions and on a 65,536-direction sphere grid, before cross-shell
rank-one processing. This constrains the signal, not individual coefficients;
a finite grid is not a continuous-sphere guarantee.

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
| Tensor constraints and `-akc_outliers` | Retained | Ignored |
| `-db0restore`, including with AKC | Supported | Not yet supported |
| `-fit_smoothing`, `-polyreg`, external `-dt` | Supported | Not yet supported |
| GNC, shell selection, multi-TE subsets | Supported | Supported |

Unsupported expectation combinations produce an error. Tensor clipping includes
b0 and is applied once. SH/SMI retain raw b0/isotropic-shell means, SH
normalization, rank-one processing and downstream SMI regression.

## Examples

```bash
# Invariants only; choose either clipping formula
tmi dwi.nii sh_max -SH -RicianClipMax -sigma sigma.nii
tmi dwi.nii sh_abs -SH -RicianClipAbs -sigma sigma.nii

# Expected Rician signal for tensors and SH
tmi dwi.nii combined -DTI -DKI -SH -RicianExp -sigma sigma.nii

# Either clipping option supports b0 restoration and AKC together
tmi dwi.nii dki_abs -DKI -RicianClipAbs -sigma sigma.nii -db0restore -akc_outliers
```

Replace `-SH` with `-SMI` to include Standard Model fitting. Add `-lmaxSH`
to set standalone SH orders independently, for example `0,2,4` for
b=0,1000,2000 at one b-shape and echo time. See [tensor details]({% link docs/TMI/DTI-DKI.md %}#rician-correction-options)
and [SH details]({% link docs/TMI/SMI.md %}#rician-bias-expectation-correction-with-a-nonnegative-signal)
for solver diagnostics.

---
layout: default
title: SMI
has_children: false
parent: TMI
nav_order: 3
---

# Standard Model Imaging
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Standard Model Imaging (SMI)

In addition to the python implementation as a part of TMI, SMI is offered as a standalone matlab toolbox. [This MATLAB toolbox](https://github.com/NYU-DiffusionMRI/SMI/blob/master/SMI.m) contains all necessary functions for parameter estimation of the Standard Model (SM) of diffusion in white matter. Check [our recent paper](https://arxiv.org/pdf/2202.02399.pdf) for details on this implementation and on the Standard Model in general. Below we provide instructions on how to run the toolbox. See the [example.m](https://github.com/NYU-DiffusionMRI/SMI/blob/master/example.m) script that performs the parameter estimation in an [example dataset](https://cai2r.net/resources/standard-model-of-diffusion-in-white-matter-the-smi-toolbox/).

TMI offers two output modes. `-SMI` estimates signal rotational invariants and
then fits the Standard Model parameters. `-SH` estimates only the signal
rotational invariants, without SMI priors, training, or model regression.
Both modes use linear least squares for SH fitting by default. For Rician bias
correction, **we recommend `-RicianExp -sigma sigma.nii` when the data include
$b\geq3000\,\mathrm{s/mm^2}$**, fitting the expected Rician magnitude. **For data
entirely below this threshold, use `-RicianClipAbs -sigma sigma.nii`.**
Rician expectation fitting is slower because it requires nonlinear optimization;
below this threshold, its additional benefit generally does not justify the
extra runtime for typical datasets.
The selected estimator applies to all shells in the call, including in mixed-shell data.
**SH fits are unconstrained by default**; add
`-nonNeg` with `-RicianExp` to constrain the underlying SH signal.
Alternatively, choose `-RicianClipMax` or `-RicianClipAbs`, with `-sigma sigma.nii`,
to fit corrected measurements with ordinary LLS.
`-RicianExp` takes precedence over either clipping choice. See
[Rician bias]({% link docs/TMI/Rician-bias.md %}) for all alternatives and their assumptions.

### Choosing an estimation mode

Recommended examples for data with a retained Rician noise floor:

| Requested output | Example command |
|---|---|
| SH invariants only, all shells below $3000\,\mathrm{s/mm^2}$ | `tmi dwi.nii sh_abs -SH -RicianClipAbs -sigma sigma.nii` |
| SMI maps and SH invariants, all shells below $3000\,\mathrm{s/mm^2}$ | `tmi dwi.nii smi_abs -SMI -RicianClipAbs -sigma sigma.nii` |
| SH invariants only, including $b\geq3000\,\mathrm{s/mm^2}$ | `tmi dwi.nii sh_exp -SH -RicianExp -sigma sigma.nii` |
| SMI maps and SH invariants, including $b\geq3000\,\mathrm{s/mm^2}$ | `tmi dwi.nii smi_exp -SMI -RicianExp -sigma sigma.nii` |
| SMI maps and SH invariants, with nonnegative SH signals | `tmi dwi.nii smi_positive -SMI -RicianExp -sigma sigma.nii -nonNeg` |
| Separate SH invariant files with the same nonnegativity option | `tmi dwi.nii sh_positive -SH -RicianExp -sigma sigma.nii -nonNeg` |

`-SMI` saves invariants in `rotinvs_smi.nii`. Run the last two commands separately
for SMI maps plus standalone `S0.nii`, `S2.nii`, … files; do not combine `-SMI` and `-SH`.

Add `-mask brain_mask.nii` to restrict fitting to a brain mask. Choose either
`-SH` or `-SMI`. In standalone `-SH` mode, an explicit `-lmaxSH` sets the actual
SH fitting order and requires enough independent directions in each directional
shell. A comma-separated list assigns an order to each shell, for example
`-lmaxSH 0,2,4,4,6`; see [per-shell fitting orders](#choose-a-fitting-order-for-each-shell).
Omit it for automatic per-shell orders. Per-shell lists are supported only
with `-SH`: with `-SMI`, `-lmax` is a single invariant order for model regression,
and SH fitting orders remain automatic. Either correction option requires an
aligned noise map in the same signal units as the DWI and data whose Rician
noise floor has not already been removed; see [parameter estimation](#parameter-estimation).

---

## Overview: The Standard Model of diffusion in white matter

Over the last 15-20 years, multiple approaches aimed to model the physics of water diffusion in white matter have relied on similar assumptions. This led to the unifying framework dubbed Standard Model (SM) of diffusion in WM as formulated in ([Novikov et al., 2018](https://doi.org/10.1016/j.neuroimage.2018.03.006)). In a nutshell, this model disentangles signal contributions from different structures, i.e. compartments, present in a white matter voxel. 

<img width="1657" alt="kernel_wEqConvolution_v2" src="https://user-images.githubusercontent.com/54751227/152564788-fc6a0fe0-1002-4354-b75e-3f962303a9ad.png">

Briefly, axons (and possibly glial processes) are represented by impermeable zero-radius cylinders (the so-called “sticks”) arranged in locally coherent fiber fascicles. The diffusion in the extra-axonal space of each fascicle is assumed to be Gaussian and described by an axially symmetric diffusion tensor. The third, optional tissue compartment is the cerebro-spinal fluid. Such multicomponent fascicles (also called kernel) are distributed in a voxel according to an arbitrary fiber orientation distribution function (ODF). All fascicles in a voxel are assumed to have the same compartment fractions and diffusivities, and differ from each other only by orientation.

The SM encompasses a number of WM models made of anisotropic Gaussian compartments with axons represented by sticks ([Kroenke et al., 2004]( https://doi.org/10.1002/mrm.20260); [Jespersen et al., 2007](https://doi.org/10.1016/j.neuroimage.2006.10.037), [2010](https://doi.org/10.1016/j.neuroimage.2009.08.053); [Fieremans et al., 2011](https://doi.org/10.1016/j.neuroimage.2011.06.006); [Zhang et al., 2012](https://doi.org/10.1016/j.neuroimage.2012.03.072); [Sotiropoulos et al., 2012](https://doi.org/10.1016/j.neuroimage.2012.01.056); [Jensen et al., 2016](https://doi.org/10.1016/j.neuroimage.2015.09.049); [Jelescu et al., 2016a](https://doi.org/10.1002/nbm.3450); [Kaden et al., 2016](https://doi.org/10.1016/j.neuroimage.2016.06.002); [Reisert et al., 2017](https://doi.org/10.1016/j.neuroimage.2016.09.058); [Novikov et al., 2018](https://doi.org/10.1016/j.neuroimage.2018.03.006); [Veraart et al., 2018](https://doi.org/10.1016/j.neuroimage.2017.09.030), to mention a few). From the SM point of view, earlier models impose constraints either on compartment parameters or the functional form of the fiber ODF; such constraints improve robustness but may introduce biases into the estimation of remaining parameters.

For more details please look at our recent publication: [Reproducibility of the Standard Model of diffusion in white matter on clinical MRI systems, Volume 257, 15 August 2022, 119290, (2022), NeuroImage](https://www.sciencedirect.com/science/article/pii/S1053811922004104), or Section 3 in this review by [Novikov et al. (2019)](https://doi.org/10.1002/nbm.3998).

### Assumptions in a nutshell
- Sufficient coarse-graining for Gaussian diffusion at the fiber bundle level (thus, no time-dependence)
- All fiber bundles inside a voxel have the same microstructural properties
- Negligible water exchange between compartments in a bundle
- No assumptions on the fibre bundle's diffusivities
- No assumptions on the functional form of the ODF

{: .note }
_Note that this model does not apply to gray matter, where exchange and the presence of soma cannot be ignored._

---

## SMI input data

### Minimal requirements for linear tensor encoding (LTE)
- A noise map, this must be a 3D array of size [nx x ny x nz]. This can be previously estimated using MPPCA, see [this Matlab implementation](https://github.com/NYU-DiffusionMRI/mppca_denoise). Or by using the `sigma` output from the designer pipeline.

{: .note }
_On the b-value units:_ We recommend microstructural units [milliseconds / (squared micrometers)].
Note that typical scanner units are b=1000 [seconds / (squared millimeters)] which equal 1 [milliseconds / (squared micrometers)], thus to convert from scanner units to microstructural units one should divide by 1000. This is done automatically as part of the TMI toolbox.

### Data with non-LTE encodings and variable echo times (TE)
Input data can also have:
- Multiple **B**-tensor shapes (specified by β). This input must be a [1 x N] vector. Only axially symmetric b-tensors are supported. β is a unitless scalar between -0.5 and 1 that indicates the **B**-tensor shape (see figure below).
- Multiple echo times (TE, assumed to be in [milliseconds]). This input must be a [1 x N] vector.

In this general scenario, each measurement is thus fully specified by: a b-value (b), a unit direction (**u**) (axis of symmetry of **B**), a b-tensor shape (β), and TE. See the figure and equation below to understand how these parameters make a b-tensor **B**:
<p align="center">
  <img width="500" alt=" AxSymB_wEq" src="https://user-images.githubusercontent.com/54751227/152437987-d79193d1-1ecc-4707-bdc3-f7cd2dec6ad6.png">
</p>

  - If no β is supplied the code assumes β=1 (linear tensor encoding, LTE, measurements).
  - If no TE is supplied, the code assumes that the TE is the same across all measurements. In this case compartmental T2 values will not be outputted and water fractions will be T2-weighted.

---

## SMI outputs
Standard Model parameters (diffusion) + compartmental T2 values (only if multiple TE data was provided). See figure below for some examples:
<img width="1690" alt="SM_maps_github" src="https://user-images.githubusercontent.com/54751227/152456997-5f24f886-03f9-4eb2-a5f8-1dd767134eae.png">
- Fractions (f, fw) and anisotropy (p2, p4) are adimensional.
- Diffusivities are returned in microstructural units [squared micrometers / milliseconds] (independently of input units).
- Compartmental T2 values are returned in [milliseconds].
- Note that compartmental T2 maps will only be outputted if variable TE data was used.
- Rotational invariants for each diffusion shell will also be an output. To normalize them and remove T2 contrast, simply divide by b0 for fixed TE data or by proton density for variable TE data.

## Signal rotational invariants without SMI

To obtain rotational invariants alone, use `-SH` instead of `-SMI`.
**`-SH` alone uses ordinary linear least squares (LLS)** for the SH fit,
with automatic fitting orders and no noise map required:

```bash
tmi dwi.nii sh_output -SH
```

**`-SH -RicianExp` fits the expected Rician magnitude without a signal
constraint by default.** Supply a noise map with `-sigma`:

```bash
tmi dwi.nii sh_rician -SH -RicianExp -sigma sigma.nii
```

Add `-nonNeg` to constrain the SH signal at acquired directions and on a dense
sphere grid. This option requires `-RicianExp`. Supplying `-sigma` alone does not
switch from LLS to Rician fitting. The same estimator choice applies to `-SMI`.
See [Rician bias expectation correction](#rician-bias-expectation-correction)
for the observation model and initialization.

For clipped correction followed by ordinary LLS, use:

```bash
tmi dwi.nii sh_clip -SH -RicianClipMax -sigma sigma.nii
tmi dwi.nii sh_abs -SH -RicianClipAbs -sigma sigma.nii
```

All four commands bypass SMI priors and regression and save separate `S0.nii`,
`S2.nii`, and higher even-order files through the maximum fitted order. Each
contains one volume per shell, identified in `SH_shells.json`. Add `-lmaxSH 6`
to fit through order 6 and save `S0.nii`, `S2.nii`, `S4.nii`, and `S6.nii`.
In `-SH` mode, `-lmaxSH` sets the SH fitting order itself. See the
[standalone SH options]({% link docs/TMI/usage.md %}#-sh).

### Choose a fitting order for each shell

Orders 2, 4, and 6 require at least 6, 15, and 28 independent directional
measurements, respectively. With `-SH`, **each shell's maximum SH order can be
specified independently** using a comma-separated list. This works with
LLS, RicianExp, RicianClipMax and RicianClipAbs, for example:

```bash
# Input containing b=0, 1000, 2000, 7000 and 25000 s/mm², at one b-shape and TE:
tmi dwi.nii sh_lls -SH -lmaxSH 0,2,4,4,6
tmi dwi.nii sh_rician -SH -RicianExp -sigma sigma.nii -lmaxSH 0,2,4,4,6
```

| Shell b-value (s/mm²) | Maximum SH fitting order | Fitted invariants |
|---|---|---|
| 0 | 0 | S0 |
| 1000 | 2 | S0, S2 |
| 2000 | 4 | S0, S2, S4 |
| 7000 | 4 | S0, S2, S4 |
| 25000 | 6 | S0, S2, S4, S6 |

The list must contain one nonnegative even order for **every shell**, including
b=0 and isotropic shells, in ascending b-value, then b-shape, then echo time.
`SH_shells.json` records this order and the order actually fitted. Use zero for
b=0/isotropic shells; these use scalar estimation with the selected correction. A single
integer applies the same order to all directional shells; omitting `-lmaxSH`
selects fitting orders automatically. Write the list as one argument, with
no spaces or trailing comma. If adapting an earlier standalone SH command,
replace `-lmax` with `-lmaxSH`; `-SH -lmax` is rejected with a migration message.

The example saves `S0.nii`, `S2.nii`, `S4.nii`, and `S6.nii`, each with five
shell volumes. Volumes above a shell's fitting order are zero placeholders,
not estimated invariants: for example, b=1000 has no fitted S4 or S6, and only
b=25000 has a fitted S6. Fitting the shells together retains cross-shell
rank-one processing: each even order uses only the shells that support it.
Changing one shell's order can therefore also change higher-order invariants
in other shells through this joint step.

`-lmaxSH` is supported only with `-SH`. For `-SMI`, the separate `-lmax` option accepts one
value (0, 2, 4, or 6) specifying the invariant order used for Standard Model
regression; the SH fitting orders remain automatic per shell. The choice
between ordinary LLS, corrected LLS and Rician-expectation fitting applies to
both output modes. Only Rician-expectation fitting supports optional `-nonNeg`.

A requested order that cannot be identified from the directions is
rejected; an unsupported higher-order invariant should not be interpreted as
a measured zero. A full-rank design can still be poorly conditioned when
directions are nearly redundant; reducing the fitting order can improve stability.

The saved invariants are in the original signal units. For figures, normalized
invariants can be computed as `S_l / mean_b0`, using the mean b=0 image at the
same echo time. When comparing estimators, use the **same b=0 reference image**
for all methods so differences in the denominator do not confound the comparison.
Apply this normalization after fitting; the Rician fit itself uses the DWI and
noise map in their original matching units.

---

## Parameter estimation
First, the rotational invariants of the diffusion signal are estimated. By default, SH coefficients are fit with ordinary linear least squares. Supplying a noise map alone does not change this estimator.

### Rician bias expectation correction

Both `-SH -RicianExp` and `-SMI -RicianExp` fit

$$
\widehat c=\arg\min_c\sum_i
\left[\mu_R((Yc)_i,\sigma)-M_i\right]^2.
$$

For the conditional least-squares framework behind this expectation fit, see
[Veraart et al. (2013)](https://doi.org/10.1002/mrm.24529).

Every nonlinear search starts from ordinary LLS on
$\sqrt{\max(0,M^2-\sigma^2)}$, including retries. The objective uses original
observations, preserving negative values after denoising. The Rician mean
is even in $Yc$, so unconstrained fits can have signed lobes.

**Nonnegativity is optional:** add `-nonNeg` to enforce
$Y_{\mathrm{grid}}c\geq0$ at acquired directions and on a 65,536-direction sphere
grid, before cross-shell rank-one processing. This constrains signals, not
individual coefficients. A finite grid does not guarantee positivity everywhere;
constraining a truncated expansion can introduce bias even for a positive
physical signal. Nonnegativity can substantially reduce variance, so users
choose this bias–variance tradeoff. At $\sigma=0$, fits reduce to ordinary LLS,
or constrained LLS with `-nonNeg`.

The solver uses an analytic Jacobian and batched fitting. With `-nonNeg`,
violated grid directions trigger constrained refinement; failures return the
best grid-feasible estimate with a warning. Initialization remains max-zero LLS.

### Clipped correction

Choose `-RicianClipMax` or `-RicianClipAbs` to fit ordinary LLS to
$\sqrt{\max(0,M^2-\sigma^2)}$ or $\sqrt{\lvert M^2-\sigma^2\rvert}$, respectively.
These modes have **no signal constraint** and cannot be combined with `-nonNeg`.
SH fits in signal space, so corrected zeros need no logarithm floor.
`-RicianExp` overrides either clipping flag and keeps max-zero initialization.

All corrections require an aligned `-sigma` map in matching signal units. Use
the original noise level responsible for the retained magnitude bias, even
after denoising. They also apply to b0 and isotropic (STE) shells: clipping
precedes averaging; expectation fitting inverts the Rician mean of the sample
average using the original sigma. Higher-order invariants of these shells are
zero. Basis normalization, rank-one processing and SMI regression are unchanged;
$S_0$ is not rank-one processed. See [Rician bias]({% link docs/TMI/Rician-bias.md %})
for equations and assumptions. No method guarantees unbiased finite-sample estimates.

### Iteration limits and diagnostics

With `-RicianExp -nonNeg`, `-RicianSH_maxiter N` caps cumulative constrained-refinement iterations per
voxel and shell, including restarts. The initial batched fit retains its
standard limit. At the cap, the best grid-feasible fit is returned with a
warning; very small caps can change estimates substantially.

When a cap is supplied, `-SH` and `-SMI` save:

- `RicianSH_early_stop.nii`: refinement reached the cap.
- `RicianSH_iterations.nii`: cumulative refinement iterations.
- `RicianSH_nonconverged.nii`: cap exhaustion or other refinement failure.
- `RicianSH_diagnostics.json`: shell metadata and counts.

Each image volume describes one shell before rank-one processing. Zero
iterations means no refinement was needed, including b0/isotropic shells.
These diagnostics do not describe the initial batched fit's convergence.

### The same correction choices for DTI and DKI

All three flags also apply to DTI/DKI/WDKI. Both clipping choices retain
existing tensor settings, including `-db0restore -akc_outliers` together.
Tensor expectation fitting is unconstrained, ignores AKC and does not yet
support b0 restoration or spatial refinements. Tensor expectation starts from
absolute-value corrected LLS; SH/SMI retain max-zero initialization. See [tensor details]({% link docs/TMI/DTI-DKI.md %}#rician-correction-options).

### Processing time and parallel workers

Use `-n_cores` for independent voxel blocks, with one BLAS thread per worker.
Clipping is generally faster than expectation fitting; runtime depends on voxel
count, SH order, SNR and active constraints. Parallelization retains the same
objective, initialization, dense-grid check and rank-one processing.

### Standard Model parameter fitting

Then, the SM parameters are estimated from the rotational invariants of the signal (up to L=4 as default option). Unlike conventional parameter estimation approaches which rely on an analytical forward model, e.g. maximum likelihood, here we use data-driven machine learning (ML) regression. This is done by applying a sufficiently flexible nonlinear regression to _training data_ generated with the forward model of interest, considering a wide distribution of model parameters, the noise level, and the protocol that was used. Then, such regression is applied to the data of interest.

For typical SNR values found in clinical dMRI experiments, the optimal regression, i.e. minimizing MSE of the training data, can be achieved already by a cubic polynomial, as it  captures all relevant degrees of freedom in the data represented by the set of rotational invariants.

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
Both modes use linear least squares for SH fitting by default. Adding
`-RicianExp -sigma sigma.nii` accounts for Rician bias through the expected
magnitude observation model and enforces a **nonnegative underlying SH signal**.
Alternatively, `-RicianClip -sigma sigma.nii` fits clipped-corrected measurements
with the same nonnegativity constraint. If both flags are present,
`-RicianExp` takes precedence for SH/SMI fitting.

### Choosing an estimation mode

| Requested output | Example command |
|---|---|
| SMI maps with the default linear SH fit | `tmi dwi.nii smi_output -SMI -sigma sigma.nii` |
| SMI maps with Rician bias expectation correction and nonnegativity | `tmi dwi.nii smi_rician -SMI -RicianExp -sigma sigma.nii` |
| SMI maps with clipped correction and signal nonnegativity | `tmi dwi.nii smi_clip -SMI -RicianClip -sigma sigma.nii` |
| Signal rotational invariants only, using default LLS | `tmi dwi.nii sh_output -SH` |
| Signal rotational invariants only, using Rician expectation and signal nonnegativity | `tmi dwi.nii sh_rician -SH -RicianExp -sigma sigma.nii` |
| Signal rotational invariants only, using clipped correction and signal nonnegativity | `tmi dwi.nii sh_clip -SH -RicianClip -sigma sigma.nii` |

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

**`-SH -RicianExp` uses the expected Rician magnitude signal and a nonnegativity
constraint on the underlying SH signal.** Supply a noise map with `-sigma`:

```bash
tmi dwi.nii sh_rician -SH -RicianExp -sigma sigma.nii
```

The constraint applies at acquired directions and on a dense sphere grid;
the default LLS fit has no such constraint. Supplying `-sigma` alone does not
switch from LLS to Rician fitting. The same estimator choice applies to `-SMI`.
See [Rician bias expectation correction](#rician-bias-expectation-correction-with-a-nonnegative-signal)
for the observation model and initialization.

For clipped correction followed by a nonnegative linear SH fit, use:

```bash
tmi dwi.nii sh_clip -SH -RicianClip -sigma sigma.nii
```

All three commands bypass SMI priors and regression and save separate `S0.nii`,
`S2.nii`, and higher even-order files through the maximum fitted order. Each
contains one volume per shell, identified in `SH_shells.json`. Add `-lmaxSH 6`
to fit through order 6 and save `S0.nii`, `S2.nii`, `S4.nii`, and `S6.nii`.
In `-SH` mode, `-lmaxSH` sets the SH fitting order itself. See the
[standalone SH options]({% link docs/TMI/usage.md %}#-sh).

### Choose a fitting order for each shell

Orders 2, 4, and 6 require at least 6, 15, and 28 independent directional
measurements, respectively. With `-SH`, **each shell's maximum SH order can be
specified independently** using a comma-separated list. This works with
LLS, RicianExp and RicianClip, for example:

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
b=0/isotropic shells; these always retain their arithmetic means. A single
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
between LLS, constrained Rician-expectation fitting and constrained
clipped-corrected LLS applies to both output modes.

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

### Rician bias expectation correction with a nonnegative signal

With `tmi -SMI -RicianExp -sigma sigma.nii`, the SH fit instead minimizes

$$\sum_i \left(\mathbb{E}[M_i \mid (Yc)_i,\sigma] - M_i\right)^2,$$

where `Yc` is the underlying SH signal, `M` is the measured magnitude, and `sigma` is the standard deviation of each Gaussian component before magnitude reconstruction. The underlying signal is constrained to be nonnegative at all acquired directions and on a deterministic 65,536-direction sphere grid. These are linear constraints on the SH coefficients, not clipping of measurements or individual coefficients. A finite grid does not certify nonnegativity everywhere on the continuous sphere.

Every nonlinear fit starts from the ordinary linear SH fit to `sqrt(max(0, M²−sigma²))`. The clipping sets negative values of the squared difference to zero before taking the square root. This correction is used only to initialize the coefficients; the nonlinear objective still uses the original magnitude measurements. The SH projection of corrected measurements can itself have negative lobes, so the constrained solver accepts an initially infeasible coefficient vector rather than shifting that starting estimate.

The Rician mean is evaluated using scaled modified Bessel functions and an analytic Jacobian. A batched Levenberg–Marquardt fit, with SciPy refinement where necessary, is retained when it already satisfies the constraints. Otherwise, SLSQP starts from the same corrected-LLS coefficients. Its objective is multiplied by a positive scalar based on the SH design matrix, with the stopping tolerance scaled consistently; this preserves the fitted objective and reduces overshooting into spurious zero stationary points. Constraint generation adds the worst violated directions until the entire dense grid is feasible. An all-zero stationary solution is retried with a trust-region solver from those same corrected-LLS coefficients, because the Rician derivative vanishes at zero even when that point is not a minimum. If constrained optimization does not converge, a warning is issued and the best grid-feasible coefficients are returned. The constraint applies to the per-shell signal fit; the existing subsequent rank-one processing of higher-order invariants is unchanged and does not alter S0.

This extends the Rician-expectation approach from S. Coelho's SHestimation with signal nonnegativity, retaining TMI's SH normalization, phase, shell grouping, and rank-one processing of invariants. The constraint prevents signed angular lobes from fitting positive magnitudes through the even Rician mean function. As in that reference, b=0 and isotropic-shell means remain unchanged. This is a fit of the expected magnitude, not a maximum-likelihood estimator or an unbiasedness guarantee. The downstream SMI regression and its noise training model remain unchanged; their statistical bias may persist, especially at low SNR.

Use a finite, nonnegative noise map aligned to the DWI in the same intensity units. At zero noise, the estimator reduces to linear least squares with the same signal constraints; already feasible linear coefficients are retained exactly. The option is intended for magnitude data that retain a Rician noise floor, not data already corrected with `designer -rician` or complex/phase-unwound denoising. When denoising reduces residual variance but retains the original Rician mean bias, keep the original sigma for this mean model, not the reduced residual SD. Match any signal rescaling and spatial alignment; other preprocessing may change the mean model itself.

The same expectation model, clipped-LLS initialization, and nonnegativity
constraint are used by `-SH -RicianExp`. In that mode, processing stops after
the signal rotational invariants have been saved.

For difficult fits, `-RicianSH_maxiter N` optionally caps the cumulative
constrained-refinement iterations per voxel and shell, including restarts
and constraint-generation rounds. The fast initial batched fit retains its
standard limit. Without this option, the existing solver limits are unchanged.
At the cap, the best grid-feasible coefficients are returned with a warning;
an early-stopped fit is not a convergence guarantee. Very small limits can
change the estimates substantially, even with clipped-LLS initialization,
because the initial SH signal need not satisfy nonnegativity.

When this option is supplied, both `-SH` and `-SMI` also save
`RicianSH_early_stop.nii`, `RicianSH_iterations.nii`, and
`RicianSH_nonconverged.nii`. Each volume describes one shell before rank-one
processing. Zero iterations means no constrained refinement was needed;
b=0 and isotropic shells also have zero counts. `RicianSH_diagnostics.json`
identifies the shells, iteration cap, and per-shell counts. The early-stop
mask marks cap exhaustion; the nonconverged mask also includes other failures
of constrained refinement. These diagnostics do not describe convergence of
the initial batched fit, whose existing warnings remain separate.

### Clipped correction with a nonnegative signal

Use `-RicianClip -sigma sigma.nii` with either `-SH` or `-SMI` to fit

$$\min_c \sum_i ((Yc)_i-\widetilde M_i)^2,\qquad
\widetilde M_i=\sqrt{\max(0,M_i^2-\sigma^2)},\qquad Y_{\rm grid}c\geq0.$$

The corrected magnitudes are the fitting target, rather than only the
initialization. This is a convex constrained linear least-squares problem;
it does not evaluate the expected Rician magnitude. The correction subtracts
**one** `sigma²`, and negative squared differences are set to zero before
taking the square root. It is an approximate correction, not an unbiasedness
guarantee.

The same acquired directions and 65,536-direction sphere grid enforce
nonnegativity. If the ordinary SH fit to the clipped data is already feasible,
its coefficients are retained exactly. Otherwise, constrained refinement
starts from those coefficients. The existing b=0 and isotropic-shell raw
arithmetic means, SH normalization, cross-shell rank-one processing and SMI
regression remain unchanged. Nonnegativity is enforced before rank-one
processing, as for `-RicianExp`.

Use uncorrected magnitude data and a matching noise map; do not apply the
clipped correction to the input yourself. When both `-RicianClip` and
`-RicianExp` are supplied, SH/SMI uses Rician-expectation fitting on the original
magnitudes, with nonnegativity; it does not fit the expectation to
already-corrected data. The optional `-RicianSH_maxiter` cap and its diagnostic
outputs also apply to this mode.

### The same correction choices for DTI and DKI

The two flags also apply to `-DTI`, `-DKI` and `-WDKI`, including in combined
runs. Both require an aligned original sigma map in matching signal units.

| Flags | DTI/DKI/WDKI estimator | SH/SMI estimator |
|---|---|---|
| Neither | Existing WLLS on original magnitudes | Ordinary LLS |
| `-RicianClip` | Existing tensor pipeline on clipped-corrected measurements | Clipped-corrected LLS with signal nonnegativity |
| `-RicianExp` | Unconstrained Rician-mean NLLS on original observations | Rician expectation with signal nonnegativity |
| Both | Same as `-RicianExp` | Same as `-RicianExp` |

`-RicianClip` corrects selected tensor measurements once, including b0, and
allows the other tensor settings: constraints, AKC handling, smoothing,
polynomial regression, b0 restoration and GNC. Clipped zeros use the existing
positive log floor and can destabilize log-domain fits.

`-RicianExp` starts tensor NLLS from clipped ordinary LLS and ignores tensor
constraints and AKC replacement. It supports GNC and multi-TE subsets.
Expectation fitting with b0 restoration, fit smoothing, polynomial regression
or external tensors is deferred and produces an explicit error. These
restrictions do not apply to `-RicianClip`.

```bash
tmi dwi.nii tensor_clip -DTI -DKI -RicianClip -sigma sigma.nii
tmi dwi.nii combined -DTI -DKI -SH -RicianExp -sigma sigma.nii
```

The tensor expectation model handles b0 observations through the same Rician
mean as the other measurements. SH/SMI retain their raw arithmetic
b0/isotropic-shell mean convention. See [tensor correction details]({% link docs/TMI/DTI-DKI.md %}#rician-correction-options)
for equations, diagnostics and supported combinations.

### Processing time and parallel workers

Both `-RicianExp` and `-RicianClip` use `-n_cores` to fit independent voxel
blocks in parallel. For example, add `-n_cores 8` to use eight worker processes,
or `-n_cores 1` for serial fitting. The default follows TMI's available-cores
minus two convention. Each worker uses one BLAS thread to avoid multiplying
threads across workers. Parallelization leaves the fitting objective,
initialization, dense nonnegativity grid and cross-shell rank-one step unchanged.

Clipped correction solves a constrained linear problem and is generally
faster than Rician-expectation fitting. Runtime depends on voxel count, angular
order, SNR, active constraints, hardware and competing jobs; a short iteration
cap is not a substitute for convergence. A full dense-grid feasibility check
is retained even though each constrained subproblem starts with fewer grid
rows and adds violated constraints as needed.

### Standard Model parameter fitting

Then, the SM parameters are estimated from the rotational invariants of the signal (up to L=4 as default option). Unlike conventional parameter estimation approaches which rely on an analytical forward model, e.g. maximum likelihood, here we use data-driven machine learning (ML) regression. This is done by applying a sufficiently flexible nonlinear regression to _training data_ generated with the forward model of interest, considering a wide distribution of model parameters, the noise level, and the protocol that was used. Then, such regression is applied to the data of interest.

For typical SNR values found in clinical dMRI experiments, the optimal regression, i.e. minimizing MSE of the training data, can be achieved already by a cubic polynomial, as it  captures all relevant degrees of freedom in the data represented by the set of rotational invariants.

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

## Standard Model Imaging (SMI) MATLAB toolbox

In addition to the python implementation as a part of TMI, SMI is offered as a standalone matlab toolbox. [This MATLAB toolbox](https://github.com/NYU-DiffusionMRI/SMI/blob/master/SMI.m) contains all necessary functions for parameter estimation of the Standard Model (SM) of diffusion in white matter. Check [our recent paper](https://arxiv.org/pdf/2202.02399.pdf) for details on this implementation and on the Standard Model in general. Below we provide instructions on how to run the toolbox. See the [example.m](https://github.com/NYU-DiffusionMRI/SMI/blob/master/example.m) script that performs the parameter estimation in an [example dataset](https://cai2r.net/resources/standard-model-of-diffusion-in-white-matter-the-smi-toolbox/).

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

---

## Parameter estimation
First, the rotational invariants of the diffusion signal are estimated. This is done using a least squares estimator. 

Then, the SM parameters are estimated from the rotational invariants of the signal (up to L=4 as default option). Unlike conventional parameter estimation approaches which rely on an analytical forward model, e.g. maximum likelihood, here we use data-driven machine learning (ML) regression. This is done by applying a sufficiently flexible nonlinear regression to _training data_ generated with the forward model of interest, considering a wide distribution of model parameters, the noise level, and the protocol that was used. Then, such regression is applied to the data of interest.

For typical SNR values found in clinical dMRI experiments, the optimal regression, i.e. minimizing MSE of the training data, can be achieved already by a cubic polynomial, as it  captures all relevant degrees of freedom in the data represented by the set of rotational invariants.

---
layout: default
title: RICE
has_children: false
parent: TMI
nav_order: 4
---

# Rotational Invariants of the Cumulant Expansion
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## RICE toolbox (Rotational Invariants of the Cumulant Expansion)

{: .note }
Rice is not yet available as part of the TMI toolbox, however its MATLAB implementation is online.

This [MATLAB toolbox](https://github.com/NYU-DiffusionMRI/RICE) contains all necessary functions for parameter estimation of the O(b^2) cumulant expansion for arbitrary axially symmetric b-tensors. Check [our recent paper](https://arxiv.org/pdf/2202.02399.pdf) for details on this implementation and [this book chapter](https://academic.oup.com/book/24921/chapter/188767614) for information on the cumulant expansion in general. Below we provide instructions on how to run the toolbox. See the [example_RICE.m](https://github.com/NYU-DiffusionMRI/RICE/blob/main/example_RICE.m) script that performs the parameter estimation in an [example dataset](https://cai2r.net/resources/standard-model-of-diffusion-in-white-matter-the-smi-toolbox/).

The toolbox also allows the parameter estimation for minimal DKI and minimal RICE protocols proposed in [this paper](dead link).

---

## Overview: The cumulant expansion in diffusion MRI
<img width="1604" alt="Ob2_cumulant_expansion_RICE" src="https://user-images.githubusercontent.com/54751227/203182843-a1097e1d-9bd8-4a88-a60e-99f0d56a5104.png">


### LTE data
For conventional dMRI data, linear tensor encoding (LTE), one can represent low-b data with the O(b) cumulant expansion as shown in Eq. (1). This is simply DTI, and it can represent DWIs up to ~b=1200 ms/mm^2.
For higher b-values (up to ~b=2500 ms/mm^2), one can represent the DWIs with the O(b^2) cumulant expansion shown in Eq. (2). This is DKI.

### Multiple b-tensors
If we consider multiple b-tensor shapes (axially symmetric) as shown in the figure below, where β parametrizes the b-tensor shape. Most common examples are: β=1 for LTE (B has only one nonzero eigenvalue), β=0 for STE (B has 3 equal nonzero eigenvalues), and β=-1/2 for PTE (B has 2 equal nonzero eigenvalues).
<img width="1206" alt="axSymB" src="https://user-images.githubusercontent.com/54751227/197211877-1d589475-8835-4bcd-861a-35ee3f9a297f.png">
We see that for O(b) signals (Eq. (3) ). This representation is still DTI.
However, for O(b^2) a new tensor shows up: the diffusion covariance tensor, C, see Eq. (4). C is more general than kurtosis, actually it contains all the information of the kurtosis tensor plus some extra.

---

## Example use cases
The [example_RICE.m](https://github.com/NYU-DiffusionMRI/RICE/blob/main/example_RICE.m) script shows some examples on how to run the full RICE fitting and also the minimal DKI and minimal RICE ones. We provide example datasets for these, [check this link](https://cai2r.net/).

Briefly, the basic usage of the code is as follows:
```
% RICE toolbox parameter estimation example

type = 'fullRICE';  %  Estimate full D and C tensors from LTE + PTE data (WLLS)
CSphase = 1;        % Use Condon-Shortley phase in spherical harmonics definition
nsl_flag = 1;       % Use local nonlinear smoothing for fitting to boost SNR

[b0, tensor_elems, RICE_maps, DIFF_maps] = RICE.fit(DWI, b, dirs, bshape, mask, CSphase, type, nls_flag)

% Compute fiber basis projections (axial and radial diffusivities and kurtosis)
DKI_maps = RICE.get_DKI_fiberBasis_maps_from_4D_DW_tensors(tensor_elems, mask, CSphase);
```
See the help in RICE.fit and RICE.get_DKI_fiberBasis_maps_from_4D_DW_tensors for more details.

The following options are available in RICE.fit for the input argument 'type': (parameter count does not include s0)
- 'minimalDTI': only MD is fit (1 elem, [D00])
- 'fullDTI': full diffusion tensor is fit (6 elem, [D00 D2m])
- 'minimalDKI': full diffusion tensor and MK are fit (7 elem, [D00 D2m S00])
- 'minimalDKI_iso': only MD  and MK are fit (2 elem, [D00 S00])
- 'fullDKI': full diffusion and kurtosis tensors are fit (21 elem, [D00 D2m S00 S2m S4m])
- 'minimalRICE': full diffusion tensor + MK + A0 are fit (8 elem, [D00 D2m S00 S2m A00])
- 'fullRICE': full diffusion and covariance tensors are fit (27 elem, [D00 D2m S00 S2m S4m A00 A2m])

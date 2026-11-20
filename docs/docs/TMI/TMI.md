---
layout: default
title: TMI
has_children: true
nav_order: 3
permalink: /docs/TMI
---

# Tissue Microstructure Imaging

`tmi` (tissue microstructure imaging) is a python tool for diffusion MRI parameter estimation. It has been designed along with the `designer` preprocessing package and they are meant to be used in conjunction - `designer` for image pre-processing and `tmi` for parameter estimation. `tmi` includes:
- [singe shell diffusion tensor estimation]({% link docs/TMI/DTI-DKI.md %}#the-diffusion-and-kurtosis-tensors)
- [multi shell kurtosis tensor estimation]({% link docs/TMI/DTI-DKI.md %}#the-diffusion-and-kurtosis-tensors)
- [Rician bias and estimator alternatives]({% link docs/TMI/Rician-bias.md %}), using `-RicianClipMax`, `-RicianClipAbs` or `-RicianExp` for DTI/DKI and SH/SMI
- [W tensor parameter extraction]({% link docs/TMI/DTI-DKI.md %}#tensor-parameters)
- [DKI outlier replacement]({% link docs/TMI/DTI-DKI.md %}#outlier-replacement)
- [Standard Model Imaging (SMI)]({% link docs/TMI/SMI.md %}), including Rician bias expectation correction with optional SH signal nonnegativity
- [Signal rotational invariants without SMI]({% link docs/TMI/SMI.md %}#signal-rotational-invariants-without-smi), using `-SH` with LLS, corrected LLS or Rician expectation fitting



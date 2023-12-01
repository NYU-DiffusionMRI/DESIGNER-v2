---
layout: default
title: TMI
has_children: true
nav_order: 1
permalink: /docs/TMI
---

# Tissue Microstructure Imaging

`tmi` (tissue microstructure imaging) is a python tool for diffusion MRI parameter estimation. It has been designed along with the `designer` preprocessing package and they are meant to be used in conjunction - `designer` for image pre-processing and `tmi` for parameter estimation. `tmi` includes:
- [singe shell diffusion tensor estimation]({% link docs/TMI/DTI-DKI.md %}#The-Diffusion-and-Kurtosis-Tensors)
- [multi shell kurtosis tensor estimation]({% link docs/TMI/DTI-DKI.md %}#The-Diffusion-and-Kurtosis-Tensors)
- [W tensor parameter extraction]({% link docs/TMI/DTI-DKI.md %}#Tensor-parameters)
- [DKI outlier replacement]({% link docs/TMI/DTI-DKI.md %}#Outlier-replacement)
- [SMI and rotational invariant estimation]({% link docs/TMI/SMI.md %})



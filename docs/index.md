---
title: Home
layout: default
nav_order: 1
description: "Temp description"
permalink: /
---

# Diffusion MRI pre-processing (DESIGNER) and parameter estimation (TMI)

Welcome to the [NYU diffusion biophysics group](https://diffusion-mri.com) documentation page. The group is led by Els Fieremans and Dmitry Novikov, and has developed a number of tools for preprocessing and extracting parameters from diffusion MRI data. Features include denoising using MPPCA [(Veraart 2016)](https://www.sciencedirect.com/science/article/pii/S1053811916303949?via%3Dihub), gibbs artifact correction using RPG [(Lee 2021)](https://onlinelibrary.wiley.com/doi/10.1002/mrm.28830), the DESIGNER preprocessing pipeline [(Ades-Aron 2018)](https://www.sciencedirect.com/science/article/pii/S1053811918306827?via%3Dihub), the standard model imaging toolbox [(Coelho 2022)](https://www.sciencedirect.com/science/article/pii/S1053811922004104), and more!

This website documents DESIGNER version 2.0 and onwards. Prior versions of the DESIGNER package are officially deprecated and are no longer being maintained. GitHub code for designer-v2 is available [here](https://github.com/NYU-DiffusionMRI/DESIGNER-v2).

Both Designer-v2 and TMI are built as [mrtrix3](https://www.mrtrix.org) external packages. Please see our [installation]({{ site.baseurl }}{% link docs/designer/installation.md %}), [usage]({{ site.baseurl }}{% link docs/designer/usage.md %}), and [examples]({{ site.baseurl }}{% link docs/designer/examples.md %}) pages for information on how to run designer.

{: .note }
This webpage is still under development! Please excuse us for any missing information and don't hesitate to reach out if you have immediate issues with package installation.


# References
All toolboxes are subject to specific referencing. Each project is accompanied by at least one article for citation.

{: .ref }

> Ades-Aron, B., Veraart, J., Kochunov, P., McGuire, S., Sherman, P., Kellner, E., Novikov, D.S. & Fieremans, E. (2018). Evaluation of the accuracy and precision of the diffusion parameter EStImation with Gibbs and NoisE removal pipeline. Neuroimage, 183, 532-543.
>
> Veraart, J., Novikov, D. S., Christiaens, D., Ades-Aron, B., Sijbers, J., & Fieremans, E. (2016). Denoising of diffusion MRI using random matrix theory. Neuroimage, 142, 394-406.
>
> Lee, H. H., Novikov, D. S., & Fieremans, E. (2021). Removal of partial Fourier‐induced Gibbs (RPG) ringing artifacts in MRI. Magnetic resonance in medicine, 86(5), 2733-2750.
>
> Coelho, S., Baete, S. H., Lemberskiy, G., Ades-Aron, B., Barrol, G., Veraart, J., Novikov, D.S. & Fieremans, E. (2022). Reproducibility of the Standard Model of diffusion in white matter on clinical MRI systems. Neuroimage. 2022 Aug 15, 257:119290
>
> Chen, J.,  Ades-Aron, B., Lee, H.H., Mehrin, S., Pang, M., Novikov, D. S., Veraart, J., & Fieremans, E. (2024). Optimization and validation of the DESIGNER preprocessing pipeline for clinical diffusion MRI in white matter aging. Imaging Neuroscience, 2 1–17


# Getting Help
The easiest source for references on these tools is right here. In addition, concerns raised in the [issues section](https://github.com/NYU-DiffusionMRI/DESIGNER-v2/issues) will reach the developers.


# Contributors
<ul class="list-style-none">
  <li class="d-inline-block mr-1">
     <a href="https://github.com/badesar1"><img src="https://avatars.githubusercontent.com/u/11949335?v=4" width="90" height="90" alt="https://github.com/badesar1"></a>
  </li>
  <li class="d-inline-block mr-1">
     <a href="https://med.nyu.edu/faculty/els-d-fieremans"><img src="https://avatars.githubusercontent.com/u/1108725?v=4" width="90" height="90" alt="https://med.nyu.edu/faculty/els-d-fieremans"></a>
  </li>
  <li class="d-inline-block mr-1">
     <a href="https://github.com/jchen33344"><img src="https://www.diffusion-mri.com/wp-content/uploads/2023/01/Screen-Shot-2023-01-07-at-10.59.32-PM.png" style="height:90px" alt="https://github.com/jchen33344"></a>
  </li>
  <li class="d-inline-block mr-1">
     <a href="https://github.com/leehhtw"><img src="https://www.diffusion-mri.com/wp-content/uploads/2021/05/Hong-Hsi-Lee.png)" style="height:90px" alt="https://github.com/leehhtw"></a>    
  </li>
  <li class="d-inline-block mr-1">
     <a href="https://github.com/rcleija"><img src="https://avatars.githubusercontent.com/u/44007271?v=4" width="90" height="90" alt="https://github.com/rcleija"></a>
  </li>
  <li class="d-inline-block mr-1">
     <a href="https://santiagocoelho.github.io"><img src="https://avatars.githubusercontent.com/u/54751227?v=4" width="90" height="90" alt="https://santiagocoelho.github.io"></a>
  </li>
  <li class="d-inline-block mr-1">
     <a href="https://github.com/d-novikov"><img src="https://avatars.githubusercontent.com/u/29991818?v=4" width="90" height="90" alt="https://github.com/d-novikov"></a>
  </li><!--
  <li class="d-inline-block mr-1">
     <a href="https://github.com/jelleveraart"><img src="https://avatars.githubusercontent.com/u/6860257?v=4" width="90" height="90" alt="https://github.com/jelleveraart"></a>
      </li>
<li class="d-inline-block mr-1">
     <a href="https://www.linkedin.com/in/gregory-lemberskiy/"><img src="https://avatars.githubusercontent.com/u/1512844?v=4" width="90" height="90" alt="https://www.linkedin.com/in/gregory-lemberskiy/"></a>
  </li>
  <li class="d-inline-block mr-1">
     <a href="https://www.linkedin.com/in/valentinmri/"><img src="https://media.licdn.com/dms/image/C4E03AQFQ9Bs9qvt2Hg/profile-displayphoto-shrink_400_400/0/1644721673667?e=1691625600&v=beta&t=qMzM60JbWuY3VqxWjF6H0K2fbw-cQRncVTfLK654qog" width="90" height="90" alt="https://www.linkedin.com/in/valentinmri/"></a>
  </li>
  <li class="d-inline-block mr-1">
     <a href="https://www.linkedin.com/in/ying-liao-nyu/"><img src="https://media.licdn.com/dms/image/C4E03AQEfT0sPV34ImQ/profile-displayphoto-shrink_400_400/0/1627317717468?e=1691625600&v=beta&t=eqAZqafA7ZusHRkVzBwwkA4r6yAcZMC0lbRhlCmU8Ig" width="90" height="90" alt="https://www.linkedin.com/in/ying-liao-nyu/"></a>
  </li>
  <li class="d-inline-block mr-1">
     <a href="https://github.com/aAbdz"><img src="https://avatars.githubusercontent.com/u/29164686?v=4" width="90" height="90" alt="https://github.com/aAbdz"></a>
  </li>-->
</ul>


# NYU License
Copyright (c) 2023 New York University
       
Permission is hereby granted, free of charge, to any non-commercial entity
('Recipient') obtaining a copy of this software and associated
documentation files (the 'Software'), to the Software solely for
non-commercial research, including the rights to use, copy and modify the
Software, subject to the following conditions: 

1. The above copyright notice and this permission notice shall be
included by Recipient in all copies or substantial portions of the
Software. 

2. THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BELIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF ORIN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE. 

3. In no event shall NYU be liable for direct, indirect, special,
incidental or consequential damages in connection with the Software.
Recipient will defend, indemnify and hold NYU harmless from any claims or
liability resulting from the use of the Software by recipient. 

4. Neither anything contained herein nor the delivery of the Software to
recipient shall be deemed to grant the Recipient any right or licenses
under any patents or patent application owned by NYU. 

5. The Software may only be used for non-commercial research and may not
be used for clinical care. 

6. Any publication by Recipient of research involving the Software shall
cite the references listed below.

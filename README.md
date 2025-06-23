# designer-v2

## Introduction
Designer is a python tool for diffusion MRI preprocessing. It includes:

* denoising using MPPCA (or optionally using patch2self through dipy)
* RPG Gibbs artifact correction
* Rician bias correction
* EPI distortion correction
* Eddy current and motion correction
* b0 normalization for multiple input series
* b1 inhomogeneity correction

For complete documentation on designer & tmi installation and usage please visit the [documentation](https://nyu-diffusionmri.github.io/DESIGNER-v2).

## Installation

designer-v2 is made available as a [python package](https://pypi.org/project/designer2/), and a [docker image](https://hub.docker.com/r/nyudiffusionmri/designer2/) bundling all dependencies. See the [Installation](https://nyu-diffusionmri.github.io/DESIGNER-v2/docs/designer/installation/) documentation for detailed steps


## Example usage for "meso" data

An example script can be found in the examples folder. It is copied here as well. Preprocessing and fitting are now split into two separate functions: designer for preprocessing and tmi for fitting.

```
#!/bin/bash

datapath=/mnt/labspace/Projects/MESO_v2.0/ALLSUBJS_2.0/M9734
meso1=M9734_074YM_DIFF_meso.nii
meso2=M9734_074YM__DIFF_meso_research.nii
pa=M9734_074YM_DIFF_meso_PA.nii

cd $datapath

designer \
-denoise -algorithm jespersen -extent 7,7,7 \
-degibbs \
-mask \
-scratch designer2_processing_test -nocleanup \
$meso1,$meso2 designer2_test.mif

tmi \
-DTI -DKI -WDKI -SMI \
-mask designer2_processing_test/brain_mask.nii \
-sigma designer2_processing_test/sigma.nii \
-nocleanup \
designer2_test.mif designer2_params_test
```

## Development Setup

### Using Dev Container

DESIGNER-v2 provides a development container configuration that ensures all developers work with the same development environment, making it easier to collaborate and avoid "it works on my machine" issues. This setup is particularly useful for ones who want to contribute to DESIGNER.

#### Prerequisites

1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop/)
2. Install [Visual Studio Code](https://code.visualstudio.com/) (or any other IDEs that support dev containers)
3. Install the [Dev Containers extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) in VS Code

#### Getting Started

1. Open the project in VS Code. When it detects the Dev Container configuration, it will show a notification. Click "Reopen in Container" or:
   - Press `Cmd/Ctrl + Shift + P` (or click "View" -> "Command Palette" in a top menu bar)
   - Type "Dev Containers: Reopen in Container" and select it

2. Wait for the container to build (this may take 10-15 minutes only for the first time). The container includes:
   - Python 3.12 environment
   - FSL, MRtrix3, and ANTs pre-installed
   - All necessary DESIGNER dependencies

3. Once inside the container, DESIGNER will be automatically installed in editable mode (`pip install -e .` defined in devcontainer.json). You can:
   - Make code changes and test them immediately
   - Use all DESIGNER commands from the VS Code integrated terminal (Press ``Ctrl + ` `` or click "View" -> "Terminal"). All changes will be reflected immediately.
     ```bash
     # Try DESIGNER commands
     designer -help
     ```

#### Troubleshooting

- If the container fails to build, try:
  ```bash
  docker system prune -a
  ```
  Then rebuild the container using VS Code.

- If you need to rebuild the container (e.g., when dependencies are updated or devcontainer.json is modified):
  - Press `Cmd/Ctrl + Shift + P` (or click "View" -> "Command Palette")
  - Type "Dev Containers: Rebuild Container" and select it

- If using SSH-configured repositories, `git push` from the dev container's terminal will fail (container cannot access local SSH config). Workarounds:
  - Use VS Code's Source Control panel, or
  - Use your local machine's terminal, or
  - Follow instructions in [VS Code Dev Container Git Credentials](https://code.visualstudio.com/remote/advancedcontainers/sharing-git-credentials)

#### Notes

- **VS Code Extensions**:
  - Default extensions are automatically installed based on devcontainer.json
  - You can install additional extensions for your personal use (Manually installed extensions will be reset when rebuilding the container)

#### Further References

For more details about the development environment:
- Container build configuration: [.devcontainer/Dockerfile](.devcontainer/Dockerfile)
- Dev container settings: [.devcontainer/devcontainer.json](.devcontainer/devcontainer.json)
- Official VS Code Dev Containers documentation: [Developing inside a Container](https://code.visualstudio.com/docs/devcontainers/containers)

# References
1. Ades-Aron, B., Veraart, J., Kochunov, P., McGuire, S., Sherman, P., Kellner, E., … & Fieremans, E. (2018). Evaluation of the accuracy and precision of the diffusion parameter EStImation with Gibbs and NoisE removal pipeline. *Neuroimage*, 183, 532-543.

2. Chen, J., Ades-Aron, B., Lee, H.H., Mehrin, S., Pang, M., Novikov, D. S., Veraart, J., & Fieremans, E. (2024). Optimization and validation of the DESIGNER preprocessing pipeline for clinical diffusion MRI in white matter aging. *Imaging Neuroscience*, 2, 1–17.

3. Tournier, J.-D., Smith, R. E., Raffelt, D., Tabbara, R., Dhollander, T., Pietsch, M., Christiaens, D., Jeurissen, B., Yeh, C.-H., & Connelly, A. (2019). MRtrix3: A fast, flexible and open software framework for medical image processing and visualisation. *NeuroImage*, 202, 116137.


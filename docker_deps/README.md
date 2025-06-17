This folder contains Docker files for DESIGNER's main dependencies. All images are hosted on the [nyudiffusionmri Docker Hub](https://hub.docker.com/u/nyudiffusionmri).

The following Docker files are included:
- `Dockerfile_mrtrix`: For MRtrix3 dependency
- `Dockerfile_fsl`: For FSL dependency
- `Dockerfile_ants`: For ANTs dependency
<br>

All Docker images are tagged with dates in the `YYYY-MM-DD` format.

Note: Prior to documenting build and version information, we used images from [Tom's Docker Hub](https://hub.docker.com/u/twom). These have been migrated to the nyudiffusionmri Docker Hub with new date-based tags. The oldest dated image for each dependency matches exactly what we used in [DESIGNER v2.0.13](https://github.com/NYU-DiffusionMRI/DESIGNER-v2/releases/tag/v2.0.13).

You can find the Docker files used to build specific images by combining a commit ID with the dependency file path in the GitHub URL. For example, to find the MRtrix Dockerfile for commit `8ea6505d1b079dfa0baa9845f5a4681411a6aa37`:
- [https://github.com/NYU-DiffusionMRI/DESIGNER-v2/blob/8ea6505d1b079dfa0baa9845f5a4681411a6aa37/docker_deps/Dockerfile_mrtrix](https://github.com/NYU-DiffusionMRI/DESIGNER-v2/blob/8ea6505d1b079dfa0baa9845f5a4681411a6aa37/docker_deps/Dockerfile_mrtrix)

Dependency versions are listed in both the Docker files and the tables below.


## MRtrix3

[Docker Hub Link](https://hub.docker.com/r/nyudiffusionmri/mrtrix3/tags)

| Docker Tag | Commit Hash for Docker File | Dependency Version |
| :--------: | :-------------------------: | :----------------: |
| 2025-06-16 |           9abe8fa           | 205dd53ef (3.4.0)  |
| 2024-02-09 |             N/A             | 205dd53ef (3.4.0)  |


Notes:
- Due to MRtrix3's API changes (and subsequent deprecation of their dev branch), we've pinned the commit for the `2025-06-16` image. The pinned version can be found in the [MRtrix3 Github](https://github.com/MRtrix3/mrtrix3/tree/205dd53ef).
- The `2024-02-09` image is identical to [twom/mrtrix3:dev-latest](https://hub.docker.com/layers/twom/mrtrix3/dev-latest/images/sha256-7630a4cd709cd7b9967f6db5dae112cd3f7be694fb5fc69c6e8ce1c0c3689d0c).


## FSL

[Docker Hub Link](https://hub.docker.com/r/nyudiffusionmri/fsl/tags)

| Docker Tag | Commit Hash for Docker File | Dependency Version |
| :--------: | :-------------------------: | :----------------: |
| 2025-06-16 |           9abe8fa           |       6.0.7        |
| 2024-02-10 |             N/A             |        N/A         |

Notes:
- The `2024-02-10` image is identical to [twom/fsl:6.0](https://hub.docker.com/layers/twom/fsl/6.0/images/sha256-4edc064ee849b8d05aaf98049f0d64c4d07dc27e9d61ad6211c1c7559625d58d).


## ANTs

[Docker Hub Link](https://hub.docker.com/r/nyudiffusionmri/ants/tags)

| Docker Tag | Commit Hash for Docker File | Dependency Version |
| :--------: | :-------------------------: | :----------------: |
| 2025-06-16 |           9abe8fa           |       2.5.4        |
| 2024-02-10 |             N/A             |       2.5.4        |

Notes:
- The `2024-02-10` image is identical to [twom/ants:v2.5.4](https://hub.docker.com/layers/twom/ants/v2.5.4/images/sha256-eb186b9a6959c60e360a4c6d38f36adbfac6709e1cc464a63b4fef4635d5fbfc).


## Updating Dependencies

After modifying a dependency's Docker file, follow these steps to update the image:

1. Build the image using:
```sh
docker build --platform=linux/amd64 -f docker_deps/$DEPENDENCY_DOCKER_FILE_NAME -t nyudiffusionmri/$DEPENDENCY_DOCKER_HUB_REPO:$NEW_DATE_TAG .
```

For example:
```sh
docker build --platform=linux/amd64 -f docker_deps/Dockerfile_mrtrix -t nyudiffusionmri/mrtrix3:2025-06-16 .
```

2. Update the `COPY` instruction in the main DESIGNER Dockerfile (located in the project root) with the new `$NEW_DATE_TAG`. Then build and test the DESIGNER image locally to ensure compatibility with the updated dependency. You can run tests similar to those in `.circleci/config.yml`:

```sh
# Build the test image with the updated dependency
docker build --platform=linux/amd64 -t designer2:test .

# Run the container and verify functionality
docker run --rm -it --platform=linux/amd64 designer2:test /bin/bash

# Inside the container, run these commands to verify the installation:
designer -version
flirt -version
mrinfo -version
N4BiasFieldCorrection --help  # Verify ANTs command
python -c "from lib import rpg; rpg.unring()"  # Verify function availability
```

3. After successful local testing, push the new dependency image to Docker Hub:
```sh
docker push $IMAGE_NAME_WITH_TAG
```

For example:
```sh
docker push nyudiffusionmri/mrtrix3:2025-06-16
```

4. Update this README by adding the new image information (tag, commit hash for Docker file, and dependency version) to the appropriate table.

5. Create a Pull Request and verify that the new build passes the CI pipeline. If the build fails, you'll need to modify the dependency Docker file and repeat the process from step 1. Once the PR is approved and merged to the main branch, the **nyudiffusionmri/designer2:main** image will be automatically built and pushed to Docker Hub.
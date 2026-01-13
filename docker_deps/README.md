This folder contains Docker files for DESIGNER's main dependencies. All images are hosted on the [nyudiffusionmri Docker Hub](https://hub.docker.com/u/nyudiffusionmri).

The following Docker files are included:
- `Dockerfile_mrtrix`: MRtrix3 dependency
- `Dockerfile_fsl`: FSL dependency
- `Dockerfile_ants`: ANTs dependency
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
| 2026-01-13 |           faa88b0           | 3.0.8              |
| 2025-06-16 |           9abe8fa           | 205dd53ef (3.0.4)  |
| 2024-02-09 |             N/A             | 205dd53ef (3.0.4)  |


Notes:
- Due to MRtrix3's API changes and subsequent deprecation of their dev branch, we've pinned the commit for the `2025-06-16` image. The pinned version can be found in the [MRtrix3 Github](https://github.com/MRtrix3/mrtrix3/tree/205dd53ef).
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
| 2025-06-20 |           26a0284           |       2.5.4        |
| 2025-06-16 |           9abe8fa           |       2.5.4        |
| 2025-05-14 |             N/A             |       2.5.4        |

Notes:
- The `2025-06-16` image has a bug with `N4BiasFieldCorrection` binary. **Do not use this image.**
- The `2025-05-14` image is identical to [twom/ants:v2.5.4](https://hub.docker.com/layers/twom/ants/v2.5.4/images/sha256-eb186b9a6959c60e360a4c6d38f36adbfac6709e1cc464a63b4fef4635d5fbfc).


## Updating Dependencies

When modifying a dependency's Docker file, follow these steps to update the image and run e2e tests with the updated dependency in dev container:

Due to the issue that with dev container cannot  `COPY --from` local image built with `--platform=linux/amd64` flag, we first need to push the image to individual account's (not **nyudiffusionmri**) Docker Hub repository. New Docker Hub account registration is required if you don't have one.

1. Build the image using:
```sh
docker build --platform=linux/amd64 -f docker_deps/$DEPENDENCY_DOCKER_FILE_NAME -t $DOCKER_HUB_USERNAME/$DEPENDENCY_DOCKER_HUB_REPO:$NEW_TAG .
```

For example:
```sh
docker build --platform=linux/amd64 -f docker_deps/Dockerfile_mrtrix -t ehddbs6425/mrtrix3:testing .
```
<br>


2. Push the image to Docker Hub, then update the `COPY` instruction in the dev container Dockerfile (`/.devcontainer/Dockerfile`) with the image just pushed.

For example:
```sh
docker push ehddbs6425/mrtrix3:testing

# Update the .devcontainer/Dockerfile
...
COPY --from=ehddbs6425/mrtrix3:testing /usr/local/mrtrix3/build /usr/local/mrtrix3_build
...
```

3. Rebuild the dev container and check if the new dependency is installed correctly. Then, run e2e tests to verify updated dependency doesn't break the existing functionality.

For example:
```sh
# Inside dev container (cwd: PROJECT ROOT), check if mrtrix version is updated:
mrinfo -version

# Run e2e tests (Refer to the main README.md for more details on how to run e2e tests locally)
```

4. Once the e2e test passes and that we can confirm the new dependency is working correctly in dev container, rename the image with the `$NEW_DATE_TAG` (e.g., `2026-01-13`) so that we can push it to **nyudiffusionmri** Docker Hub later:

For example:
```sh
# Tag the image with the new date tag
docker tag ehddbs6425/mrtrix3:testing nyudiffusionmri/mrtrix3:2026-01-13
```

5. Update the `COPY` instruction of both dev container Dockerfile (`/.devcontainer/Dockerfile`) and main DESIGNER Dockerfile (located in the project root) with the new `$NEW_DATE_TAG`. Then, build the production (not dev container) DESIGNER image locally and perform bare minimum verification as shown below (these verifications are the same as in `.circleci/config.yml`):

For example:
```sh
# .devcontainer/Dockerfile
...
COPY --from=nyudiffusionmri/mrtrix3:2026-01-13 /usr/local/mrtrix3/build /usr/local/mrtrix3_build
...

# Dockerfile (project root)
...
COPY --from=nyudiffusionmri/mrtrix3:2026-01-13 /usr/local/mrtrix3/build /usr/local/mrtrix3_build
...

# cwd: PROJECT ROOT
# Build the test production image with the updated dependency
docker build --platform=linux/amd64 --target production -t designer2:local .

# Run the container
docker run --rm -it --platform=linux/amd64 designer2:local /bin/bash

# Inside the container, run these commands to verify the installation:
designer -version
flirt -version  # verify FSL installation
mrinfo -version
N4BiasFieldCorrection --help  # Verify ANTs command
python -c "from lib import rpg; rpg.unring()"  # Verify function availability
```
<br>


6. Finally, push the new dependency image to **nyudiffusionmri** Docker Hub:

For example:
```sh
docker push nyudiffusionmri/mrtrix3:2026-01-13
```
<br>

4. Update this README by adding the new image information (tag, commit hash for Docker file, and dependency version) to the appropriate table.
<br>

5. Push to GitHub and verify that the new build passes the CI pipeline (`build_on_commit`). If the build fails, modify the dependency Docker file and repeat the process from step 1. Otherwise, create a Pull Request. Once the PR is approved and merged to the main branch, the **nyudiffusionmri/designer2:main** image will be automatically built and pushed to Docker Hub.

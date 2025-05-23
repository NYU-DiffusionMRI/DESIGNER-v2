---
layout: default
title: Installation
parent: designer
has_children: false
nav_order: 3
---

# Installation
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---
## Running DESIGNER v2 using Docker (most convenient)

We provide convenient docker images to quickly get started with Designer without the need to install and its dependencies. <b>Official versions will be tagged</b> (nyudiffusionmri/designer2:\<tag\>; latest tag is v2.0.13), while we also maintain images following our main development branch (nyudiffusionmri/designer2:main). See [Docker Hub](https://hub.docker.com/repository/docker/nyudiffusionmri/designer2/tags?page=1&ordering=last_updated) for an overview of all releases.

To get started, install [Docker Desktop](https://www.docker.com/products/docker-desktop/) on your machine. If you are using an Apple Silicon device, navigate to the Docker app's general settings and enable the "Use Rosetta for x86/amd64 emulation on Apple Silicon" option for compatibility. We recommend setting the memory limit to at least 15GB for optimal performance.

Verifying can be done by running `docker run nyudiffusionmri/designer2:<tag> designer`

You can then start an interactive terminal in your container by running `docker run -it nyudiffusionmri/designer2:<tag> /bin/bash`

Files can be passed by mounting a local folder into the container: `docker run -it -v <path to local folder>:/data nyudiffusionmri/designer2:<tag> /bin/bash`. In this case you can run commands in the container directly.

 You can also run designer directly on the command line. Any files can be referenced by their absolute path in the container (e.g. `/data/<your file>`). 
`docker run -it -v <path to local data folder>:/data nyudiffusionmri/designer2:<tag> designer ...`

Note that you can replace `nyudiffusionmri/designer2:<tag>` by any other tag. 

For users who are working on HCP servers and cannot use Docker, DESIGNER may also be used with [Singularity](https://docs.sylabs.io/guides/2.6/user-guide/installation.html). The following can be run to singularity-ize our Docker image: `singularity pull docker://nyudiffusionmri/designer2:<tag>`

We can then run DESIGNER by running: `singularity run designer2_<tag>.sif designer...`


## Installing DESIGNER v2 (if you know what you are doing ;))

If running DESIGNER v2 using Docker is not your preferred method, you can try to install DESIGNER v2 as a [Python package](https://pypi.org/project/designer2/). This method requires installation of fsl and mrtrix3 on your local machine.

### PIP install

The new version of designer consists of a number of updates and improvements to the original implementation. Designer version 2 is written in python and c++ and is an external [Mrtrix3](https://www.mrtrix.org) project.

Designer can be installed in one line using pip: 
```
pip install designer2
```

To upgrade designer to the newest version after it has already been installed, just use the command:
```
pip install --upgrade designer2
```


### Requirements and Dependencies
Currently, installing designer locally requires a source installation of mrtrix3. A future update of the mrtrix3 precompiled binaries may eliminate these requirements, however as of today, users must follow the steps [here](https://mrtrix.readthedocs.io/en/latest/installation/build_from_source.html) to configure and build the source mrtrix3 package. 

Installing the `dev` branch is needed, as there are additional features in `dwifslpreproc` that have no yet made it to the main codebase. The source mrtrix3 installation can be completed by running the following lines after ensuring that mrtrix3 dependencies are installed properly.
```
git clone https://github.com/MRtrix3/mrtrix3.git
cd mrtrix3/
cmake -B build -DCMAKE_INSTALL_PREFIX=/path/to/installation/
cmake --build build
cmake --install build
```

After mrtrix3 has been successfully installed, users should create the `PYTHONPATH` environment variable and link it against the mrtrix3 python libraries:
```
export PYTHONPATH=</path/to/mrtrix3/lib>
```

For example, if a user clones the mrtrix3 github repository to the directory `/usr/local/mrtrix3`, they would run the command `export PYTHONPATH=/usr/local/mrtrix3/lib`. We recommend that users add this line to their `.bashrc` or `.bash_profile` so that designer is permanently configured properly in the users shell environment.

---

Designer also relies on FSL tools for EPI distortion correction and eddy current and motion correction. Users can install FSL [here](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation) and follow the steps outlined in the FSL documentation.

---

In addition to Mrtrix3 and FSL, designer should automatically install a number of additional python dependencies including:
- numpy (>=1.21.0,<2.0.0)
- scipy (>=1.9.0)
- numpy_groupies (>=0.9.0)
- antspyx (>=0.3)
- dipy (>=1.4)
- tqdm (>=4.62.3)
- joblib (>=1.2)
- cvxpy (>=1.2)
- pandas (>=1.5)

We recommend that users keep python and pip as up-to-date as possible to ensure that Designer runs smoothly. Designer was originally written in Python 3.9 and is known to be compatible with versions as low as Python 3.7 and up to 3.12. We also recommend that users run Designer within a python environment such as [conda](https://www.anaconda.com), [pyenv](https://github.com/pyenv/pyenv), or whichever you prefer. The choice of environment is up to the user, however environments are useful to manage potential dependency conflict.

{: .warning }
There is a known [issue](https://github.com/conda/conda/issues/12051) in conda environments on windows and Ubuntu that can prevent the installation of `designer`. The fix is to install `numpy-base` before running the `pip install designer2` command. Additionally, there is a known conflict with numpy 2.0.0 and cvxpy. If you have numpy 2.0.0 installed, try uninstalling it and installing numpy 1.21.0 before running the `pip install designer2` command.




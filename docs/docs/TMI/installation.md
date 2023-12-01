---
layout: default
title: Installation
parent: TMI
has_children: false
nav_order: 1
---

# Installation
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Installing TMI

TMI is installed alongside Designer. TMI is written in python and and is an external [Mrtrix3](https://www.mrtrix.org) project.

Please follow the [installation instructions]({{ site.baseurl }}{% link docs/designer/installation.md %}) for Designer.

---

## Requirements and Dependencies
Currently, TMI requires a source installation of mrtrix3 in order to run. A future update of the mrtrix3 precompiled binaries may eliminate these requirements, however as of today, users must follow the steps [here](https://mrtrix.readthedocs.io/en/latest/installation/build_from_source.html) to configure and build the source mrtrix3 package. After mrtrix3 has been successfully installed, users should create the `PYTHONPATH` environment variable and link it against the mrtrix3 python libraries:
```
export PYTHONPATH=</path/to/mrtrix3/lib>
```

For example, if a user clones the mrtrix3 github repository to the directory `/usr/local/mrtrix3`, they would run the command `export PYTHONPATH=/usr/local/mrtrix3/lib`. We recommend that users add this line to their `.bashrc` or `.bash_profile` so that designer is permanently configured properly in the users shell environment.

---

In addition to Mrtrix3, TMI should automatically install a number of additional python dependencies including:
- numpy
- scipy
- numpy_groupies
- antspyx (for MRI image i/o)
- dipy (diffusion denoising using patch2self)
- tqdm
- joblib
- cvxpy
- pandas

We recommend that users keep python and pip as up-to-date as possible to ensure that Designer runs smoothly. Designer was originally written in Python 3.9 and is known to be compatible with versions as low as Python 3.7 and up to 3.11. We also strongly recommend that users run TMI within a python environment such as [conda](https://www.anaconda.com), [pyenv](https://github.com/pyenv/pyenv), or whichever you prefer.
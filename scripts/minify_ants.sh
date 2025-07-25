# this script may not work properly.

# neurodocker minify is not supported on amd64.
# tried with Rosetta2, but 'neurodocker minify' command failed with internal retrozip error.

#!/bin/bash

set -euo pipefail

mkdir -p tmp

# 1. Convert to .mif and embed FSL gradients
mrconvert tests/data/D3/meso_slice.nii.gz tmp/in.mif \
  -fslgrad tests/data/D3/meso_slice.bvec tests/data/D3/meso_slice.bval

# 2. Extract b=0 volumes & average
dwiextract tmp/in.mif - -bzero \
  | mrmath - mean tmp/mean_bzero.mif -axis 3

# 3. Generate brain mask
dwi2mask mean tmp/in.mif tmp/mask.mif

# 4. Convert to NIfTI with canonical strides for ANTs
mrconvert tmp/mean_bzero.mif tmp/mean_bzero.nii \
  -strides +1,+2,+3
mrconvert tmp/mask.mif tmp/mask.nii \
  -strides +1,+2,+3


neurodocker generate docker \
    --pkg-manager apt \
    --base-image debian:bookworm \
    --ants version=2.4.3 \
    | docker build --platform linux/amd64 -t ants:2.4.3 -

docker run --rm -itd \
    --cap-add=SYS_PTRACE \
    --privileged \
    -v ./tmp:/data \
    --platform linux/amd64 \
    --name ants-container \
    ants:2.4.3

cmd="N4BiasFieldCorrection \
  -d 3 \
  -i /data/mean_bzero.nii \
  -w /data/mask.nii \
  -o [ /data/corrected.nii, /data/init_bias.nii ] \
  -s 4 \
  -b [ 100, 3 ] \
  -c [ 1000, 0.0 ]"

neurodocker minify \
    --container ants-container \
    --dir /opt \
    "$cmd"

docker export ants-container | docker import - ants:2.4.3-min-n4biascorr

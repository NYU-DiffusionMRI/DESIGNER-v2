#!/bin/bash

set -euo pipefail

FS_TGZ="fs.tgz"
FS_URL="https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.3.2/freesurfer-linux-ubuntu18_amd64-7.3.2.tar.gz"
FS_MODEL_PATH="freesurfer/models/synthstrip.1.pt"

MODEL_DIR="tests/models"
MODEL_OUT="${MODEL_DIR}/synthstrip_v7.3.2_test.pt"

wget -O "${FS_TGZ}" "${FS_URL}"
mkdir -p "${MODEL_DIR}"

tar -xzf "${FS_TGZ}" \
    --strip-components=1 \
    --wildcards \
    "${FS_MODEL_PATH}" \
    -O > "${MODEL_OUT}"

rm -f "${FS_TGZ}"

echo "Model extracted to ${MODEL_OUT}"

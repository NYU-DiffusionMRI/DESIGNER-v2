FROM --platform=linux/amd64 python:3.12.1-bookworm

COPY --from=twom/fsl:6.0 /usr/local/fsl /usr/local/fsl
COPY --from=twom/mrtrix3:dev-latest /usr/local/mrtrix3 /usr/local/mrtrix3
COPY --from=twom/ants:v2.5.4 /usr/local/ants /usr/local/ants


RUN apt-get -qq update \
    && apt-get install -yq --no-install-recommends \
    libgl1-mesa-dev \
    cmake \
    && rm -rf /var/lib/apt/lists/*


ENV FSLDIR=/usr/local/fsl
ENV FSLOUTPUTTYPE=NIFTI_GZ
ENV PATH="${PATH}:/usr/local/fsl/bin:/usr/local/mrtrix3/build/bin:/usr/local/ants/bin"
ENV LD_LIBRARY_PATH="/usr/local/mrtrix3/build/lib:/usr/local/ants/lib"
ENV LD_LIBRARY_PATH="/app/rpg_cpp/fftw-3.3.10/build/lib:${LD_LIBRARY_PATH}"
ENV PYTHONPATH="/usr/local/mrtrix3/python/lib:$PYTHONPATH"

RUN sh -c "$(wget -O- https://github.com/deluan/zsh-in-docker/releases/download/v1.1.5/zsh-in-docker.sh)" -- \
    -p git -p ssh-agent -p 'history-substring-search' \
    -a 'bindkey "\$terminfo[kcuu1]" history-substring-search-up' \
    -a 'bindkey "\$terminfo[kcud1]" history-substring-search-down'





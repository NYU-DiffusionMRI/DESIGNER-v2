FROM --platform=linux/amd64 python:3.12.1-bookworm 

COPY --from=twom/fsl:latest /usr/local/fsl /usr/local/fsl
COPY --from=twom/mrtrix3:dev-latest /usr/local/mrtrix3 /usr/local/mrtrix3


RUN apt-get -qq update \
    && apt-get install -yq --no-install-recommends \
    libgl1-mesa-dev \
    cmake \
    && rm -rf /var/lib/apt/lists/*


ENV FSLDIR=/usr/local/fsl
ENV PATH="${PATH}:/usr/local/fsl/bin:/usr/local/mrtrix3/bin:/root/.local/bin"
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/mrtrix3/src:/usr/local/mrtrix3/core
RUN echo ". /usr/local/fsl/etc/fslconf/fsl.sh" >> /root/.bashrc
ENV PYTHONPATH=/usr/local/mrtrix3/build/lib

RUN sh -c "$(wget -O- https://github.com/deluan/zsh-in-docker/releases/download/v1.1.5/zsh-in-docker.sh)" -- \
    -p git -p ssh-agent -p 'history-substring-search' \
    -a 'bindkey "\$terminfo[kcuu1]" history-substring-search-up' \
    -a 'bindkey "\$terminfo[kcud1]" history-substring-search-down'





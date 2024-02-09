# Use an official Python runtime as a parent image
FROM python:3.12.1-bookworm as base
FROM python:3.12.1-bookworm AS base-builder



FROM base-builder AS mrtrix3-builder

ARG MRTRIX3_GIT_COMMITISH="dev"
RUN apt-get -qq update \
    && apt-get install -yq --no-install-recommends \
    libeigen3-dev \
    libfftw3-dev \
    libgl1-mesa-dev \
    libpng-dev \
    libqt5opengl5-dev \
    libqt5svg5-dev \
    zlib1g-dev \
    ninja-build \
    ccache \
    cmake \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt
RUN git clone https://github.com/MRtrix3/mrtrix3.git
WORKDIR /opt/mrtrix3
RUN git checkout dev
ENV CMAKE_GENERATOR=Ninja
RUN cmake -B build -DCMAKE_INSTALL_PREFIX=/opt/mrtrix3
RUN cmake --build build
RUN cmake --install build


FROM python:3.12.1-bookworm AS fsl-installer

WORKDIR /opt
RUN wget https://fsl.fmrib.ox.ac.uk/fsldownloads/fslconda/releases/fslinstaller.py -O fslinstaller.py -O fslinstaller.py
RUN python fslinstaller.py -d /opt/fsl



FROM python:3.12.1-bookworm as designer-builder

RUN apt-get -qq update \ 
    && apt-get install -yq --no-install-recommends \
    libopenblas-dev \
    libhdf5-dev \
    cmake

COPY requirements.txt .
# Install any needed packages specified in requirements.txt
RUN pip install --user -r requirements.txt



# Set the working directory in the container to /app
WORKDIR /app


# Add the current directory contents into the container at /app
ADD . /app

# Run setup.py
RUN python setup.py install --user


FROM python:3.12.1-bookworm as final
# Install runtime system dependencies.
RUN apt-get -qq update \
    && apt-get install -yq --no-install-recommends \
    binutils \
    dc \
    less \
    libfftw3-dev \
    libgl1-mesa-glx \
    libgomp1 \
    liblapack3 \
    libpng16-16 \
    libqt5core5a \
    libqt5gui5 \
    libqt5network5 \
    libqt5svg5 \
    libqt5widgets5 \
    libquadmath0 \
    python3-distutils \
    && rm -rf /var/lib/apt/lists/*


COPY --from=fsl-installer /opt/fsl /opt/fsl
COPY --from=mrtrix3-builder /opt/mrtrix3 /opt/mrtrix3
COPY --from=designer-builder /root/.local /root/.local

ENV FSLDIR=/usr/local/fsl
ENV PATH="${PATH}:/usr/local/fsl/bin:/opt/mrtrix3/bin:/root/.local/bin"
RUN echo ". /usr/local/fsl/etc/fslconf/fsl.sh" >> /root/.bashrc
ENV PYTHONPATH=/usr/local/mrtrix3/lib

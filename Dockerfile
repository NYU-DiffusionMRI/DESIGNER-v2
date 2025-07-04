# Use an official Python runtime as a parent image
FROM python:3.12.4-bookworm AS base

# Copy dependencies from other images
COPY --from=nyudiffusionmri/fsl:2025-06-16 /usr/local/fsl /usr/local/fsl
COPY --from=nyudiffusionmri/mrtrix3:2025-06-16 /usr/local/mrtrix3/build /usr/local/mrtrix3_build
COPY --from=nyudiffusionmri/ants:2025-06-20 /usr/local/ants /usr/local/ants

# Install common dependencies
RUN apt-get -qq update \
    && apt-get install -yq --no-install-recommends \
    libgl1-mesa-dev \
    cmake \
    make \
    gcc \
    g++ \
    wget \
    hdf5-tools \
    libhdf5-serial-dev \
    liblapack-dev \
    libopenblas-dev \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# FFTW build stage
FROM base AS fftw-builder

WORKDIR /tmp
# Download FFTW (this layer will be cached unless the URL changes)
RUN wget http://www.fftw.org/fftw-3.3.10.tar.gz \
    && tar xzf fftw-3.3.10.tar.gz \
    && cd fftw-3.3.10 \
    && ./configure \
        --enable-shared \
        --disable-static \
        --enable-threads \
        --disable-doc \
        --enable-sse2 \
        --enable-avx \
        CFLAGS="-O3 -fPIC" \
    && make -j$(nproc) \
    && make install \
    && cd / \
    && rm -rf /tmp

# Final stage
FROM base
COPY --from=fftw-builder /usr/local/lib/libfftw* /usr/local/lib/
COPY --from=fftw-builder /usr/local/include/fftw3* /usr/local/include/

# Link and cache FFTW library
RUN ldconfig

WORKDIR /app

COPY requirements.txt .
RUN pip install -r requirements.txt --no-cache-dir

COPY . .

# Run setup.py
RUN python -m pip install .

ENV FSLDIR=/usr/local/fsl
ENV FSLOUTPUTTYPE=NIFTI_GZ
ENV PATH="${PATH}:/usr/local/fsl/bin:/usr/local/mrtrix3_build/bin:/usr/local/ants/bin"
ENV LD_LIBRARY_PATH="/usr/local/mrtrix3_build/src:/usr/local/mrtrix3_build/core:/usr/local/ants/lib:/usr/local/lib"
ENV PYTHONPATH="/usr/local/mrtrix3_build/lib"
RUN echo ". /usr/local/fsl/etc/fslconf/fsl.sh" >> /root/.bashrc

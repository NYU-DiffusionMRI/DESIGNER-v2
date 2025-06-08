# Use an official Python runtime as a parent image
FROM python:3.12.4-bookworm AS base

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
# Copy dependencies from other images
COPY --from=twom/fsl:6.0 /usr/local/fsl /usr/local/fsl
COPY --from=twom/mrtrix3:dev-latest /usr/local/mrtrix3/build /usr/local/mrtrix3_build
COPY --from=twom/ants:v2.5.4 /usr/local/ants /usr/local/ants

# Copy FFTW files from builder stage
COPY --from=fftw-builder /usr/local/lib/libfftw* /usr/local/lib/
COPY --from=fftw-builder /usr/local/include/fftw3* /usr/local/include/

# Verify FFTW installation and set up library cache
RUN ldconfig \
    && ldconfig -p | grep fftw

# Install Python dependencies
COPY requirements.txt .
RUN pip install -r requirements.txt --no-cache-dir

# Set the working directory in the container to /app
WORKDIR /app

# Add the current directory contents into the container at /app
COPY . /app

# Run setup.py
RUN python -m pip install .

ENV FSLDIR=/usr/local/fsl
ENV FSLOUTPUTTYPE=NIFTI_GZ
ENV PATH="${PATH}:/usr/local/fsl/bin:/usr/local/mrtrix3_build/bin:/usr/local/ants/bin"
ENV LD_LIBRARY_PATH="/usr/local/mrtrix3_build/src:/usr/local/mrtrix3_build/core:/usr/local/ants/lib:/usr/local/lib"
ENV PYTHONPATH="/usr/local/mrtrix3_build/lib"
RUN echo ". /usr/local/fsl/etc/fslconf/fsl.sh" >> /root/.bashrc

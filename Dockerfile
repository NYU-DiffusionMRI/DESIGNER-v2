# Use an official Python runtime as a parent image
FROM python:3.12.4-bookworm

# Multi stage build from existing dependencies
COPY --from=twom/fsl:6.0 /usr/local/fsl /usr/local/fsl
COPY --from=twom/mrtrix3:dev-latest /usr/local/mrtrix3/build /usr/local/mrtrix3
COPY --from=antsx/ants:latest /opt/ants /usr/local/ants

# RUN apt-get -qq update \
#     && apt-get install -yq --no-install-recommends \
#     libgl1-mesa-dev \
#     cmake \
#     && rm -rf /var/lib/apt/lists/*

RUN apt-get -qq update && \
    apt-get install -yq --no-install-recommends \
    libgl1-mesa-dev \
    cmake \
    gcc \
    g++ \
    libopenblas-dev \
    liblapack-dev \
    pkg-config \
    libhdf5-serial-dev \
    hdf5-tools && \
    rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
# Install any needed packages specified in requirements.txt
RUN pip install -r requirements.txt --no-cache-dir



# Set the working directory in the container to /app
WORKDIR /app


# Add the current directory contents into the container at /app
ADD . /app

# RUN /app/rpg_cpp/compile.sh 
# Run setup.py
RUN python -m pip install .


ENV FSLDIR=/usr/local/fsl
ENV FSLOUTPUTTYPE=NIFTI_GZ
ENV PATH="${PATH}:/usr/local/fsl/bin:/usr/local/mrtrix3/bin"
ENV LD_LIBRARY_PATH="/usr/local/mrtrix3/src:/usr/local/mrtrix3/core"
ENV LD_LIBRARY_PATH="/app/rpg_cpp/fftw-3.3.10/build/lib:${LD_LIBRARY_PATH}"
RUN echo ". /usr/local/fsl/etc/fslconf/fsl.sh" >> /root/.bashrc
ENV PYTHONPATH=/usr/local/mrtrix3/lib
ENV PATH="/usr/local/ants/bin:$PATH"
ENV LD_LIBRARY_PATH="/usr/local/ants/lib:$LD_LIBRARY_PATH"
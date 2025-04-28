# Use an official Python runtime as a parent image
FROM python:3.12.4-bookworm

# Multi stage build from existing dependencies
COPY --from=twom/fsl:6.0 /usr/local/fsl /usr/local/fsl
COPY --from=twom/mrtrix3:dev-latest /usr/local/mrtrix3/build /usr/local/mrtrix3
COPY --from=twom/ants:v2.5.4 /usr/local/ants /usr/local/ants


RUN apt-get -qq update && \
    apt-get install -yq --no-install-recommends \
    cmake \
    g++ \
    gcc \
    hdf5-tools \
    libgl1-mesa-dev \
    libhdf5-serial-dev \
    liblapack-dev \
    libopenblas-dev \
    pkg-config && \
    rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
# Install any needed packages specified in requirements.txt
RUN pip install -r requirements.txt --no-cache-dir



# Set the working directory in the container to /app
WORKDIR /app


# Add the current directory contents into the container at /app
ADD . /app

# Run setup.py
RUN python -m pip install .


ENV FSLDIR=/usr/local/fsl
ENV FSLOUTPUTTYPE=NIFTI_GZ
ENV PATH="${PATH}:/usr/local/fsl/bin:/usr/local/mrtrix/bin:/usr/local/ants/bin"
ENV LD_LIBRARY_PATH="/usr/local/mrtrix/lib:/usr/local/ants/lib"
ENV LD_LIBRARY_PATH="/app/rpg_cpp/fftw-3.3.10/build/lib:${LD_LIBRARY_PATH}"
ENV PYTHONPATH="/usr/local/mrtrix/lib:$PYTHONPATH"
RUN echo ". /usr/local/fsl/etc/fslconf/fsl.sh" >> /root/.bashrc
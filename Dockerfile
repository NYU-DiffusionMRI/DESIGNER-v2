# Use an official Python runtime as a parent image
FROM twom/designer_deps:latest

# Multi stage build from existing dependencies
COPY --from=twom/fsl:6.0 /usr/local/fsl /usr/local/fsl
COPY --from=twom/mrtrix3:dev-latest /usr/local/mrtrix3/build /usr/local/mrtrix3

RUN apt-get -qq update \
    && apt-get install -yq --no-install-recommends \
    libgl1-mesa-dev \
    cmake \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
# Install any needed packages specified in requirements.txt
RUN pip install -r requirements.txt



# Set the working directory in the container to /app
WORKDIR /app


# Add the current directory contents into the container at /app
ADD . /app

# Run setup.py
RUN python setup.py install

ENV FSLDIR=/usr/local/fsl
ENV PATH="${PATH}:/usr/local/fsl/bin:/usr/local/mrtrix3/bin"
RUN echo ". /usr/local/fsl/etc/fslconf/fsl.sh" >> /root/.bashrc
ENV PYTHONPATH=/usr/local/mrtrix3/lib

# Use an official Python runtime as a parent image
FROM python:3.12.1-bookworm 

COPY --from=twom/fsl:latest /usr/local/fsl /usr/local/fsl
COPY --from=twom/mrtrix3:dev-latest /usr/local/mrtrix3 /usr/local/mrtrix3

RUN apt-get -qq update \
    && apt-get install -yq --no-install-recommends \
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

ENV FSLDIR=/usr/local/fsl
ENV PATH="${PATH}:/usr/local/fsl/bin:/usr/local/mrtrix3/bin:/root/.local/bin"
RUN echo ". /usr/local/fsl/etc/fslconf/fsl.sh" >> /root/.bashrc
ENV PYTHONPATH=/usr/local/mrtrix3/lib

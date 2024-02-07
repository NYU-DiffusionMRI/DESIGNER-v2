# Use an official Python runtime as a parent image
FROM --platform=linux/amd64 python:3.12.1-bookworm

RUN apt-get -y update
RUN apt-get -y install git g++ libeigen3-dev zlib1g-dev libqt5opengl5-dev libqt5svg5-dev libgl1-mesa-dev libfftw3-dev libtiff5-dev libpng-dev libopenblas-dev libhdf5-dev cmake

WORKDIR /usr/local
RUN wget https://fsl.fmrib.ox.ac.uk/fsldownloads/fslconda/releases/fslinstaller.py -O fslinstaller.py
RUN python fslinstaller.py -d /usr/local/fsl

WORKDIR  /usr/local
RUN git clone https://github.com/MRtrix3/mrtrix3.git
WORKDIR /usr/local/mrtrix3
RUN ./configure && ./build


ENV PYTHONPATH=/usr/local/mrtrix3/bin




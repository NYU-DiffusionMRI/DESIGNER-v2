#!/bin/bash

# Compile fftw, nifti libraries and rpg cpp code for command line use
# By Ricardo Coronado-Leija 17-Feb-2023

if [ $# -eq 0 ]; then

echo "Compile third party librares"

# get current path and append rpg_cpp
rpg_path=`pwd`/rpg_cpp

# Entering folder with compressed third party (fftw and nifti) libraries
cd $rpg_path/thirdparty/

# 1) === COMPILING FFTW ===
# Remove folders fftw and fftw-3.3.10 if found
if [ -d fftw-3.3.10 ]; then
rm -rf fftw-3.3.10
fi
if [ -d fftw ]; then
rm -rf fftw
fi
# Uncompress fftw-3.3.10
tar xvzf fftw-3.3.10.tar.gz
# Create folder for fftw compilation
mkdir fftw
# Obtain the absolute path to the folder for fftw compilation
#export FFTW_DIR_NAME=`readlink -f thirdparty/fftw`
export FFTW_DIR_NAME=$(pwd)/fftw
echo "#######*********INSTALLATION FOLDER: ${FFTW_DIR_NAME}*************######"
# move inside fftw folder
cd fftw-3.3.10
# Configure fftw
./configure --prefix=${FFTW_DIR_NAME}
# Install fftw
make
make install
# Exit fftw folder
cd ..


# 2) === COMPILING NIFTI ===
# Remove nifticlib-2.0.0 folder if found
if [ -d nifticlib-2.0.0 ]; then
rm -rf nifticlib-2.0.0
fi
# Uncompress nifticlib-2.0.0
tar xvzf nifticlib-2.0.0.tar.gz
# Remove Makefile using csh
rm nifticlib-2.0.0/Makefile
# Copy Makefile using bash to nifticlib-2.0.0  
cp Makefile nifticlib-2.0.0
# Entering nifti folder
cd nifticlib-2.0.0
# insall nifti
make all
# Exiting nifti folder
cd ..

# Exiting folder with compressed third party (fftw and nifti) libraries
cd ..

fi

echo "Compile RPG"
# 3) === COMPILING RPG ===
if [ -f $rpg_path/unring_rpg.cpp ]; then
rm -rf $rpg_path/rpg
fi
g++ unring_rpg.cpp -o rpg -O3 -Lthirdparty/nifticlib-2.0.0/lib -Lthirdparty/fftw/lib -Ithirdparty/nifticlib-2.0.0/include -Ithirdparty/fftw/include  -lnifticdf -lniftiio -lfftw3 -lznz -lz


chmod -R 0777 $rpg_path/thirdparty/fftw-3.3.10
chmod -R 0777 $rpg_path/thirdparty/nifticlib-2.0.0

rm -rf $rpg_path/thirdparty/fftw
rm -rf $rpg_path/thirdparty/fftw-3.3.10
rm -rf $rpg_path/thirdparty/nifticlib-2.0.0
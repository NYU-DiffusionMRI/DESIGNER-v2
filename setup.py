from setuptools import setup, find_packages, Extension
from glob import glob
import sys, os, stat
import subprocess

with open('requirements.txt') as f:
    requirements = f.readlines()

long_description = '''designer and TMI package for use with brain 
    diffusion MRI processing. Designer is used for image prepocessing 
    including denoising, partial-fourier gibbs correction, epi eddy 
    current and motion, and normalization. TMI is used for dti/dti/wmti/smi
    along with outlier corection.'''

def change_permissions_recursive(path, mode):
    for root, dirs, files in os.walk(path, topdown=False):
        for dir in [os.path.join(root, d) for d in dirs]:
            os.chmod(dir, mode)
    for file in [os.path.join(root, f) for f in files]:
            os.chmod(file, mode)

#change_permissions_recursive('rpg_cpp', 0o777)
#subprocess.run(['./rpg_cpp/compile.sh'], shell=True)

setup(
        name ='designer',
        version ='0.0.1',
        author ='Benjamin Ades-Aron',
        author_email ='benjamin.ades-aron@nyulangone.org',
        url ='https://github.com/badesar1/designer_v2_dev.git',
        description ='Test for designerV2',
        long_description = long_description,
        long_description_content_type ="text/markdown",
        license ='NYU',
        #scripts=['bin/designer', 'bin/tmi'],
        packages = find_packages(),
        data_files=[('constant', glob('constant/*')), ('rpg_cpp', ['rpg_cpp/rpg'])],
        include_package_data=True,
        entry_points = {
            'console_scripts': [
                'designer = src.designer:main',
                'tmi = src.tmi:main'
            ]
        },
        classifiers = [
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
        ],
        keywords ='python diffusion mri mppca gibbs preprocessing dki smi',
        install_requires = requirements,
        zip_safe = False,
        #ext_modules=[module1]
)
    
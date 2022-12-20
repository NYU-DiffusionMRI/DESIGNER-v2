from setuptools import setup, find_packages
from glob import glob
import sys
import subprocess

with open('requirements.txt') as f:
    requirements = f.readlines()

long_description = '''designer and TMI package for use with brain 
    diffusion MRI processing. Designer is used for image prepocessing 
    including denoising, partial-fourier gibbs correction, epi eddy 
    current and motion, and normalization. TMI is used for dti/dti/wmti/smi
    along with outlier corection.'''


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
        #packages = ['bin'],
        #package_data = {'mrtrix3': ['bin/mrtrix3.py']},
        data_files=[('constant', glob('constant/*'))],
        include_package_data=True,
        entry_points = {
            'console_scripts': [
                'designer = bin.designer:main',
                'tmi = bin.tmi:main'
            ]
        },
        classifiers = [
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
        ],
        keywords ='python diffusion mri mppca gibbs preprocessing dki smi',
        install_requires = requirements,
        zip_safe = False
)
    
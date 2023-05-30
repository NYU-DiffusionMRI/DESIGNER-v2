from setuptools import setup, find_packages, Extension
from glob import glob
import os
import subprocess

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

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

change_permissions_recursive('rpg_cpp', 0o777)
subprocess.run(['./rpg_cpp/compile.sh'], shell=True)

setup(
        name ='designer2',
        version ='0.0.9',
        author ='Benjamin Ades-Aron',
        author_email ='benjamin.ades-aron@nyulangone.org',
        url ='https://github.com/badesar1/designer_v2_dev.git',
        description ='Test for designerV2',
        long_description = long_description,
        long_description_content_type ="text/markdown",
        license ='NYU',
        packages = find_packages(),
        data_files=[('constant', glob('constant/*')), ('rpg_cpp', ['rpg_cpp/rpg'])],
        include_package_data=True,
        entry_points = {
            'console_scripts': [
                'designer = designer2.designer:main',
                'tmi = designer2.tmi:main'
            ]
        },
        classifiers = [
            "Programming Language :: Python :: 3",
            "License :: OSI Approved",
            "Operating System :: OS Independent",
        ],
        keywords = 'python diffusion mri mppca gibbs preprocessing dki smi DESIGNER TMI',
        python_requires = '>=3',
        install_requires = ['numpy>=1.23',
                            'scipy>=1.9', 
                            'numpy_groupies>=0.9',
                            'antspyx>=0.3',
                            'scikit-image>=0.19', 
                            'dipy>=1.4', 
                            'tqdm', 
                            'joblib>=1.2', 
                            'cvxpy>=1.2', 
                            'pandas>=1.5'],
        zip_safe = False,
)
    
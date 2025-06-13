import os
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from glob import glob

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

long_description = '''designer and TMI package for use with brain
    diffusion MRI processing. Designer is used for image preprocessing
    including denoising, partial-fourier gibbs correction, epi eddy
    current and motion, and normalization. TMI is used for dti/dki/wmti/smi
    along with outlier correction.'''

class CustomBuildExt(build_ext):
    def build_extension(self, ext):
        # Use system FFTW installation
        ext.include_dirs.append('/usr/local/include')
        ext.library_dirs.append('/usr/local/lib')
        super().build_extension(ext)

conda_prefix = os.environ.get('CONDA_PREFIX')

class get_pybind_include:
    """Helper class to determine the pybind11 include path"""

    def __str__(self):
        import pybind11
        return pybind11.get_include()

# System or Conda FFTW headers
include_path = '/usr/local/include' if not conda_prefix else os.path.join(conda_prefix, 'include')
# System or Conda FFTW libraries
lib_path = '/usr/local/lib' if not conda_prefix else os.path.join(conda_prefix, 'lib')

ext_modules = [
    Extension(
        "lib.rpg",
        ["rpg_cpp/unring_rpg_pybind.cpp"],
        include_dirs=[
            get_pybind_include(),
            include_path,
        ],
        library_dirs=[
            lib_path,
        ],
        libraries=["fftw3", "fftw3_threads", "m"],
        extra_compile_args=["-std=c++11"],
    ),
]

setup(
        name ='designer2',
        version ='2.0.12',
        author ='Benjamin Ades-Aron',
        author_email ='benjamin.ades-aron@nyulangone.org',
        url ='https://github.com/NYU-DiffusionMRI/DESIGNER-v2',
        description ='designerV2',
        long_description = long_description,
        long_description_content_type ="text/markdown",
        license ='NYU',
        packages = find_packages(),
        data_files=[('constant', glob('constant/*'))],
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
        ext_modules = ext_modules,
        cmdclass={"build_ext": CustomBuildExt},
        keywords = 'python diffusion mri mppca gibbs preprocessing dki smi DESIGNER TMI',
        python_requires = '>=3.9',
        install_requires = requirements,
        zip_safe = False,
        verbose = True,
)

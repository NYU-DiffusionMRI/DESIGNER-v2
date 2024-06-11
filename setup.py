import os
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import subprocess
from glob import glob

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

long_description = '''designer and TMI package for use with brain 
    diffusion MRI processing. Designer is used for image preprocessing 
    including denoising, partial-fourier gibbs correction, epi eddy 
    current and motion, and normalization. TMI is used for dti/dki/wmti/smi
    along with outlier correction.'''

class CustomBuildExt(build_ext):
    def run(self):
        self.build_fftw()
        super().run()

    def build_fftw(self):
        fftw_source_dir = os.path.abspath("rpg_cpp/fftw-3.3.10")  # Adjusted path
        fftw_build_dir = os.path.join(fftw_source_dir, "build")

        if not os.path.exists(fftw_build_dir):
            os.makedirs(fftw_build_dir)

        env = os.environ.copy()
        env['CFLAGS'] = '-fPIC'
        subprocess.check_call(
            ["./configure", "--prefix=" + fftw_build_dir, "--enable-shared", "CFLAGS=-fPIC"],
            cwd=fftw_source_dir,
            env=env
        )
        subprocess.check_call(["make"], cwd=fftw_source_dir)
        subprocess.check_call(["make", "install"], cwd=fftw_source_dir)

        self.fftw_include_dir = os.path.join(fftw_build_dir, "include")
        self.fftw_lib_dir = os.path.join(fftw_build_dir, "lib")

    def build_extension(self, ext):
        ext.include_dirs.append(self.fftw_include_dir)
        ext.library_dirs.append(self.fftw_lib_dir)
        super().build_extension(ext)

class get_pybind_include:
    """Helper class to determine the pybind11 include path"""

    def __str__(self):
        import pybind11
        return pybind11.get_include()

ext_modules = [
    Extension(
        "lib.rpg",
        ["rpg_cpp/unring_rpg_pybind.cpp"],
        include_dirs=[
            get_pybind_include(),
            # FFTW include directory will be added during the build
        ],
        library_dirs=[
            # FFTW library directory will be added during the build
        ],
        libraries=["fftw3"],
        extra_compile_args=["-std=c++11"],
    ),
]

setup(
        name ='designer2',
        version ='2.0.8',
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
        python_requires = '>=3',
        install_requires = ['numpy>=1.21.0,<2.0.0',
                            'scipy>=1.9', 
                            'numpy_groupies>=0.9',
                            'antspyx>=0.3',
                            'scikit-image>=0.19', 
                            'dipy>=1.5', 
                            'tqdm', 
                            'joblib>=1.2', 
                            'cvxpy>=1.2', 
                            'pandas>=1.5',
                            'pybind11>=2.12.0'],
        zip_safe = False,
)
    
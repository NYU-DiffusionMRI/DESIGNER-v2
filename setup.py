import os
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

# Dependencies are now defined in pyproject.toml
# This setup.py only handles the C++ extension

class CustomBuildExt(build_ext):
    def build_extension(self, ext):
        # Use system FFTW installation
        ext.include_dirs.append('/usr/local/include')
        ext.library_dirs.append('/usr/local/lib')
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
            '/usr/local/include',  # System FFTW headers
        ],
        library_dirs=[
            '/usr/local/lib',  # System FFTW libraries
        ],
        libraries=["fftw3", "fftw3_threads", "m"],
        extra_compile_args=["-std=c++11"],
    ),
]

setup(
    ext_modules=ext_modules,
    cmdclass={"build_ext": CustomBuildExt},
)

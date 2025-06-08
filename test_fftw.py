#!/usr/bin/env python3

import os
import sys
import ctypes
from ctypes.util import find_library

def test_fftw_library():
    """Test if FFTW library is properly installed and accessible."""
    # Try to find FFTW library
    fftw_lib = find_library('fftw3')
    if not fftw_lib:
        print("ERROR: FFTW library not found!")
        return False

    print(f"Found FFTW library at: {fftw_lib}")

    # Try to load the library
    try:
        _ = ctypes.CDLL(fftw_lib)
        print("Successfully loaded FFTW library")
        return True
    except Exception as e:
        print(f"ERROR: Failed to load FFTW library: {e}")
        return False

def test_python_imports():
    """Test if our package can import and use FFTW."""
    try:
        # Try to import our package
        import lib.rpg
        print("Successfully imported lib.rpg")
        return True
    except ImportError as e:
        print(f"ERROR: Failed to import lib.rpg: {e}")
        return False
    except Exception as e:
        print(f"ERROR: Unexpected error importing lib.rpg: {e}")
        return False

def main():
    print("Testing FFTW installation...")
    print("\n1. Testing FFTW library:")
    lib_ok = test_fftw_library()

    print("\n2. Testing Python imports:")
    import_ok = test_python_imports()

    print("\nSummary:")
    print(f"FFTW library test: {'PASS' if lib_ok else 'FAIL'}")
    print(f"Python import test: {'PASS' if import_ok else 'FAIL'}")

    return 0 if (lib_ok and import_ok) else 1

if __name__ == "__main__":
    sys.exit(main())
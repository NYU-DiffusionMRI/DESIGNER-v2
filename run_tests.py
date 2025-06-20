#!/usr/bin/env python3
"""
Test runner script for DESIGNER-v2 integration tests.

This script provides a convenient way to run tests with different configurations
and generate comprehensive reports.
"""

import argparse
import sys
import subprocess
import os
from pathlib import Path
import time
import json


def run_command(cmd, capture_output=True, check=True):
    """Run a command and return the result."""
    print(f"Running: {' '.join(cmd)}")
    try:
        result = subprocess.run(
            cmd, 
            capture_output=capture_output, 
            text=True, 
            check=check
        )
        return result
    except subprocess.CalledProcessError as e:
        print(f"Command failed with exit code {e.returncode}")
        if e.stdout:
            print("STDOUT:", e.stdout)
        if e.stderr:
            print("STDERR:", e.stderr)
        raise


def install_dependencies():
    """Install test dependencies."""
    print("Installing test dependencies...")
    run_command([sys.executable, "-m", "pip", "install", "-r", "requirements-test.txt"])
    run_command([sys.executable, "-m", "pip", "install", "-e", "."])


def run_unit_tests(verbose=True, coverage=True):
    """Run unit tests."""
    print("\n" + "="*50)
    print("RUNNING UNIT TESTS")
    print("="*50)
    
    cmd = [sys.executable, "-m", "pytest", "tests/test_integration.py"]
    
    if verbose:
        cmd.append("-v")
    
    if coverage:
        cmd.extend([
            "--cov=designer2",
            "--cov=lib", 
            "--cov-report=term-missing",
            "--cov-report=html:htmlcov",
            "--cov-report=xml:coverage.xml"
        ])
    
    cmd.extend(["-m", "not (performance or stress)"])
    
    start_time = time.time()
    result = run_command(cmd, capture_output=False, check=False)
    end_time = time.time()
    
    print(f"\nUnit tests completed in {end_time - start_time:.2f} seconds")
    return result.returncode == 0


def run_performance_tests(verbose=True):
    """Run performance tests."""
    print("\n" + "="*50)
    print("RUNNING PERFORMANCE TESTS")
    print("="*50)
    
    cmd = [
        sys.executable, "-m", "pytest", 
        "tests/test_performance.py",
        "-m", "performance"
    ]
    
    if verbose:
        cmd.append("-v")
    
    start_time = time.time()
    result = run_command(cmd, capture_output=False, check=False)
    end_time = time.time()
    
    print(f"\nPerformance tests completed in {end_time - start_time:.2f} seconds")
    return result.returncode == 0


def run_stress_tests(verbose=True, timeout=300):
    """Run stress tests."""
    print("\n" + "="*50)
    print("RUNNING STRESS TESTS")
    print("="*50)
    
    cmd = [
        sys.executable, "-m", "pytest",
        "tests/test_performance.py",
        "-m", "stress",
        f"--timeout={timeout}"
    ]
    
    if verbose:
        cmd.append("-v")
    
    start_time = time.time()
    result = run_command(cmd, capture_output=False, check=False)
    end_time = time.time()
    
    print(f"\nStress tests completed in {end_time - start_time:.2f} seconds")
    return result.returncode == 0


def run_linting():
    """Run linting checks."""
    print("\n" + "="*50)
    print("RUNNING LINTING CHECKS")
    print("="*50)
    
    success = True
    
    # Black formatting check
    print("Checking code formatting with Black...")
    try:
        run_command([sys.executable, "-m", "black", "--check", "designer2", "lib", "tests"])
        print("✓ Black formatting check passed")
    except subprocess.CalledProcessError:
        print("✗ Black formatting check failed")
        success = False
    
    # Flake8 linting
    print("Running flake8 linting...")
    try:
        run_command([sys.executable, "-m", "flake8", "designer2", "lib", "tests"])
        print("✓ Flake8 linting passed")
    except subprocess.CalledProcessError:
        print("✗ Flake8 linting failed")
        success = False
    
    # MyPy type checking
    print("Running mypy type checking...")
    try:
        run_command([sys.executable, "-m", "mypy", "designer2", "lib", "--ignore-missing-imports"])
        print("✓ MyPy type checking passed")
    except subprocess.CalledProcessError:
        print("✗ MyPy type checking failed")
        success = False
    
    return success


def generate_test_report():
    """Generate a comprehensive test report."""
    print("\n" + "="*50)
    print("GENERATING TEST REPORT")
    print("="*50)
    
    report = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "python_version": sys.version,
        "platform": os.name,
        "tests": {}
    }
    
    # Check for test artifacts
    artifacts = {
        "coverage_xml": "coverage.xml",
        "coverage_html": "htmlcov/index.html", 
        "test_results": "test-results.xml"
    }
    
    for name, path in artifacts.items():
        if os.path.exists(path):
            report["tests"][name] = {"exists": True, "path": path}
            print(f"✓ {name}: {path}")
        else:
            report["tests"][name] = {"exists": False, "path": path}
            print(f"✗ {name}: {path} (not found)")
    
    # Save report
    with open("test-report.json", "w") as f:
        json.dump(report, f, indent=2)
    
    print(f"\nTest report saved to: test-report.json")
    
    if os.path.exists("htmlcov/index.html"):
        print(f"Coverage report available at: file://{os.path.abspath('htmlcov/index.html')}")


def main():
    """Main test runner function."""
    parser = argparse.ArgumentParser(description="DESIGNER-v2 Test Runner")
    parser.add_argument("--install-deps", action="store_true", 
                       help="Install test dependencies before running")
    parser.add_argument("--unit", action="store_true", 
                       help="Run unit tests")
    parser.add_argument("--performance", action="store_true",
                       help="Run performance tests")
    parser.add_argument("--stress", action="store_true",
                       help="Run stress tests")
    parser.add_argument("--lint", action="store_true",
                       help="Run linting checks")
    parser.add_argument("--all", action="store_true",
                       help="Run all tests")
    parser.add_argument("--fast", action="store_true",
                       help="Run only fast unit tests")
    parser.add_argument("--no-coverage", action="store_true",
                       help="Skip coverage reporting")
    parser.add_argument("--timeout", type=int, default=300,
                       help="Timeout for stress tests (seconds)")
    parser.add_argument("--quiet", "-q", action="store_true",
                       help="Reduce output verbosity")
    
    args = parser.parse_args()
    
    # Set defaults if no specific tests requested
    if not any([args.unit, args.performance, args.stress, args.lint, args.all, args.fast]):
        args.fast = True
    
    verbose = not args.quiet
    coverage = not args.no_coverage
    
    success = True
    
    print("DESIGNER-v2 Test Runner")
    print("="*50)
    
    # Install dependencies if requested
    if args.install_deps:
        try:
            install_dependencies()
        except subprocess.CalledProcessError:
            print("Failed to install dependencies")
            return 1
    
    # Run tests
    if args.all or args.fast or args.unit:
        if not run_unit_tests(verbose=verbose, coverage=coverage):
            success = False
    
    if args.all or args.performance:
        if not run_performance_tests(verbose=verbose):
            success = False
    
    if args.all or args.stress:
        if not run_stress_tests(verbose=verbose, timeout=args.timeout):
            success = False
    
    if args.all or args.lint:
        if not run_linting():
            success = False
    
    # Generate report
    generate_test_report()
    
    print("\n" + "="*50)
    if success:
        print("ALL TESTS PASSED ✓")
        return 0
    else:
        print("SOME TESTS FAILED ✗")
        return 1


if __name__ == "__main__":
    sys.exit(main())

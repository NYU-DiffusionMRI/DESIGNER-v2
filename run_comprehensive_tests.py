#!/usr/bin/env python3
"""
Comprehensive test runner for DESIGNER-v2 with support for both mocked and real dependency tests.

This script provides different testing modes:
1. Mock tests: Run with mocked dependencies (fast, no external requirements)
2. Real tests: Run with actual dependencies (slower, requires MRtrix3/FSL/etc.)
3. Coverage tests: Run with coverage reporting
4. Performance tests: Run performance benchmarks

Usage:
    python run_comprehensive_tests.py --help
    python run_comprehensive_tests.py --mock     # Run with mocked dependencies
    python run_comprehensive_tests.py --real     # Run with real dependencies
    python run_comprehensive_tests.py --all      # Run all tests
"""

import argparse
import sys
import subprocess
import os
from pathlib import Path
import time
import json


def check_dependency(cmd: str) -> bool:
    """Check if a command-line dependency is available."""
    try:
        result = subprocess.run([cmd, '--help'], 
                              capture_output=True, 
                              timeout=10)
        return result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return False


def check_dependencies():
    """Check availability of external dependencies."""
    deps = {
        'mrtrix3': check_dependency('mrinfo'),
        'fsl': check_dependency('fslinfo'),
        'ants': check_dependency('antsRegistration'),
        'python': check_dependency('python'),
        'pytest': check_dependency('pytest')
    }
    return deps


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


def run_mock_tests(verbose=True, coverage=True):
    """Run tests with mocked dependencies."""
    print("\n" + "="*60)
    print("RUNNING MOCK TESTS (NO EXTERNAL DEPENDENCIES REQUIRED)")
    print("="*60)
    
    cmd = [sys.executable, "-m", "pytest", "tests/test_integration.py", "tests/test_units.py"]
    
    if verbose:
        cmd.append("-v")
    
    if coverage:
        cmd.extend([
            "--cov=designer2",
            "--cov=lib", 
            "--cov-report=term-missing",
            "--cov-report=html:htmlcov_mock",
            "--cov-report=xml:coverage_mock.xml"
        ])
    
    # Exclude real dependency tests
    cmd.extend(["-m", "not real_deps"])
    
    start_time = time.time()
    result = run_command(cmd, capture_output=False, check=False)
    end_time = time.time()
    
    print(f"\nMock tests completed in {end_time - start_time:.2f} seconds")
    return result.returncode == 0


def run_real_tests(verbose=True, coverage=True):
    """Run tests with real dependencies."""
    print("\n" + "="*60)
    print("RUNNING REAL INTEGRATION TESTS (REQUIRES EXTERNAL DEPENDENCIES)")
    print("="*60)
    
    # Check dependencies first
    deps = check_dependencies()
    missing_deps = [name for name, available in deps.items() if not available and name != 'ants']
    
    if missing_deps:
        print(f"WARNING: Missing dependencies: {missing_deps}")
        print("Some tests may be skipped.")
    
    cmd = [sys.executable, "-m", "pytest", "tests/test_real_integration.py"]
    
    if verbose:
        cmd.append("-v")
    
    if coverage:
        cmd.extend([
            "--cov=designer2",
            "--cov=lib", 
            "--cov-report=term-missing",
            "--cov-report=html:htmlcov_real",
            "--cov-report=xml:coverage_real.xml"
        ])
    
    # Include only real dependency tests
    cmd.extend(["-m", "real_deps"])
    
    start_time = time.time()
    result = run_command(cmd, capture_output=False, check=False)
    end_time = time.time()
    
    print(f"\nReal integration tests completed in {end_time - start_time:.2f} seconds")
    return result.returncode == 0


def run_performance_tests(verbose=True):
    """Run performance tests."""
    print("\n" + "="*60)
    print("RUNNING PERFORMANCE TESTS")
    print("="*60)
    
    cmd = [
        sys.executable, "-m", "pytest", 
        "tests/test_performance.py",
        "tests/test_real_integration.py",
        "-m", "performance"
    ]
    
    if verbose:
        cmd.append("-v")
    
    start_time = time.time()
    result = run_command(cmd, capture_output=False, check=False)
    end_time = time.time()
    
    print(f"\nPerformance tests completed in {end_time - start_time:.2f} seconds")
    return result.returncode == 0


def run_linting():
    """Run linting checks."""
    print("\n" + "="*60)
    print("RUNNING LINTING CHECKS")
    print("="*60)
    
    success = True
    
    # Check if linting tools are available
    linting_tools = ['black', 'flake8']
    available_tools = [tool for tool in linting_tools if check_dependency(tool)]
    
    if not available_tools:
        print("No linting tools available. Installing...")
        try:
            run_command([sys.executable, "-m", "pip", "install", "black", "flake8"])
        except subprocess.CalledProcessError:
            print("Failed to install linting tools. Skipping linting.")
            return True
    
    # Black formatting check
    print("Checking code formatting with Black...")
    try:
        run_command([sys.executable, "-m", "black", "--check", "designer2", "lib", "tests"])
        print("✓ Black formatting check passed")
    except subprocess.CalledProcessError:
        print("✗ Black formatting check failed")
        print("Run 'black designer2 lib tests' to fix formatting")
        success = False
    
    # Flake8 linting
    print("Running flake8 linting...")
    try:
        run_command([sys.executable, "-m", "flake8", "designer2", "lib", "tests", 
                    "--max-line-length=88", "--extend-ignore=E203,W503"])
        print("✓ Flake8 linting passed")
    except subprocess.CalledProcessError:
        print("✗ Flake8 linting failed")
        success = False
    
    return success


def generate_comprehensive_report():
    """Generate a comprehensive test report."""
    print("\n" + "="*60)
    print("GENERATING COMPREHENSIVE TEST REPORT")
    print("="*60)
    
    report = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "python_version": sys.version,
        "platform": os.name,
        "dependencies": check_dependencies(),
        "test_artifacts": {}
    }
    
    # Check for test artifacts
    artifacts = {
        "mock_coverage_xml": "coverage_mock.xml",
        "mock_coverage_html": "htmlcov_mock/index.html",
        "real_coverage_xml": "coverage_real.xml", 
        "real_coverage_html": "htmlcov_real/index.html",
        "test_results": "test-results.xml"
    }
    
    for name, path in artifacts.items():
        if os.path.exists(path):
            report["test_artifacts"][name] = {"exists": True, "path": path}
            print(f"✓ {name}: {path}")
        else:
            report["test_artifacts"][name] = {"exists": False, "path": path}
            print(f"✗ {name}: {path} (not found)")
    
    # Save report
    with open("comprehensive-test-report.json", "w") as f:
        json.dump(report, f, indent=2)
    
    print(f"\nComprehensive test report saved to: comprehensive-test-report.json")
    
    # Print summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)
    print(f"Dependencies available: {sum(report['dependencies'].values())}/{len(report['dependencies'])}")
    print(f"Test artifacts generated: {sum(1 for a in report['test_artifacts'].values() if a['exists'])}/{len(report['test_artifacts'])}")
    
    if os.path.exists("htmlcov_mock/index.html"):
        print(f"Mock test coverage: file://{os.path.abspath('htmlcov_mock/index.html')}")
    if os.path.exists("htmlcov_real/index.html"):
        print(f"Real test coverage: file://{os.path.abspath('htmlcov_real/index.html')}")


def main():
    """Main test runner function."""
    parser = argparse.ArgumentParser(description="DESIGNER-v2 Comprehensive Test Runner")
    parser.add_argument("--install-deps", action="store_true", 
                       help="Install test dependencies before running")
    parser.add_argument("--mock", action="store_true", 
                       help="Run tests with mocked dependencies")
    parser.add_argument("--real", action="store_true",
                       help="Run tests with real dependencies")
    parser.add_argument("--performance", action="store_true",
                       help="Run performance tests")
    parser.add_argument("--lint", action="store_true",
                       help="Run linting checks")
    parser.add_argument("--all", action="store_true",
                       help="Run all available tests")
    parser.add_argument("--no-coverage", action="store_true",
                       help="Skip coverage reporting")
    parser.add_argument("--quiet", "-q", action="store_true",
                       help="Reduce output verbosity")
    parser.add_argument("--check-deps", action="store_true",
                       help="Only check dependencies and exit")
    
    args = parser.parse_args()
    
    # Set defaults if no specific tests requested
    if not any([args.mock, args.real, args.performance, args.lint, args.all]):
        args.mock = True  # Default to mock tests
    
    verbose = not args.quiet
    coverage = not args.no_coverage
    
    success = True
    
    print("DESIGNER-v2 Comprehensive Test Runner")
    print("="*60)
    
    # Check dependencies
    if args.check_deps:
        deps = check_dependencies()
        print("Dependency Check Results:")
        for name, available in deps.items():
            status = "✓" if available else "✗"
            print(f"  {status} {name}")
        return 0 if all(deps.values()) else 1
    
    # Install dependencies if requested
    if args.install_deps:
        try:
            install_dependencies()
        except subprocess.CalledProcessError:
            print("Failed to install dependencies")
            return 1
    
    # Run tests
    if args.all or args.mock:
        if not run_mock_tests(verbose=verbose, coverage=coverage):
            success = False
    
    if args.all or args.real:
        if not run_real_tests(verbose=verbose, coverage=coverage):
            success = False
    
    if args.all or args.performance:
        if not run_performance_tests(verbose=verbose):
            success = False
    
    if args.all or args.lint:
        if not run_linting():
            success = False
    
    # Generate comprehensive report
    generate_comprehensive_report()
    
    print("\n" + "="*60)
    if success:
        print("ALL TESTS PASSED ✓")
        return 0
    else:
        print("SOME TESTS FAILED ✗")
        return 1


if __name__ == "__main__":
    sys.exit(main())

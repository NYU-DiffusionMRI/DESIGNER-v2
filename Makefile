# Makefile for DESIGNER-v2 testing and development

.PHONY: help install install-dev test test-fast test-integration test-performance test-stress coverage lint format clean docs

help:  ## Show this help message
	@echo "Available commands:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

install:  ## Install package and dependencies
	pip install -e .

install-dev:  ## Install package with development dependencies
	pip install -e .
	pip install -r requirements-test.txt

test:  ## Run all tests
	pytest tests/ -v

test-fast:  ## Run fast unit tests only
	pytest tests/test_integration.py -v -m "not (performance or stress or slow)"

test-integration:  ## Run integration tests
	pytest tests/test_integration.py -v -m "integration"

test-performance:  ## Run performance tests
	pytest tests/test_performance.py -v -m "performance"

test-stress:  ## Run stress tests
	pytest tests/test_performance.py -v -m "stress" --timeout=300

coverage:  ## Run tests with coverage report
	pytest tests/ -v --cov=designer2 --cov=lib --cov-report=term-missing --cov-report=html --cov-report=xml
	@echo "Coverage report generated in htmlcov/index.html"

coverage-upload:  ## Upload coverage to codecov
	codecov

lint:  ## Run linting checks
	black --check designer2 lib tests
	flake8 designer2 lib tests
	mypy designer2 lib --ignore-missing-imports

format:  ## Format code with black
	black designer2 lib tests
	isort designer2 lib tests

security:  ## Run security checks
	bandit -r designer2 lib

clean:  ## Clean build artifacts and cache
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf htmlcov/
	rm -rf .coverage
	rm -rf coverage.xml
	rm -rf test-results.xml
	rm -rf .pytest_cache/
	rm -rf .tox/
	find . -type d -name __pycache__ -delete
	find . -type f -name "*.pyc" -delete

docs:  ## Build documentation
	cd docs && sphinx-build -W -b html . _build/html

install-hooks:  ## Install pre-commit hooks
	pre-commit install

test-all:  ## Run all tests including performance and stress tests
	pytest tests/ -v --cov=designer2 --cov=lib --cov-report=term-missing

ci:  ## Run CI pipeline locally
	make lint
	make test-fast
	make coverage

benchmark:  ## Run benchmark tests
	pytest tests/test_performance.py -v -m "performance" --benchmark-only

# Environment setup
setup-dev:  ## Set up development environment
	pip install -e .
	pip install -r requirements-test.txt
	pip install pre-commit
	pre-commit install
	@echo "Development environment setup complete!"

# Release helpers
check-release:  ## Check if package is ready for release
	python setup.py check --strict --metadata
	twine check dist/*

build:  ## Build distribution packages
	python -m build

upload-test:  ## Upload to test PyPI
	twine upload --repository testpypi dist/*

upload:  ## Upload to PyPI
	twine upload dist/*

# Docker helpers
docker-test:  ## Run tests in Docker container
	docker build -t designer2-test .
	docker run --rm designer2-test pytest tests/ -v

# Parallel testing
test-parallel:  ## Run tests in parallel
	pytest tests/ -v -n auto

# Test with specific Python version using tox
test-py38:  ## Test with Python 3.8
	tox -e py38

test-py39:  ## Test with Python 3.9
	tox -e py39

test-py310:  ## Test with Python 3.10
	tox -e py310

test-py311:  ## Test with Python 3.11
	tox -e py311

# Generate test data
generate-test-data:  ## Generate synthetic test data
	python -c "from tests.conftest import *; import tempfile; import shutil; d = tempfile.mkdtemp(); print(f'Test data in: {d}')"

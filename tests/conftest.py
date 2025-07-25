import pytest
import yaml

def pytest_addoption(parser):
    parser.addoption(
        "--no-cleanup",
        action="store_true",
        default=False,
        help="Skip cleanup of DESIGNER's scratch (processing) directory and TMI's params directory.",
    )


@pytest.fixture(scope="session")
def tolerance_config() -> dict:
    with open("tests/tolerance_config.yaml") as f:
        config = yaml.safe_load(f)
    
    return config

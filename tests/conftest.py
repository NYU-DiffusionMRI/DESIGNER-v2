import pytest
import numpy as np
import random

# @pytest.fixture(autouse=True)
# def set_random_seed():
#     SEED = 42
#     random.seed(SEED)
#     np.random.seed(SEED)

@pytest.fixture(scope="session")
def cli_seed():
    return "42"


import shutil
from pathlib import Path

@pytest.fixture(scope="session", autouse=True)
def cleanup_tmp_dir():
    yield
    tmp_dir = Path("tests/tmp")
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)

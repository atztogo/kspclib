import os
import pytest
import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture(scope='session')
def return_None():
    return None

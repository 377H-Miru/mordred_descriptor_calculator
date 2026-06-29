import os
import sys

# Ensure src layout is on sys.path for pytest discovery in all CI environments
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))

import pytest
from mordred_descriptor_calculator.compat import patch_numpy_for_mordred_compat

patch_numpy_for_mordred_compat()

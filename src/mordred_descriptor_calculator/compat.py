import sys
import numpy as np

def patch_numpy_for_mordred_compat():
    """
    Mordred may reference deprecated NumPy aliases in some environments.
    This compatibility shim maps np.product to np.prod, np.float to float, etc.
    """
    if not hasattr(np, 'product'):
        np.product = np.prod
    if not hasattr(np, 'float'):
        np.float = float
    if not hasattr(np, 'int'):
        np.int = int
    if not hasattr(np, 'bool'):
        np.bool = bool
    if not hasattr(np, 'typeDict'):
        np.typeDict = np.sctypeDict
    if 'numpy.core.fromnumeric' in sys.modules:
        try:
            setattr(sys.modules['numpy.core.fromnumeric'], 'product', np.prod)
        except Exception:
            pass

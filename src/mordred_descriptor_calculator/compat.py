import sys
import numpy as np

def patch_numpy_for_mordred_compat():
    """
    Mordred may reference deprecated NumPy aliases in some environments.
    This compatibility shim maps np.product to np.prod when needed.
    """
    if not hasattr(np, 'product'):
        np.product = np.prod
    if 'numpy.core.fromnumeric' in sys.modules:
        try:
            setattr(sys.modules['numpy.core.fromnumeric'], 'product', np.prod)
        except Exception:
            pass

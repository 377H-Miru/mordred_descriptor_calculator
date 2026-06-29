import sys
import numpy as np

def patch_numpy_for_mordred_compat():
    """
    Mordred references deprecated NumPy aliases (specifically np.product and np.typeDict) removed in modern NumPy versions.
    This compatibility shim maps missing aliases to their modern equivalents safely without corrupting Pandas C-extensions.
    """
    if not hasattr(np, 'product'):
        try:
            setattr(np, 'product', np.prod)
        except Exception:
            pass
    if not hasattr(np, 'cumproduct'):
        try:
            setattr(np, 'cumproduct', np.cumprod)
        except Exception:
            pass
    if not hasattr(np, 'typeDict'):
        try:
            setattr(np, 'typeDict', np.sctypeDict)
        except Exception:
            pass

    if 'numpy.core.fromnumeric' in sys.modules:
        try:
            setattr(sys.modules['numpy.core.fromnumeric'], 'product', getattr(np, 'prod', None))
        except Exception:
            pass

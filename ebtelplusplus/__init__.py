"""
ebtelplusplus: Zero-dimensional hydrodynamics of coronal loops
"""
from ebtelplusplus.high_level import EbtelResult, run

try:
    from ebtelplusplus._version import __version__
except ImportError:
    __version__ = "unknown"

__all__ = ["run", "EbtelResult", "__version__"]

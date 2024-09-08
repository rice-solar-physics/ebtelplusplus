"""
ebtelplusplus: Zero-dimensional hydrodynamics of coronal loops
"""
from ebtelplusplus.high_level import build_configuration, run

try:
    from ebtelplusplus._version import __version__
except ImportError:
    __version__ = "unknown"

__all__ = ["run", "build_configuration", "__version__"]

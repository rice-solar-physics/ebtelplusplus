"""
ebtelplusplus: Zero-dimensional hydrodynamics of coronal loops
"""
from ebtelplusplus.high_level import run, build_configuration

try:
    from ebtelplusplus.version import __version__
except ImportError:
    __version__ = "unknown"

__all__ = ["run", "build_configuration", "__version__"]

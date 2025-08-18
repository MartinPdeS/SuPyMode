from .utils import load_superset  # noqa: F401, W292


try:
    from ._version import version as __version__  # noqa: F401

except ImportError:
    __version__ = "0.0.0"

# -

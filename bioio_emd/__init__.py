"""BioIO reader and writer plugin for EMD images."""

try:
    from ._version import __version__  # type: ignore
except Exception:
    __version__ = "0.0.0"

from .reader import Reader
from .reader_metadata import ReaderMetadata

__all__ = ["Reader", "ReaderMetadata"]

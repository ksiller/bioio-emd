"""BioIO reader and writer plugin for EMD images."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("bioio-emd")
except PackageNotFoundError:
    __version__ = "uninstalled"

from .reader import Reader
from .reader_metadata import ReaderMetadata

__all__ = ["Reader", "ReaderMetadata"]

from typing import List

import bioio_base.reader
import bioio_base.reader_metadata


class ReaderMetadata(bioio_base.reader_metadata.ReaderMetadata):
    """Metadata for the EMD reader plugin."""

    @staticmethod
    def get_supported_extensions() -> List[str]:
        """Return the file extensions this plugin reads."""
        return [".emd"]

    @staticmethod
    def get_reader() -> bioio_base.reader.Reader:
        """Return the Reader class for this plugin."""
        from .reader import Reader

        return Reader

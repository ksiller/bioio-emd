from typing import Any

from bioio_base import dimensions, types
from bioio_base.writer import Writer


class EmdWriter(Writer):
    """Write image data to EMD files."""

    @staticmethod
    def save(
        data: types.ArrayLike,
        uri: types.PathLike,
        dim_order: str = dimensions.DEFAULT_DIMENSION_ORDER,
        **kwargs: Any,
    ) -> None:
        """Write an array to an EMD file."""
        raise NotImplementedError()

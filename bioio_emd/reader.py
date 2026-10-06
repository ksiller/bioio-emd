from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import dask.array as da
import h5py
import numpy as np
import rsciio.emd
import xarray as xr
from bioio_base import constants
from bioio_base.dimensions import DimensionNames
from bioio_base.io import pathlike_to_fs
from bioio_base.reader import Reader as BaseReader
from bioio_base.types import PhysicalPixelSizes
from fsspec.spec import AbstractFileSystem

_STANDARD_DIMS = {
    "t": DimensionNames.Time,
    "c": DimensionNames.Channel,
    "z": DimensionNames.SpatialZ,
    "thickness": DimensionNames.SpatialZ,
    "y": DimensionNames.SpatialY,
    "x": DimensionNames.SpatialX,
}
_TARGET_DIMS = (
    DimensionNames.Time,
    DimensionNames.Channel,
    DimensionNames.SpatialZ,
    DimensionNames.SpatialY,
    DimensionNames.SpatialX,
)


class Reader(BaseReader):
    """Read EMD files and expose them as TCZYX arrays."""

    def __init__(
        self,
        image: Any,
        fs_kwargs: Optional[Dict[str, Any]] = None,
        **kwargs: Any,
    ) -> None:
        """Open an EMD file and index its datasets."""
        self._fs, self._path = pathlike_to_fs(
            image,
            enforce_exists=True,
            fs_kwargs=fs_kwargs or {},
        )
        self._h5_file: Optional[h5py.File] = None
        self._raw_datasets = self._load_datasets()
        self._scenes = tuple(
            str(dataset.get("name") or f"Dataset_{index}")
            for index, dataset in enumerate(self._raw_datasets)
        )
        self._current_scene_index = 0

    @staticmethod
    def _is_supported_image(
        fs: AbstractFileSystem,
        path: str,
        **kwargs: Any,
    ) -> bool:
        """Return whether the path has an .emd suffix."""
        return Path(path).suffix.lower() == ".emd"

    @property
    def scenes(self) -> Tuple[str, ...]:
        """Return the dataset names stored in the file."""
        return self._scenes

    @property
    def physical_pixel_sizes(self) -> PhysicalPixelSizes:
        """Return Z, Y, and X pixel scales in micrometers."""
        axes = self._current_dataset().get("axes", [])
        return PhysicalPixelSizes(
            self._scale_for(axes, "z", "thickness"),
            self._scale_for(axes, "y"),
            self._scale_for(axes, "x"),
        )

    def _read_delayed(self) -> xr.DataArray:
        """Return the current scene as a dask-backed TCZYX array."""
        image = self._as_tczyx()
        if isinstance(image.data, da.Array):
            return image
        return image.copy(data=da.from_array(np.asarray(image.data)))

    def _read_immediate(self) -> xr.DataArray:
        """Return the current scene as an in-memory TCZYX array."""
        image = self._as_tczyx()
        data = image.data.compute() if isinstance(image.data, da.Array) else image.data
        return image.copy(data=np.asarray(data))

    def _load_datasets(self) -> List[Dict[str, Any]]:
        """Load datasets with RosettaSciIO, falling back to classic EMD."""
        try:
            datasets = rsciio.emd.file_reader(self._path, lazy=True)
        except OSError as error:
            if "not a supported EMD file" not in str(error):
                raise
            datasets = self._read_classic_emd()
        if not datasets:
            raise ValueError(f"No datasets found in {self._path}.")
        return datasets

    def _read_classic_emd(self) -> List[Dict[str, Any]]:
        """Read Berkeley EMD files that store version as a single root attribute."""
        self._h5_file = h5py.File(self._path, "r")
        data_group = self._h5_file.get("data")
        if not isinstance(data_group, h5py.Group):
            raise ValueError(f"{self._path} has no EMD data group.")

        datasets = []
        for name, group in data_group.items():
            if not isinstance(group, h5py.Group) or "data" not in group:
                continue
            pixels = group["data"]
            axes = [
                self._classic_axis(group.get(f"dim{index + 1}"), size)
                for index, size in enumerate(pixels.shape)
            ]
            datasets.append(
                {
                    "name": name,
                    "data": da.from_array(pixels, chunks="auto"),
                    "axes": axes,
                    "metadata": self._classic_metadata(name),
                }
            )
        return datasets

    def _classic_axis(
        self,
        dataset: Optional[h5py.Dataset],
        size: int,
    ) -> Dict[str, Any]:
        """Build an axis description from a classic EMD dimension dataset."""
        if dataset is None:
            return {"name": "", "size": size, "navigate": True, "scale": 1.0}

        values = dataset[...]
        name = _as_text(dataset.attrs.get("name", ""))
        units = _normalize_emd_units(_as_text(dataset.attrs.get("units", "")))
        if values.dtype.kind in {"S", "U", "O"}:
            return {
                "name": name,
                "units": units,
                "size": size,
                "navigate": False,
                "scale": 1.0,
                "labels": [_as_text(value) for value in values],
            }

        coords = np.asarray(values, dtype=float).reshape(-1)
        if coords.size > 1:
            scale = float(np.median(np.diff(coords)))
            offset = float(coords[0])
        else:
            scale = 1.0
            offset = float(coords[0]) if coords.size else 0.0
        return {
            "name": name,
            "units": units,
            "size": size,
            "navigate": True,
            "scale": scale,
            "offset": offset,
        }

    def _classic_metadata(self, dataset_name: str) -> Dict[str, Any]:
        """Collect version and microscope, sample, and user attributes."""
        assert self._h5_file is not None
        metadata: Dict[str, Any] = {"dataset": dataset_name}
        version = self._h5_file.attrs.get("version")
        if version is not None:
            metadata["version"] = _as_text(version)
        for group_name in ("microscope", "sample", "user"):
            group = self._h5_file.get(group_name)
            if isinstance(group, h5py.Group):
                metadata[group_name] = {
                    key: _as_text(group.attrs[key]) for key in group.attrs
                }
        return metadata

    def _current_dataset(self) -> Dict[str, Any]:
        """Return the dictionary for the active scene."""
        try:
            return self._raw_datasets[self._current_scene_index]
        except IndexError as error:
            raise IndexError(
                f"Scene index {self._current_scene_index} is not present in "
                f"{self._path}."
            ) from error

    def _as_tczyx(self) -> xr.DataArray:
        """Map the current dataset onto a TCZYX array."""
        dataset = self._current_dataset()
        data = dataset["data"]
        axes = list(dataset.get("axes", []))
        source_names = self._unique_axis_names(axes, data.ndim)
        rename, channel_axis = self._map_axes(source_names, axes)

        image = xr.DataArray(data, dims=source_names).rename(rename)
        for dim in _TARGET_DIMS:
            if dim not in image.dims:
                image = image.expand_dims({dim: 1})
        image = image.transpose(*_TARGET_DIMS)

        channel_size = int(image.sizes[DimensionNames.Channel])
        image = image.assign_coords(
            {DimensionNames.Channel: self._channel_names(channel_axis, channel_size)}
        )
        image.attrs[constants.METADATA_UNPROCESSED] = dataset.get("metadata", {})
        return image

    @staticmethod
    def _unique_axis_names(axes: List[Dict[str, Any]], ndim: int) -> List[str]:
        """Return unique axis names, one per data dimension."""
        raw_names = []
        for index in range(ndim):
            if index < len(axes) and axes[index].get("name"):
                raw_names.append(str(axes[index]["name"]))
            else:
                raw_names.append(f"dim{index}")

        seen: Dict[str, int] = {}
        names = []
        for name in raw_names:
            count = seen.get(name.lower(), 0)
            seen[name.lower()] = count + 1
            names.append(name if count == 0 else f"{name}_{count}")
        return names

    @staticmethod
    def _map_axes(
        source_names: List[str],
        axes: List[Dict[str, Any]],
    ) -> Tuple[Dict[str, str], Optional[Dict[str, Any]]]:
        """Map source axis names onto TCZYX and identify the channel axis."""
        rename: Dict[str, str] = {}
        used = set()
        channel_axis: Optional[Dict[str, Any]] = None

        for index, name in enumerate(source_names):
            target = _STANDARD_DIMS.get(name.lower())
            if target is not None and target not in used:
                rename[name] = target
                used.add(target)
                if target == DimensionNames.Channel and index < len(axes):
                    channel_axis = axes[index]

        for index, name in enumerate(source_names):
            if name in rename:
                continue
            axis = axes[index] if index < len(axes) else {}
            prefers_channel = axis.get("navigate", True) is False
            if prefers_channel and DimensionNames.Channel not in used:
                rename[name] = DimensionNames.Channel
                used.add(DimensionNames.Channel)
                channel_axis = axis
                continue
            for target in _TARGET_DIMS:
                if target not in used:
                    rename[name] = target
                    used.add(target)
                    break
            else:
                raise ValueError(
                    "EMD data has more than 5 dimensions and cannot be stored as TCZYX. "
                    f"Axes: {source_names}."
                )

        return rename, channel_axis

    @staticmethod
    def _channel_names(axis: Optional[Dict[str, Any]], size: int) -> List[str]:
        """Return labels for the channel coordinate."""
        labels = None if axis is None else axis.get("labels")
        if labels and len(labels) == size:
            return [str(label) for label in labels]
        if axis is None:
            if size == 1:
                return ["EM_Signal"]
            return [f"channel_{index}" for index in range(size)]
        axis_name = axis.get("name", "signal")
        return [f"{axis_name}_channel_{index}" for index in range(size)]

    @staticmethod
    def _axis_named(axes: List[Dict[str, Any]], name: str) -> Dict[str, Any]:
        """Return the axis dictionary with the given name."""
        for axis in axes:
            if str(axis.get("name", "")).lower() == name:
                return axis
        return {}

    def _scale_for(self, axes: List[Dict[str, Any]], *names: str) -> Optional[float]:
        """Return the first matching axis scale in micrometers."""
        for name in names:
            axis = self._axis_named(axes, name)
            if axis:
                return self._scale_in_micrometers(axis)
        return None

    @staticmethod
    def _scale_in_micrometers(axis: Dict[str, Any]) -> Optional[float]:
        """Convert an axis scale to micrometers."""
        if not axis:
            return None
        scale = float(axis.get("scale", 1.0))
        unit = str(axis.get("units", axis.get("unit", ""))).lower()
        if unit in {"m", "meter", "meters"}:
            return scale * 1e6
        if unit in {"nm", "nanometer", "nanometers"}:
            return scale * 1e-3
        if unit in {"mm", "millimeter", "millimeters"}:
            return scale * 1e3
        return scale


def _as_text(value: Any) -> str:
    """Decode bytes and NumPy values to text."""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.ndarray):
        if value.shape == ():
            return _as_text(value.item())
        return ",".join(_as_text(item) for item in value)
    if isinstance(value, np.generic):
        return _as_text(value.item())
    return str(value)


def _normalize_emd_units(units: str) -> str:
    """Convert bracketed EMD unit strings such as ``[n_m]`` to ``nm``."""
    text = units.strip().strip("[]")
    for prefix, replacement in (("n_", "n"), ("u_", "u"), ("p_", "p"), ("k_", "k")):
        if text.startswith(prefix):
            return replacement + text[len(prefix) :]
    return text

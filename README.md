# bioio-emd

A BioIO plugin for reading and writing `.emd` images. It follows the [BioIO](https://github.com/bioio-devs/bioio) reader and writer plugin APIs so files can be opened through `bioio.BioImage` alongside other format plugins.

## Reading .emd files

```python
from bioio import BioImage
import bioio_emd

img = BioImage("file.emd", reader=bioio_emd.Reader)
img.data
```

## Accessing metadata

Source axes are mapped onto `TCZYX`. Known names such as `x`, `y`, `z`, and `thickness` keep their spatial axes. A non-navigation axis, such as `complex`, becomes the channel axis. Remaining axes fill unused slots in the order `T`, `C`, `Z`, `Y`, `X`.

```python
img.dims.order  # "TCZYX"
img.dims.T, img.dims.C, img.dims.Z, img.dims.Y, img.dims.X
img.channel_names
img.physical_pixel_sizes  # Z, Y, X in micrometers
```

`img.metadata` is the file metadata dictionary. Classic Berkeley EMD files include the root version and the `microscope`, `sample`, and `user` groups:

```python
img.metadata["sample"]["material"]
img.metadata["microscope"]["voltage"]
```

# Writing .emd files

```python
from bioio.writers import EmdWriter

EmdWriter.save(image, "file.emd", dim_order="ZYX")
```

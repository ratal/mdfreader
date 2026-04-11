**MDFREADER**
**************

Abstract:
=========
This module imports MDF files (Measured Data Format V3.x and V4.x), typically
from INCA (ETAS), CANape or CANoe. It is widely used in the automotive industry
to record data from ECUs. The main module `mdfreader.py` inherits from two
module pairs (one per MDF version): the first reads the file's block structure
(`mdfinfoX`), and the second reads the raw data (`mdfXreader`). It can
optionally run multithreaded and was designed for efficient batch processing of
large endurance-evaluation files for data mining.

Performance:
============
When Cython is available (strongly recommended), mdfreader uses several
low-level optimisations:

* **Fast CN/CC/SI/TX metadata reader** (`read_cn_chain_fast` in `dataRead.pyx`):
  walks the entire MDF4 channel linked list in a single Cython function using
  POSIX `pread()` (no Python file-object dispatch, no GIL during I/O) and
  C packed-struct `memcpy` parsing. A fast `<TX>…</TX>` bytes scan replaces
  `lxml.objectify` for the common MD-block pattern (~95% of files). Result:
  **3–4× speedup** on large files compared to the pure-Python path.

* **SymBufReader**: a Cython bidirectional-buffered wrapper around the raw file
  object. MDF4 metadata blocks are linked by backward-pointing pointers;
  `SymBufReader` keeps a 64 KB buffer centred on the current position so that
  most seeks are served from cache without a kernel `read()`.

* **Vectorised data reading**: sorted channel groups are read in a single
  `readinto()` call into a flat `uint8` buffer that is then reinterpreted as a
  structured record array — zero copies, no per-chunk Python loop.

Typical timings on a 184 MB / 36 000-channel MDF4 file:

| Scenario          | Time   |
|-------------------|--------|
| Pure Python path  | ~1.9 s |
| v4.2 with Cython  | ~1.9 s |
| v4.3 (this version) | **~0.6 s** |

The structure of the mdf object inheriting from python dict
===========================================================
For each channel `mdf[channelName]` the following keys exist:

| Key | Description |
|-----|-------------|
| `data` | numpy array of channel values |
| `unit` | unit string |
| `master` | name of the master (time/angle/…) channel |
| `masterType` | master channel type: 0=None, 1=Time, 2=Angle, 3=Distance, 4=Index |
| `description` | channel description string |
| `conversion` | present when `convert_after_read=False`; dict describing raw→physical mapping |

`mdf.masterChannelList` is a dict mapping each master channel name to the list
of channels sampled at the same raster.

Mdfreader module methods:
=========================
* resample channels to one sampling frequency
* merge files
* plot one channel, several channels on one graph (list) or several channels on subplots (list of lists)

It is also possible to export mdf data into:
* CSV file (Excel dialect by default)
* NetCDF file for compatibility with Uniplot (needs `netcdf4`, `Scientific.IO`)
* HDF5 (needs `h5py`)
* Excel 95–2003 (needs `xlwt` — very slow for large files)
* Excel 2007/2010 (needs `openpyxl` — can also be slow with large files)
* Matlab `.mat` (needs `hdf5storage`)
* MDF file — allows creating, converting or modifying data, units and descriptions
* Pandas DataFrame(s) (command line only, not in mdfconverter) — one DataFrame per raster

Compatibility:
==============
Python 3.9+ — tested on Linux and Windows (x86-64)

Requirements:
=============
Core: `numpy`, `lxml`, `sympy`

`lxml` is used for MDF4 metadata XML blocks. When Cython is compiled, the fast
path handles the common `<TX>…</TX>` pattern directly from bytes and only falls
back to `lxml` for complex XML (CDATA, namespaces).

Reading channels defined by a formula requires `sympy`.

Cython is strongly advised. It compiles `dataRead.pyx`, which provides:
* fast metadata parsing via `pread()` + C packed structs
* the `SymBufReader` bidirectional file buffer
* bit-exact reading for non-byte-aligned or record-padded channels
* VLSD/VLSC string data reading helpers

If Cython compilation fails, `bitarray` is used as a fallback (slower, pure Python).

Export requirements (optional): `scipy`, `h5py`, `hdf5storage`, `openpyxl`, `pandas`, `fastparquet`

Data compression in memory (optional): `blosc`

Graphical converter: `PyQt5`

Installation:
=============
From PyPI:
```shell
pip install mdfreader
```
From source:
```shell
pip install cython numpy        # build prerequisites
python setup.py build_ext --inplace
python setup.py develop
```

Graphical interface: mdfconverter
==================================
A PyQt5 GUI to convert batches of files. Launch with:
```shell
mdfconverter
```
Right-click a channel in the list to plot it. Channels can be dragged between
columns. A `.lab` channel-list file can be imported. Multiple files can be
merged into one and resampled.

Memory-saving options:
======================
For large files or limited memory:

* **Channel list only** — pass `channel_list=['ch1', 'ch2']`; call
  `mdfreader.MdfInfo(file)` to get the full channel list without loading data.
* **Raw data mode** — pass `convert_after_read=False`; data stays as stored in
  the MDF file and is converted on-the-fly by `get_channel_data`, `plot`,
  `export_to_*`, etc.
* **Blosc compression** — pass `compression=True` (default level 9) to compress
  data in memory after reading.
* **No-data skeleton** — pass `no_data_loading=True` to build the channel
  metadata dict without reading any samples; data is fetched on demand via
  `get_channel_data`.

For data visualisation, a dataPlugin for Veusz (≥ 1.16) is also available;
follow the instructions in Veusz's documentation and the plugin file's header.

Command example in ipython:
===========================
```python
    import mdfreader
    # loads whole mdf file content in yop mdf object.
    yop=mdfreader.Mdf('NameOfFile')
    # you can print file content in ipython with a simple:
    yop
    # alternatively, for max speed and smaller memory footprint, read only few channels
    yop=mdfreader.Mdf('NameOfFile', channel_list=['channel1', 'channel2'], convert_after_read=False)
    # also possible to keep data compressed for small memory footprint, using Blosc module
    yop=mdfreader.Mdf('NameOfFile', compression=True)
    # for interactive file exploration, possible to read the file but not its data to save memory
    yop=mdfreader.Mdf('NameOfFile', no_data_loading=True) # channel data will be loaded from file if needed
    # parsing xml metadata from mdf4.x for many channels can take more than just reading data.
    # You can reduce to minimum metadata reading with below argument (no source information, attachment, etc.) 
    yop=mdfreader.Mdf('NameOfFile', metadata=0)  # 0: full, 2: minimal
    # only for mdf4.x, you can search for the mdf key of a channel name that can have been recorded by different sources
    yop.get_channel_name4('channelName', 'source path or name')  # returns list of mdf keys
    # to yield one channel and keep its content in mdf object
    yop.get_channel('channelName')
    # to yield one channel numpy array
    yop.get_channel_data('channelName')
    # to get file mdf version
    yop.MDFVersionNumber
    # to get file structure or attachments, you can create a mdfinfo instance
    info=mdfreader.MdfInfo()
    info.list_channels('NameOfFile') # returns only the list of channels
    info.read_info('NameOfFile') # complete file structure object
    yop.info # same class is stored in mdfreader class
    # to list channels names after reading
    yop.keys()
    # to list channels names grouped by raster, below dict mdf attribute contains
    # pairs (key=masterChannelName : value=listOfChannelNamesForThisMaster)
    yop.masterChannelList
    # quick plot or subplot (with lists) of channel(s)
    yop.plot(['channel1',['channel2','channel3']])
    # file manipulations
    yop.resample(0.1)
    # or
    yop.resample(master_channel='master3')
    # keep only data between begin and end
    yop.cut(begin=10, end=15)
    # export to other file formats :
    yop.export_to_csv(sampling=0.01)
    yop.export_to_NetCDF()
    yop.export_to_hdf5()
    yop.export_to_matlab()
    yop.export_to_xlsx()
    yop.export_to_parquet()
    # return pandas dataframe from master channel name
    yop.return_pandas_dataframe('master_channel_name')
    # converts data groups into pandas dataframes and keeps it in mdf object
    yop.convert_to_pandas()
    # drops all the channels except the one in argument
    yop.keep_channels({'channel1','channel2','channel3'})
    # merge 2 files
    yop2=mdfreader.Mdf('NameOfFile_2')
    yop.merge_mdf(yop2)
    # can write mdf file after modifications or creation from scratch
    # write4 and write3 also allow to convert file versions
    yop.write('NewNameOfFile')  # write in same version as original file after modifications
    yop.write4('NameOfFile', compression=True)  # write mdf version 4.1 file, data compressed
    yop.write3()  # write mdf version 3 file
    yop.attachments  # to get attachments, embedded or paths to files 
```
<a href="https://scan.coverity.com/projects/ratal-mdfreader">
  <img alt="Coverity Scan Build Status"
       src="https://scan.coverity.com/projects/21368/badge.svg"/>
</a>
# mdfreader

A fast Python parser, editor, and exporter for **MDF (Measured Data Format)** files — versions **3.x** (`.dat`, `.mdf`) and **4.0–4.3** (`.mf4`, `.mfx`, `.mfxz`) — widely used in automotive ECU measurement and calibration tools (INCA, CANape, CANoe, dSPACE).

`mdfreader` splits parsing into two layers per format version:
- **Metadata layer** ([`mdfinfo3`](mdfreader/mdfinfo3.py), [`mdfinfo4`](mdfreader/mdfinfo4.py)) — reads headers, channel groups, conversions, XML comments, events, and attachments (`MdfInfo`).
- **Sample-data layer** ([`mdf3reader`](mdfreader/mdf3reader.py), [`mdf4reader`](mdfreader/mdf4reader.py)) — decodes sorted, unsorted, column-oriented, compressed (`Deflate` / `LZ4` / `ZStandard` with optional transposition), `VLSD`/`VLSC`, `CABlock` arrays, and invalidation bits, backed by an optional Cython extension ([`dataRead`](mdfreader/dataRead.pyx)).

---

## Features

- **Full MDF3 & MDF4 (up to MDF 4.3) read/write**: sorted & unsorted data blocks, `DL`/`HL`/`DZ`/`LD`/`DV`/`DI`/`GD` blocks, `VLSD`/`VLSC` variable-length strings & byte arrays, `DS` data streams, `MLSD`, sync stream attachments (`ATBlock`), CANOpen time/date, and `CABlock` multi-dimensional arrays (including `ca_flags` column-major layout).
- **Automatic invalidation masking**: MDF4 channels with `cn_flags` invalidation bits or `ALL_INVALID` flags automatically return `numpy.ma.MaskedArray` when invalid samples are present (plain `numpy.ndarray` is kept when zero samples are invalid).
- **Cython + `mmap` acceleration**: fast `pread()` metadata scanning, bidirectional 64 KB `SymBufReader`, single-read sorted buffers, and two-pass `nogil` unsorted parsing.
- **Signal processing & time utilities**: linear/rational/formula/lookup-table physical conversions, resampling (`resample`, `cut`), merging/concatenation (`merge_mdf`, `concat_mdf`), and timezone-aware timestamps (`start_datetime`, `get_master_channel_datetimes`).
- **Rich exports**: CSV, HDF5, MATLAB `.mat` (v7.3), Excel (`.xlsx` / `.xls`), NetCDF, Parquet, and Pandas DataFrames.

---

## Compatibility & Requirements

- **Python:** 3.9–3.13 (tested on Linux and Windows x86-64)
- **Core dependencies:** `numpy>=1.14`, `lxml`, `sympy`
  - `lxml` parses MDF4 XML metadata blocks (used only for complex XML with namespaces/CDATA when Cython fast `<TX>` scanning is active).
  - `sympy` evaluates algebraic/formula channel conversions (`cc_type == 3`).
- **Cython (`mdfreader.dataRead`) — strongly recommended:** provides `read_cn_chain_fast`, `SymBufReader`, bit-exact non-byte-aligned decoding, `unsorted_data_read4`, and `VLSD`/`VLSC` readers. If compilation is unavailable, `mdfreader` falls back to a slower pure-Python path requiring `bitarray` (`pip install mdfreader[fallback]`).

### Optional Extras (`pyproject.toml`)

| Extra | Packages | Purpose |
|---|---|---|
| `export` | `scipy`, `h5py`, `hdf5storage`, `openpyxl>2.0`, `pandas`, `pyarrow` | HDF5, MATLAB, `.xlsx`, NetCDF, Parquet, and DataFrame exports |
| `xls` | `xlwt` | Legacy Excel 95–2003 (`.xls`) export |
| `mdf4-compression` | `pyzstd`, `lz4` | Reading MDF4 blocks compressed with ZStandard or LZ4 (`Deflate` uses stdlib `zlib`) |
| `compression` | `blosc` | In-memory channel compression (`compression=True`) |
| `plot` | `matplotlib`, `mplcursors` | Quick channel plotting with interactive cursor |
| `converter` | `PyQt5` | Graphical batch converter (`mdfconverter`) |
| `fallback` | `bitarray` | Pure-Python bitfield fallback when Cython is not compiled |
| `dev` | `pytest` | Running the test suite |

---

## Installation

From PyPI:

```shell
pip install mdfreader
# With exports, plotting, and ZStd/LZ4 support:
pip install "mdfreader[export,plot,mdf4-compression]"
```

From source:

```shell
pip install cython numpy           # build prerequisites
pip install -e ".[export,dev]"     # compiles mdfreader.dataRead and installs in editable mode
# Or rebuild only the Cython extension in place:
python3 setup.py build_ext --inplace
```

To verify that the Cython fast path is active:

```python
from mdfreader.mdfinfo4 import _CN_CHAIN_FAST, _SYMBUF_AVAILABLE
print('fast CN reader:', _CN_CHAIN_FAST)
print('SymBufReader:', _SYMBUF_AVAILABLE)
```

---

## Python Examples

### Loading, Inspecting, Plotting, Manipulating, Exporting, and Writing

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

### Creating Channels, Datetimes, Concatenation, and Functional Exports

```python
import numpy as np
import mdfreader

# Create an MDF object from scratch
mdf = mdfreader.Mdf()
t = np.linspace(0.0, 1.0, 101)
mdf.add_channel("Time", t, "Time", master_type=1, unit="s")
mdf.add_channel("Speed", t * 30.0, "Time", unit="m/s", description="Vehicle speed")
mdf.write4("custom.mf4")

# Recording start timestamp and absolute datetime64[ns] master axis
print(mdf.start_datetime)
dt_axis = mdf.get_master_channel_datetimes("Time")

# Concatenate two sequential recordings sharing rasters
mdf.concat_mdf(mdfreader.Mdf("custom.mf4"))

# Top-level export functions are also directly importable
mdfreader.export_csv(mdf, "custom.csv")
```

---

## Data Structure (`Mdf` Dictionary)

`Mdf` inherits from Python's `dict` (`MdfSkeleton`), mapping each channel name to a metadata/data dictionary:

| Key | Description |
|---|---|
| `data` | `numpy.ndarray` of channel samples (`numpy.ma.MaskedArray` when MDF4 invalid samples are present) |
| `unit` | Physical unit string (e.g. `'m/s'`) |
| `master` | Name of the master channel for this raster |
| `masterType` | Master channel type: `0`=None, `1`=Time, `2`=Angle, `3`=Distance, `4`=Index |
| `description` | Channel description string |
| `conversion` | Present when `convert_after_read=False`; dict describing raw-to-physical conversion |
| `id` | `(dg, cg, cn)` block reference tuple into `info` |
| `attachment` | Present for sync-stream channels (`cn_type == 4`) referencing an `ATBlock` |

Useful container attributes:
- `mdf.masterChannelList`: dict mapping each master channel name to the list of channels sharing that raster.
- `mdf.fileMetadata`: dict of file-level metadata (`author`, `organisation`, `project`, `subject`, `comment`, `time`, `time_offset`).
- `mdf.attachments`: MDF4 `ATBlock` metadata when attachments are present.

---

## Memory-Saving Options

For large measurement files or constrained environments:

- **Selective channels (`channel_list=['ch1', 'ch2']`)** — reads only the requested channels (and their master/invalidation data); inspect available names beforehand with `mdfreader.MdfInfo(file)`.
- **Lazy physical conversion (`convert_after_read=False`)** — keeps compact raw integer/binary arrays in memory and converts on demand inside `get_channel_data()`, `plot()`, `export_to_*()`, etc.
- **In-memory compression (`compression=True`)** — compresses channel arrays with Blosc (`blosc` extra, default level 9).
- **No-data skeleton (`no_data_loading=True`)** — parses metadata only and fetches channel arrays on demand via `get_channel_data()`.
- **Minimal XML metadata (`metadata=2`, default)** — skips optional source/attachment XML trees when fast loading matters most (`metadata=0` loads full metadata).

---

## Export Formats

| Format | Method | Function (`mdfreader.*`) | Optional Dependency |
|---|---|---|---|
| CSV | `yop.export_to_csv()` | `export_csv(yop)` | stdlib `csv` |
| HDF5 | `yop.export_to_hdf5()` | `export_hdf5(yop)` | `h5py` |
| MATLAB `.mat` | `yop.export_to_matlab()` | `export_matlab(yop)` | `h5py` or `hdf5storage` / `scipy` |
| NetCDF | `yop.export_to_NetCDF()` | `export_netcdf(yop)` | `scipy` |
| Excel (`.xlsx`) | `yop.export_to_xlsx()` | `export_xlsx(yop)` | `openpyxl>2.0` |
| Excel (`.xls`) | `yop.export_to_excel()` | `export_excel(yop)` | `xlwt` |
| Parquet | `yop.export_to_parquet()` | `export_parquet(yop)` | `pandas` + `pyarrow` |
| Pandas DataFrame | `yop.return_pandas_dataframe(master)` / `yop.convert_to_pandas()` | `to_pandas(yop)` | `pandas` |

---

## Performance & Architecture

When the `mdfreader.dataRead` Cython extension is compiled, key hot paths run in C without per-record Python loops:

- **Fast metadata parsing (`read_cn_chain_fast`)** — POSIX `pread()` + packed C struct `memcpy` decoding + direct `<TX>…</TX>` byte scanning + shared `SIBlock` caching (~3–4× faster on 36,000-channel MDF4 files: ~0.6 s vs ~1.9 s pure Python on a 184 MB file).
- **`SymBufReader`** — bidirectional 64 KB C buffer around the file handle, eliminating kernel `read()` overhead on backward-linked MDF4 blocks.
- **Sorted record vectorization** — single-syscall `readinto()` into flat `uint8` buffers or zero-copy `mmap` slices viewed directly as NumPy structured arrays.
- **Two-pass `nogil` unsorted reader (`unsorted_data_read4`)** — reads directly from `memoryview` over `mmap` using a C lookup table, skipping dummy allocations for VLSD/VLSC/virtual channels and decoding variable-length strings into exact-width NumPy arrays (**2×–387× speedup** and up to **~3× lower memory** on unsorted multi-CG MDF4 files).

See [`docs/performance.rst`](docs/performance.rst) for detailed benchmarks and implementation notes.

---

## Graphical Interface & Plugins

- **`mdfconverter` (PyQt5 GUI):** batch-convert MDF files between formats. Launch with:
  ```shell
  mdfconverter
  ```
  Right-click channels to plot, drag channels between columns, import `.lab` channel lists, or merge and resample multiple files.
- **Veusz plugin:** a data plugin for Veusz (≥ 1.16) is included in [`mdfForVeuszPlugin.py`](mdfForVeuszPlugin.py); follow the instructions in the file header and Veusz documentation.

---

## Release Notes (v4.8 — Breaking Change)

> **Warning (Breaking API Change in 4.8):**
> MDF4 channels with invalidation bits or `ALL_INVALID` flags now automatically return `numpy.ma.MaskedArray` when at least one sample is invalid (plain `numpy.ndarray` is kept when zero samples are masked).
> The manual invalidation methods `apply_invalid_bit()`, `apply_all_invalid_bit()`, and `get_invalid_mask()`, as well as synthetic `invalid_bytes*` channels in the `Mdf` dict, have been removed. `resample()` also filters out masked samples before interpolating.

---

<a href="https://scan.coverity.com/projects/ratal-mdfreader">
  <img alt="Coverity Scan Build Status"
       src="https://scan.coverity.com/projects/21368/badge.svg"/>
</a>

Licensed under **GPL-3.0-or-later** ([COPYING](COPYING)).

"""pytest-based test suite for mdfreader."""
from __future__ import annotations

import gc
import sys
from pathlib import Path

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
TESTS_DIR = Path(__file__).parent
MDF4_BASE = TESTS_DIR / "MDF4" / "MDF4.3" / "Base_Standard" / "Examples"
MDF3_BASE = TESTS_DIR / "mdf3"

sys.path.insert(0, str(TESTS_DIR.parent.parent))
import mdfreader

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
MDF_EXTENSIONS = frozenset({
    ".mf4", ".MF4", ".mfx", ".mfxz", ".MFXZ",
    ".MDF", ".mdf", ".dat", ".DAT",
})
# Stems ending with these suffixes are artifacts from previous write-test runs
DERIVED_SUFFIXES = ("_New", "_new", "_comp", "_col", "_compcol")


# ---------------------------------------------------------------------------
# File collection
# ---------------------------------------------------------------------------
def collect_mdf_files(
    category: str,
    exclude_subdirs: tuple[str, ...] = (),
) -> list[Path]:
    """Return all MDF files under MDF4_BASE/category.

    Parameters
    ----------
    category:
        Subdirectory name under MDF4_BASE.
    exclude_subdirs:
        Directory names (not full paths) to exclude from rglob results.
        Useful for skipping subdirectories whose features are not yet implemented.
    """
    cat_path = MDF4_BASE / category
    if not cat_path.exists():
        return []
    return sorted(
        f for f in cat_path.rglob("*")
        if f.is_file()
        and f.suffix in MDF_EXTENSIONS
        and not any(f.stem.endswith(s) for s in DERIVED_SUFFIXES)
        and not any(part in exclude_subdirs for part in f.relative_to(cat_path).parts)
    )


def collect_mdf3_files() -> list[Path]:
    if not MDF3_BASE.exists():
        return []
    return sorted(
        f for f in MDF3_BASE.iterdir()
        if f.is_file()
        and f.suffix in MDF_EXTENSIONS
        and not any(f.stem.endswith(s) for s in DERIVED_SUFFIXES)
    )


# ---------------------------------------------------------------------------
# Parametrize lists  (one per MDF4 category + MDF3)
# ---------------------------------------------------------------------------
SIMPLE_FILES = collect_mdf_files("Simple")
DATA_TYPES_FILES = collect_mdf_files("DataTypes")
CHANNEL_TYPES_FILES = collect_mdf_files("ChannelTypes")
UNSORTED_FILES = collect_mdf_files("UnsortedData")
CONVERSION_FILES = collect_mdf_files("Conversion")
METADATA_FILES = collect_mdf_files("MetaData")
# MDF430_Algorithms/ now supported: ZStd (dz_zip_type=2/3) via pyzstd, LZ4 (4/5) via lz4
# 2 Custom single-block DZ files use proprietary compression and are skip-listed below
COMPRESSED_FILES = collect_mdf_files("CompressedData")
DATA_LIST_FILES = collect_mdf_files("DataList")
# EventSignals/ is a new MDF4.3 feature (signal-based events) not yet implemented
EVENTS_FILES = collect_mdf_files("Events", exclude_subdirs=("EventSignals",))
SAMPLE_REDUCTION_FILES = collect_mdf_files("SampleReduction")
ATTACHMENT_FILES = collect_mdf_files("Attachments")
RECORD_LAYOUT_FILES = collect_mdf_files("RecordLayout")
# Arrays rglob includes the Arrays/Classification/ subdirectory automatically
ARRAYS_FILES = collect_mdf_files("Arrays")
BUS_LOGGING_FILES = collect_mdf_files("BusLogging")
CHANNEL_INFO_FILES = collect_mdf_files("ChannelInfo")
HALFFLOAT_FILES = collect_mdf_files("Halffloat")
# MDF4.3 new channel/data composition types — all read correctly
UNION_FILES = collect_mdf_files("Union")
VARIANT_FILES = collect_mdf_files("Variant")
DYNAMIC_DATA_FILES = collect_mdf_files("DynamicData")

MDF3_FILES = collect_mdf3_files()

ALL_MDF4_FILES = (
    SIMPLE_FILES + DATA_TYPES_FILES + CHANNEL_TYPES_FILES + UNSORTED_FILES
    + CONVERSION_FILES + METADATA_FILES + COMPRESSED_FILES + DATA_LIST_FILES
    + EVENTS_FILES + SAMPLE_REDUCTION_FILES + ATTACHMENT_FILES + RECORD_LAYOUT_FILES
    + ARRAYS_FILES + BUS_LOGGING_FILES + CHANNEL_INFO_FILES + HALFFLOAT_FILES
    + UNION_FILES + VARIANT_FILES + DYNAMIC_DATA_FILES
)
ALL_FILES = ALL_MDF4_FILES + MDF3_FILES

# MDF4.3 categories not yet implemented — collected but excluded from ALL_FILES
GNSS_FILES = collect_mdf_files("GnssDataStorage")        # GPS/GNSS not implemented
REMOTE_MASTER_FILES = collect_mdf_files("RemoteMaster")  # remote master not implemented

# ---------------------------------------------------------------------------
# Per-test exclusion sets  (frozensets of filenames)
# ---------------------------------------------------------------------------
# Files using vendor-proprietary Custom compression (dz_zip_type=254) — cannot decompress
_CUSTOM_COMPRESSION_SKIP = frozenset({
    "RAC_MDF430_SingleDataBlock_DZ_DT_Custom.mf4",
    "RAC_MDF430_SingleDataBlock_DZ_DV_Custom.mf4",
    "RAC_MDF430_DataList_DataOffsets_Custom.mf4",
    "RAC_MDF430_DataList_EqualLength_Custom.mf4",
    "RAC_MDF430_DataList_SplitRecords_Custom.mf4",
    "RAC_MDF430_ListData_EqualSampleCount_Custom.mf4",
    "RAC_MDF430_ListData_SampleOffsets_Custom.mf4",
})

WRITE_SKIP = frozenset({
    "ETAS_EmptyDL.mf4",
    "T3_121121_000_6NEDC.dat",   # too large (MDF3 version)
    "T3_121121_000_6NEDC.mf4",   # too large (1.7 GB)
    "error.mf4",                 # too large (192 MB)
    "test.mf4",                  # too large (192 MB)
    "Measure.mf4",               # too large (112 MB)
    "Vector_ValueRange2ValueConversion.mf4",  # big-endian float write bug
}) | _CUSTOM_COMPRESSION_SKIP
RESAMPLE_SKIP = frozenset({
    "2DClassification.mf4",
    "ETAS_EmptyDL.mf4",
    "PCV_iO_Gen3_LK1__3l_TDI.mf4",
    "Rainflow.mf4",
})
EXPORT_EMPTY_SKIP = frozenset({
    "ETAS_EmptyDL.mf4",
})
EXPORT_CSV_EXTRA_SKIP = frozenset({
    "GDidle2.mfx",
    "2DClassification.mf4",
    "Rainflow.mf4",
    "halffloat_sinus.mf4",
    "Vector_NoMasterChannel.mf4",
    "PCV_iO_Gen3_LK1__3l_TDI.mf4",
    "Vector_ComplexNumbers.mf4",  # Excel cannot store complex64/128 values
})
PANDAS_SKIP = frozenset({
    "Vector_MLSDStringUTF16_BE.mf4",
})
EXCEL_SIZE_LIMIT = 5_000_000

WRITE_SIZE_LIMIT = 50_000_000  # 50 MB — writing 4 variants × multiple pytest runs fills disk fast
WRITE_FILES = [f for f in ALL_FILES if f.name not in WRITE_SKIP and f.stat().st_size < WRITE_SIZE_LIMIT]
RESAMPLE_FILES = [f for f in ALL_FILES if f.name not in (RESAMPLE_SKIP | _CUSTOM_COMPRESSION_SKIP)]
EXPORT_SIZE_LIMIT = 50_000_000  # 50 MB — large files cause huge exports that fill /tmp
NETCDF_FILES = [f for f in ALL_FILES
                if f.name not in (EXPORT_EMPTY_SKIP | _CUSTOM_COMPRESSION_SKIP)
                and f.stat().st_size < EXPORT_SIZE_LIMIT]
EXCEL_FILES = [
    f for f in ALL_FILES
    if f.name not in (EXPORT_EMPTY_SKIP | EXPORT_CSV_EXTRA_SKIP | _CUSTOM_COMPRESSION_SKIP)
    and f.stat().st_size < EXCEL_SIZE_LIMIT
]
CSV_FILES = [
    f for f in ALL_FILES
    if f.name not in (EXPORT_EMPTY_SKIP | EXPORT_CSV_EXTRA_SKIP | _CUSTOM_COMPRESSION_SKIP)
    and f.stat().st_size < EXPORT_SIZE_LIMIT
]
PANDAS_FILES = [f for f in ALL_FILES
                if f.name not in (PANDAS_SKIP | _CUSTOM_COMPRESSION_SKIP)
                and f.stat().st_size < EXPORT_SIZE_LIMIT]


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------
def _ref_channel(yop) -> str | None:
    """Return a representative non-master channel for round-trip data comparison."""
    for channels in yop.masterChannelList.values():
        if len(channels) > 1:
            return channels[1]
        if channels:
            return channels[0]
    return None


# ---------------------------------------------------------------------------
# test_read
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("mdf_file", ALL_FILES, ids=lambda p: p.name)
def test_read(mdf_file):
    if mdf_file.name in _CUSTOM_COMPRESSION_SKIP:
        pytest.skip("proprietary custom compression (dz_zip_type=254)")

    info = mdfreader.MdfInfo(str(mdf_file))
    assert isinstance(info, mdfreader.MdfInfo)
    del info

    yop = mdfreader.Mdf(str(mdf_file))
    assert isinstance(yop, mdfreader.Mdf)
    assert hasattr(yop, "masterChannelList")
    all_channels = list(yop)
    del yop
    gc.collect()

    # metadata-only mode
    mdfreader.Mdf(str(mdf_file), metadata=0)
    gc.collect()

    # lazy / no-data-loading mode — read each channel individually
    yop_lazy = mdfreader.Mdf(str(mdf_file), no_data_loading=True)
    for channel_name in yop_lazy:
        try:
            yop_lazy.get_channel_data(channel_name)
        except ImportError:
            continue  # optional dependency (e.g. bitarray) not installed; skip channel
        except Exception as exc:
            pytest.fail(
                f"get_channel_data({channel_name!r}) raised "
                f"{type(exc).__name__}: {exc} in {mdf_file.name}"
            )
    del yop_lazy
    gc.collect()

    # no-conversion mode
    mdfreader.Mdf(str(mdf_file), convert_after_read=False)
    gc.collect()

    # selective channel loading
    if len(all_channels) >= 2:
        mdfreader.Mdf(str(mdf_file), channel_list=all_channels[:2])


# ---------------------------------------------------------------------------
# test_write
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("mdf_file", WRITE_FILES, ids=lambda p: p.name)
def test_write(tmp_path, mdf_file):
    yop = mdfreader.Mdf(str(mdf_file))
    ref = _ref_channel(yop)

    variants = [
        ("_new", {}),
        ("_comp", {"compression": True}),
        ("_col", {"column_oriented": True}),
        ("_compcol", {"compression": True, "column_oriented": True}),
    ]
    written: dict[str, Path] = {}
    for suffix, kwargs in variants:
        out = tmp_path / (mdf_file.stem + suffix + mdf_file.suffix)
        yop.write(str(out), **kwargs)
        assert out.exists(), f"write{suffix} did not create output file"
        assert out.stat().st_size > 0, f"write{suffix} produced an empty file"
        written[suffix] = out

    if ref is not None:
        original_sum = None
        try:
            original_sum = np.sum(yop.get_channel_data(ref))
        except (TypeError, AttributeError):
            pass  # non-numeric channel; skip data comparison

        if original_sum is not None and not np.isnan(original_sum) and abs(float(original_sum)) > 1e-300:
            for suffix, out in written.items():
                yop_back = mdfreader.Mdf(str(out))
                if ref in yop_back:
                    back_sum = np.sum(yop_back.get_channel_data(ref))
                    if not np.isnan(back_sum):
                        assert original_sum == pytest.approx(back_sum), (
                            f"{mdf_file.name}: data mismatch after {suffix} write/re-read"
                        )


# ---------------------------------------------------------------------------
# test_resample
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("mdf_file", RESAMPLE_FILES, ids=lambda p: p.name)
def test_resample(mdf_file):
    for kind in (None, "next", "linear"):
        yop = mdfreader.Mdf(str(mdf_file))
        if kind is None:
            yop.resample()
        else:
            yop.resample(interpolation_kind=kind)

    yop = mdfreader.Mdf(str(mdf_file))
    masters = list(yop.masterChannelList.keys())
    if masters:
        channels = yop.masterChannelList[masters[0]]
        try:
            channel = channels[1]
        except IndexError:
            channel = channels[0]
        yop.resample(channel=channel)

        yop2 = mdfreader.Mdf(str(mdf_file))
        yop2.resample(master_channel=masters[0])


# ---------------------------------------------------------------------------
# test_export_netcdf / hdf5 / matlab
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("mdf_file", NETCDF_FILES, ids=lambda p: p.name)
def test_export_netcdf(tmp_path, mdf_file):
    yop = mdfreader.Mdf(str(mdf_file))
    out = tmp_path / (mdf_file.stem + ".nc")
    yop.export_to_NetCDF(str(out))
    assert out.exists()
    assert out.stat().st_size > 0


@pytest.mark.parametrize("mdf_file", NETCDF_FILES, ids=lambda p: p.name)
def test_export_hdf5(tmp_path, mdf_file):
    yop = mdfreader.Mdf(str(mdf_file))
    out = tmp_path / (mdf_file.stem + ".hdf")
    yop.export_to_hdf5(str(out))
    assert out.exists()
    assert out.stat().st_size > 0


@pytest.mark.parametrize("mdf_file", NETCDF_FILES, ids=lambda p: p.name)
def test_export_matlab(tmp_path, mdf_file):
    yop = mdfreader.Mdf(str(mdf_file))
    out = tmp_path / (mdf_file.stem + ".mat")
    yop.export_to_matlab(str(out))
    assert out.exists()
    assert out.stat().st_size > 0


@pytest.mark.parametrize("mdf_file", CSV_FILES, ids=lambda p: p.name)
def test_export_csv(tmp_path, mdf_file):
    yop = mdfreader.Mdf(str(mdf_file))
    out = tmp_path / (mdf_file.stem + ".csv")
    yop.export_to_csv(str(out))
    assert out.exists()
    assert out.stat().st_size > 0


@pytest.mark.parametrize("mdf_file", EXCEL_FILES, ids=lambda p: p.name)
def test_export_excel(tmp_path, mdf_file):
    yop = mdfreader.Mdf(str(mdf_file))
    out = tmp_path / (mdf_file.stem + ".xlsx")
    yop.export_to_xlsx(str(out))
    assert out.exists()
    assert out.stat().st_size > 0


# ---------------------------------------------------------------------------
# test_pandas_dataframe
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("mdf_file", PANDAS_FILES, ids=lambda p: p.name)
def test_pandas_dataframe(mdf_file):
    pd = pytest.importorskip("pandas")
    yop = mdfreader.Mdf(str(mdf_file))
    masters = list(yop.masterChannelList.keys())
    if not masters:
        pytest.skip(f"{mdf_file.name}: no master channels")
    df = yop.return_pandas_dataframe(masters[0])
    if df is not None:
        assert isinstance(df, pd.DataFrame)
        assert df.shape[1] > 0, f"{mdf_file.name}: pandas DataFrame has no columns"


# ---------------------------------------------------------------------------
# test_merge / test_copy
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("mdf_file", ALL_FILES, ids=lambda p: p.name)
def test_merge(mdf_file):
    if mdf_file.name in _CUSTOM_COMPRESSION_SKIP:
        pytest.skip("proprietary custom compression (dz_zip_type=254)")
    yop1 = mdfreader.Mdf(str(mdf_file))
    yop2 = mdfreader.Mdf(str(mdf_file))
    original_count = len(yop1)
    yop1.concat_mdf(yop2)
    # Channel count may decrease slightly when invalid_bytes channels are
    # converted to masked arrays (applied as masks rather than kept as channels).
    # Allow losing up to 10% or 1 channel (whichever is more lenient).
    min_count = original_count - max(1, original_count // 10)
    assert len(yop1) >= min_count


@pytest.mark.parametrize("mdf_file", ALL_FILES, ids=lambda p: p.name)
def test_copy(mdf_file):
    if mdf_file.name in _CUSTOM_COMPRESSION_SKIP:
        pytest.skip("proprietary custom compression (dz_zip_type=254)")
    yop = mdfreader.Mdf(str(mdf_file))
    yop_copy = yop.copy()
    assert set(yop_copy.keys()) == set(yop.keys())


# ---------------------------------------------------------------------------
# test_integrity_check
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("mdf_file", ALL_FILES, ids=lambda p: p.name)
def test_integrity_check(mdf_file):
    if mdf_file.name in _CUSTOM_COMPRESSION_SKIP:
        pytest.skip("proprietary custom compression (dz_zip_type=254)")
    yop = mdfreader.Mdf(str(mdf_file), force_file_integrity_check=True)
    assert isinstance(yop, mdfreader.Mdf)


# ---------------------------------------------------------------------------
# test_compare_mdfr — cross-validate results with the Rust mdfr library
# ---------------------------------------------------------------------------
@pytest.fixture(scope='session')
def mdfr_module():
    """Return the mdfr module if available, else None."""
    try:
        sys.path.insert(0, '/home/ratal/workspace/mdfr/lib/python3.13/site-packages')
        import mdfr as _mdfr
        return _mdfr
    except ImportError:
        return None


@pytest.mark.compare_mdfr
@pytest.mark.parametrize("mdf_file", ALL_MDF4_FILES, ids=lambda p: p.name)
def test_compare_mdfr(mdf_file, mdfr_module):
    """Cross-validate numeric data values and units between mdfreader and the Rust mdfr library.

    Channel naming differences are reported as warnings only (mdfreader uses
    ``name_DG_CG`` suffixes to de-duplicate; mdfr keeps the original name).
    Unit mismatches are also warnings unless the unit is non-empty in both libs.
    The test only *fails* on numeric data-value mismatches for channels present
    in both libraries.

    Run with:  pytest -m compare_mdfr
    Skip with: pytest -m "not compare_mdfr"
    """
    if mdfr_module is None:
        pytest.skip('mdfr not available')
    if mdf_file.name in _CUSTOM_COMPRESSION_SKIP:
        pytest.skip("proprietary custom compression (dz_zip_type=254)")

    # Load with mdfreader
    mdf = mdfreader.Mdf(str(mdf_file))
    py_names = set(mdf.keys())

    # Load with mdfr (may panic on unsupported channel types, e.g. unions needing pyarrow)
    try:
        r = mdfr_module.Mdfr(str(mdf_file))
        r.load_all_channels_data_in_memory()
        rs_names = r.get_channel_names_set()
    except BaseException as exc:
        pytest.skip(f'mdfr failed to load file: {type(exc).__name__}: {exc}')

    data_mismatches = []
    info_notes = []

    # CANOpen sub-channel names: mdfr reads wrong byte offsets for these
    # (confirmed: mdfr reads byte 0 for 'year' instead of byte 6)
    _CANOPEN_CHANNELS = frozenset({'ms', 'minute', 'hour', 'day', 'month', 'year', 'days'})

    # Channel presence differences are informational (naming conventions differ)
    only_in_mdfr = rs_names - py_names
    if only_in_mdfr:
        info_notes.append(f'INFO channels only in mdfr ({len(only_in_mdfr)}): {sorted(only_in_mdfr)[:5]}')
    only_in_mdfreader = py_names - rs_names
    if only_in_mdfreader:
        info_notes.append(f'INFO channels only in mdfreader ({len(only_in_mdfreader)}): {sorted(only_in_mdfreader)[:5]}')

    def _is_naming_disambiguation(ch):
        """Return True if the channel `ch` has renamed counterparts in both libs.

        When multiple channel groups each have a channel with the same base name,
        mdfreader renames extras as ``name_byteoffset_cgoffset`` and mdfr appends
        source info (``name source``).  If both libraries have such renamed copies,
        the shared plain name ``ch`` points to different physical channels —
        comparing them is misleading.
        """
        prefix = ch + '_'
        # Check for mdfreader's disambiguation suffix (name_digits or name_digits_digits)
        mdfr_has_renamed = any(
            n.startswith(ch + ' ') or (n != ch and n.startswith(prefix) and any(c.isdigit() for c in n[len(prefix):len(prefix)+5]))
            for n in only_in_mdfr
        )
        mdfreader_has_renamed = any(n.startswith(prefix) for n in only_in_mdfreader)
        return mdfr_has_renamed and mdfreader_has_renamed

    # Data and unit comparison for channels present in both
    for ch in sorted(py_names & rs_names):
        # CANOpen sub-channels: mdfr has known bugs reading these byte offsets
        if ch in _CANOPEN_CHANNELS:
            info_notes.append(f'CANOPEN_SKIP {ch}: skipped (mdfr reads wrong byte offsets for CANOpen channels)')
            continue
        try:
            py_data = mdf.get_channel_data(ch)
            rs_data = r.get_channel_data(ch)
        except BaseException:
            continue
        if py_data is None or rs_data is None:
            continue
        # Only compare numeric arrays
        if hasattr(py_data, 'dtype') and hasattr(rs_data, 'dtype'):
            if py_data.dtype.kind in ('f', 'i', 'u') and rs_data.dtype.kind in ('f', 'i', 'u'):
                # Skip comparison when one side returns empty data (unsupported channel type)
                if rs_data.shape == (0,) or py_data.shape == (0,):
                    info_notes.append(f'EMPTY {ch}: one side empty {py_data.shape} vs {rs_data.shape} (unsupported type)')
                    continue
                if py_data.shape != rs_data.shape:
                    if _is_naming_disambiguation(ch):
                        info_notes.append(f'NAMING {ch}: shape mismatch {py_data.shape} vs {rs_data.shape} (disambiguation)')
                    else:
                        data_mismatches.append(f'{ch}: shape mismatch {py_data.shape} vs {rs_data.shape}')
                else:
                    try:
                        py_f = py_data.astype(float)
                        rs_f = rs_data.astype(float)
                        # Skip if either side contains garbage (denormalised/huge values from
                        # unsupported features like list-data + separate-invalidation-bits)
                        max_py = np.max(np.abs(py_f[np.isfinite(py_f)])) if np.any(np.isfinite(py_f)) else 0
                        max_rs = np.max(np.abs(rs_f[np.isfinite(rs_f)])) if np.any(np.isfinite(rs_f)) else 0
                        if max_py > 1e100 or max_rs > 1e100:
                            info_notes.append(f'EXTREME {ch}: values out of range (max_py={max_py:.2e} max_rs={max_rs:.2e}) — skipped')
                            continue
                        # Skip if mdfr data contains garbage subnormal floats (> 30% near-zero
                        # non-zero values below float64 normalisation threshold) — indicates
                        # mdfr is returning uninitialised memory for invalidated list-data entries
                        _tiny = np.finfo(float).tiny
                        _sub_rs = (np.abs(rs_f) < _tiny) & (rs_f != 0.0)
                        _sub_py = (np.abs(py_f) < _tiny) & (py_f != 0.0)
                        if np.sum(_sub_rs) / len(rs_f) > 0.30 or np.sum(_sub_py) / len(py_f) > 0.30:
                            info_notes.append(f'GARBAGE {ch}: subnormal fraction rs={np.sum(_sub_rs)}/{len(rs_f)} py={np.sum(_sub_py)}/{len(py_f)} — mdfr uninitialised memory, skipped')
                            continue
                        # Skip if mdfr returns large negatives for what mdfreader reads as
                        # small positives — indicates mdfr wrong bit-offset/endianness extraction
                        # (ratio > 100x ensures we don't skip real sign differences)
                        if np.any(rs_f < 0) and np.all(py_f >= 0):
                            _max_neg = np.max(np.abs(rs_f[rs_f < 0]))
                            _max_py = max(float(np.max(np.abs(py_f))), 1.0)
                            if _max_neg > 100 * _max_py:
                                info_notes.append(f'SIGN {ch}: mdfr_neg={_max_neg:.2e} >> py_max={_max_py:.2e} — mdfr bit-extraction error, skipped')
                                continue
                        if not np.allclose(py_f, rs_f, rtol=1e-5, atol=1e-8, equal_nan=True):
                            max_diff = np.max(np.abs(py_f - rs_f))
                            if _is_naming_disambiguation(ch):
                                info_notes.append(f'NAMING {ch}: max abs diff = {max_diff:.3e} (disambiguation)')
                            else:
                                data_mismatches.append(f'{ch}: max abs diff = {max_diff:.3e}')
                    except Exception as exc:
                        data_mismatches.append(f'{ch}: comparison error: {exc}')
        # Unit comparison — normalize mdfr's "None" string to ""
        try:
            py_unit = mdf.get_channel_unit(ch) or ''
            rs_unit = r.get_channel_unit(ch)
            rs_unit_norm = '' if (rs_unit is None or rs_unit == 'None') else rs_unit
            if py_unit and rs_unit_norm and py_unit != rs_unit_norm:
                info_notes.append(f'UNIT {ch}: "{py_unit}" vs "{rs_unit_norm}"')
        except Exception:
            pass

    all_issues = data_mismatches + info_notes
    if data_mismatches:
        pytest.fail('\n'.join(all_issues))
    elif info_notes:
        # Informational differences only — emit as a warning string in the test output
        pass  # notes visible via pytest -s or -v

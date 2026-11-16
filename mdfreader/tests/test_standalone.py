"""Standalone unit and integration tests for mdfreader that require no sample files."""
import numpy as np
import pytest

import mdfreader
from mdfreader.mdf3reader import _formula_conversion as formula_conv3
from mdfreader.mdf4reader import _formula_conversion as formula_conv4


@pytest.mark.parametrize("formula", [
    "sin(X)",
    "sin(x)",
    "SIN(X)",
    "cos(X)",
    "2*x + 1",
    "exp(X)",
    "1 / (X + 1)",
])
def test_formula_conversion_mdf4_and_mdf3(formula):
    vec = np.array([0.0, 1.0, 2.0])
    res4 = formula_conv4(vec, formula)
    res3 = formula_conv3(vec, {"textFormula": formula})
    assert np.allclose(res4, res3, atol=1e-5)


def test_mdf_roundtrip_write_and_read(tmp_path):
    """Create an MDF in memory from raw numpy arrays, write to disk, read back."""
    m = mdfreader.Mdf()
    t = np.linspace(0.0, 1.0, 51)
    sig = np.sin(2 * np.pi * t)
    m.add_channel("Time", t, "Time", master_type=1, unit="s", description="Time master")
    m.add_channel("SineWave", sig, "Time", unit="V", description="Output sine wave")

    out_file = str(tmp_path / "synthetic_test.mf4")
    m.write(out_file)

    # Read back
    m_read = mdfreader.Mdf(out_file)
    assert "Time" in m_read
    assert "SineWave" in m_read
    assert np.allclose(m_read.get_channel_data("Time"), t, atol=1e-5)
    assert np.allclose(m_read.get_channel_data("SineWave"), sig, atol=1e-5)
    assert m_read.get_channel_unit("Time") == "s"
    assert m_read.get_channel_unit("SineWave") == "V"


def test_mdf_resample_and_filtering(tmp_path):
    m = mdfreader.Mdf()
    t = np.linspace(0.0, 1.0, 101)
    sig1 = t * 2.0
    sig2 = t ** 2
    m.add_channel("Time", t, "Time", master_type=1, unit="s")
    m.add_channel("Linear", sig1, "Time", unit="m")
    m.add_channel("Quad", sig2, "Time", unit="m/s")

    # Resample to 0.05s
    m.resample(0.05)
    t_res = m.get_channel_data("Time")
    assert len(t_res) == 20

    # Keep only Linear
    m.keep_channels(["Linear"])
    assert "Linear" in m
    assert "Quad" not in m


def test_export_csv_and_xlsx_with_special_characters(tmp_path):
    m = mdfreader.Mdf()
    t = np.linspace(0.0, 1.0, 11)
    sig = np.array([10.0] * 11)
    # Channel name containing null bytes and control chars (as found in some legacy files)
    bad_name = "ENGSPD_\x00\x00A5E2"
    m.add_channel("Time", t, "Time", master_type=1, unit="s")
    m.add_channel(bad_name, sig, "Time", unit="rpm\x01")

    csv_path = str(tmp_path / "export.csv")
    m.export_to_csv(csv_path)
    assert (tmp_path / "export.csv").stat().st_size > 0

    xlsx_path = str(tmp_path / "export.xlsx")
    m.export_to_xlsx(xlsx_path)
    assert (tmp_path / "export.xlsx").stat().st_size > 0


def test_mdf_single_inheritance_and_mro():
    """Verify Mdf uses single inheritance from MdfSkeleton and eliminates diamond inheritance."""
    from mdfreader.mdf import MdfSkeleton
    from mdfreader.mdf3reader import Mdf3
    from mdfreader.mdf4reader import Mdf4
    from mdfreader.mdfreader import Mdf

    # Must inherit directly from MdfSkeleton
    assert issubclass(Mdf, MdfSkeleton)
    assert not issubclass(Mdf, Mdf3)
    assert not issubclass(Mdf, Mdf4)

    # Flattened MRO: no diamond inheritance
    assert [c.__name__ for c in Mdf.__mro__] == ['Mdf', 'MdfSkeleton', 'dict', 'object']

    # All format-specific methods are delegated
    for method in ['read3', 'write3', '_get_channel_data3', '_convert_all_channel3',
                   '_convert_channel3', '_convert3']:
        assert hasattr(Mdf, method), f"Missing MDF3 delegation: {method}"

    for method in ['read4', 'write4', '_get_channel_data4', '_convert_all_channel4',
                   '_convert_channel4', '_convert_channel_data4', '_write4_column',
                   '_write4_non_column', 'get_channel_name4']:
        assert hasattr(Mdf, method), f"Missing MDF4 delegation: {method}"


def test_pandas_parquet_and_datetime_delegations(tmp_path):
    """Test dataframe, parquet export, and datetime master conversions delegated to mdfreader.export."""
    m = mdfreader.Mdf()
    t = np.linspace(0.0, 1.0, 21)
    sig = 5.0 * t
    m.add_channel("Time", t, "Time", master_type=1, unit="s")
    m.add_channel("Voltage", sig, "Time", unit="V")

    # Master datetimes
    dt = m.get_master_channel_datetimes("Time")
    assert dt is not None
    assert len(dt) == len(t)
    assert str(dt.dtype).startswith("datetime64")

    # Return pandas DataFrame
    df = m.return_pandas_dataframe("Time")
    assert df is not None
    assert "Voltage" in df.columns
    assert len(df) == 21

    # Parquet export (if optional dependency pyarrow is available)
    try:
        import pyarrow  # noqa: F401
        has_pyarrow = True
    except ImportError:
        has_pyarrow = False

    pq_path = str(tmp_path / "test.parquet")
    m.export_to_parquet(pq_path)
    if has_pyarrow:
        assert (tmp_path / "test.parquet").stat().st_size > 0

    # In-place pandas conversion
    m.convert_to_pandas()
    assert "Time_group" in m
    assert m._pandasframe is True


def test_functional_export_api(tmp_path):
    """Verify top-level export functions are directly callable."""
    m = mdfreader.Mdf()
    t = np.linspace(0.0, 1.0, 11)
    m.add_channel("Time", t, "Time", master_type=1, unit="s")
    m.add_channel("Speed", t * 10, "Time", unit="km/h")

    csv_path = str(tmp_path / "functional_export.csv")
    out = mdfreader.export_csv(m, csv_path)
    assert out == csv_path
    assert (tmp_path / "functional_export.csv").stat().st_size > 0


def test_resample_with_non_interpolatable_channels():
    """Regression test for issue #223: removing string/array channels during resample_group
    must not skip subsequent channels in masterChannelList."""
    m = mdfreader.Mdf()
    t = np.linspace(0.0, 1.0, 101)
    m.add_channel("Time", t, "Time", master_type=1, unit="s")
    m.add_channel("StrChan1", np.array(["a"] * 101, dtype="U4"), "Time")
    m.add_channel("FloatAfterStr1", t * 2.0, "Time", unit="V")
    m.add_channel("Array2DChan", np.zeros((101, 3), dtype=np.float64), "Time")
    m.add_channel("FloatAfterArray", t * 3.0, "Time", unit="A")

    m.resample(sampling=0.1)

    assert "StrChan1" not in m
    assert "Array2DChan" not in m
    assert "FloatAfterStr1" in m
    assert "FloatAfterArray" in m
    channel_lengths = {len(a["data"]) for a in m.values()}
    assert len(channel_lengths) == 1
    assert next(iter(channel_lengths)) == 10


def test_mask_integer_bits_big_endian():
    """Regression test for big-endian bitfields where aligned container exceeds span."""
    from mdfreader.mdf4reader import _mask_integer_bits

    # 3-byte span (bit_offset=5, bit_count=12) viewed in 4-byte >i4 container
    # Bytes 0..2 = 0xa0 0x24 0x80 -> 0x124 = 292; byte 3 = 0xfc (trailing neighbor byte)
    arr_pos = np.frombuffer(b"\xa0\x24\x80\xfc", dtype=">i4").copy()
    res_pos = _mask_integer_bits(arr_pos, sig_dt=3, bit_count=12, bit_offset=5, in_place=True)
    assert res_pos.tolist() == [292]

    # Negative 12-bit BE value (-1 = 0xfff << 5 in bytes 0..2 = 0x01ffe0), byte 3 = 0xab
    arr_neg = np.frombuffer(b"\x01\xff\xe0\xab", dtype=">i4")  # read-only view, in_place=False
    res_neg = _mask_integer_bits(arr_neg, sig_dt=3, bit_count=12, bit_offset=5, in_place=False)
    assert res_neg.tolist() == [-1]

    # 32-bit BE signed at bit_offset=4 in 8-byte >i8 (actual=5 bytes, extra=3 bytes)
    raw = ((-12345 & 0xFFFFFFFF) << 28).to_bytes(8, "big")
    arr_64 = np.frombuffer(raw, dtype=">i8")
    res_64 = _mask_integer_bits(arr_64, sig_dt=3, bit_count=32, bit_offset=4, in_place=False)
    assert res_64.tolist() == [-12345]


def test_apply_ca_inverse_layout():
    """Verify ca_flags bit 6 (inverse / column-major layout) reorders per-sample dimensions."""
    from mdfreader.mdf4reader import _apply_ca_inverse_layout

    # Two samples of (3, 4) stored column-major
    expected = np.array(
        [[[0.1 * i + 0.01 * j + 0.001 * s for j in range(4)] for i in range(3)] for s in range(2)]
    )
    raw_bytes = b"".join(expected[s].tobytes(order="F") for s in range(2))
    raw_c = np.frombuffer(raw_bytes, dtype="<f8").reshape((2, 3, 4))
    reordered = _apply_ca_inverse_layout(raw_c, {"CABlock": {"ca_dim_size": (3, 4), "ca_flags": 64}})
    assert reordered.flags["C_CONTIGUOUS"]
    np.testing.assert_allclose(reordered, expected)

    # Without bit 6, data stays unchanged
    unchanged = _apply_ca_inverse_layout(raw_c, {"CABlock": {"ca_dim_size": (3, 4), "ca_flags": 0}})
    assert unchanged is raw_c

    # 1D array with bit 6 is a no-op
    arr_1d = np.arange(8, dtype="<f8").reshape((2, 4))
    assert _apply_ca_inverse_layout(arr_1d, {"CABlock": {"ca_dim_size": 4, "ca_flags": 64}}) is arr_1d


def test_extract_invalid_mask_and_conversion_preservation():
    """Verify _extract_invalid_mask (bits 0/1, multi-byte, zero-invalid -> None) and CC mask preservation."""
    from mdfreader.mdf4reader import Mdf4, _extract_invalid_mask
    from mdfreader.mdf import conversionField, dataField

    # Byte 0: 0b01000001 (bits 0 and 6), Byte 1: 0b00000010 (bit 9)
    # Sample 1 byte 0 = 0b00111110 (62): bit 6 is 0 (must NOT match 6 & 0x07!)
    inv = np.array([[0x41, 0x00], [0x3E, 0x02], [0x01, 0x00], [0x40, 0x02]], dtype=np.uint8)

    mask0 = _extract_invalid_mask(inv, 0, cn_flags=0b10, n_records=4)
    mask6 = _extract_invalid_mask(inv, 6, cn_flags=0b10, n_records=4)
    mask9 = _extract_invalid_mask(inv, 9, cn_flags=0b10, n_records=4)
    assert mask0.dtype == bool and mask0.tolist() == [True, False, True, False]
    assert mask6.dtype == bool and mask6.tolist() == [True, False, False, True]
    assert mask9.dtype == bool and mask9.tolist() == [False, True, False, True]

    # Zero masked samples -> None (keep plain ndarray, do not allocate MaskedArray)
    assert _extract_invalid_mask(inv, 7, cn_flags=0b10, n_records=4) is None

    # cn_flags bit 0 (ALL_INVALID) masks all samples even without invalidation bytes
    mask_all = _extract_invalid_mask(None, 0, cn_flags=0b01, n_records=4)
    assert mask_all.tolist() == [True, True, True, True]

    # _convert_channel_data4 preserves MaskedArray mask across table/interpolation conversions
    raw_masked = np.ma.array([0.0, 9999.0, 2.0, 3.0], mask=[False, True, False, False])
    chan_dict = {
        dataField: raw_masked,
        conversionField: {"type": 4, "parameters": {"cc_val": [0.0, 0.0, 3.0, 30.0]}},
    }
    converted = Mdf4._convert_channel_data4(chan_dict, "sig", convert_tables=True)["sig"]
    assert isinstance(converted, np.ma.MaskedArray)
    assert converted.mask.tolist() == [False, True, False, False]
    np.testing.assert_allclose(converted.compressed(), [0.0, 20.0, 30.0])


def test_resample_filters_masked_samples():
    """Verify Mdf.resample skips masked placeholder values when interpolating."""
    m = mdfreader.Mdf()
    t = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    m.add_channel("Time", t, "Time", master_type=1, unit="s")
    sig = np.ma.array([0.0, 1e6, 20.0, 30.0, 40.0, 50.0, 60.0],
                      mask=[False, True, False, False, False, False, False])
    m.add_channel("Sig", sig, "Time", unit="V")
    all_bad = np.ma.array(np.ones(7), mask=np.ones(7, dtype=bool))
    m.add_channel("AllBad", all_bad, "Time", unit="V")

    m.resample(sampling=1.0)
    res_sig = m.get_channel_data("Sig")
    assert not isinstance(res_sig, np.ma.MaskedArray)
    np.testing.assert_allclose(res_sig, [0.0, 10.0, 20.0, 30.0, 40.0, 50.0])

    res_bad = m.get_channel_data("AllBad")
    assert isinstance(res_bad, np.ma.MaskedArray)
    assert res_bad.mask.all()


def test_write4_and_read_masked_array_roundtrip(tmp_path):
    """Verify MaskedArray round-trips via write4 (row/column, uncompressed/compressed) and lazy/selective reads."""
    m = mdfreader.Mdf()
    t = np.array([0.0, 0.1, 0.2, 0.3], dtype=np.float64)
    sig_masked = np.ma.array([1.0, 99.0, 3.0, 4.0], dtype=np.float64, mask=[False, True, False, False])
    sig_valid = np.array([10.0, 20.0, 30.0, 40.0], dtype=np.float64)
    m.add_channel("Time", t, "Time", master_type=1, unit="s")
    m.add_channel("SigMasked", sig_masked, "Time", unit="V")
    m.add_channel("SigValid", sig_valid, "Time", unit="A")

    for suffix, kwargs in [
        ("row.mf4", {"column_oriented": False, "compression": False}),
        ("row_comp.mf4", {"column_oriented": False, "compression": True}),
        ("col.mf4", {"column_oriented": True, "compression": False}),
        ("col_comp.mf4", {"column_oriented": True, "compression": True}),
    ]:
        out = str(tmp_path / suffix)
        m.write4(out, **kwargs)

        # Full read
        r_full = mdfreader.Mdf(out)
        assert not any(k.startswith("invalid_bytes") for k in r_full)
        d_m = r_full.get_channel_data("SigMasked")
        assert isinstance(d_m, np.ma.MaskedArray)
        assert d_m.mask.tolist() == [False, True, False, False]
        np.testing.assert_allclose(d_m.compressed(), [1.0, 3.0, 4.0])
        assert type(r_full.get_channel_data("SigValid")) is np.ndarray

        # Selective read (channel_list)
        r_sel = mdfreader.Mdf(out, channel_list=["SigMasked"])
        d_sel = r_sel.get_channel_data("SigMasked")
        assert isinstance(d_sel, np.ma.MaskedArray)
        assert d_sel.mask.tolist() == [False, True, False, False]

        # Lazy read (no_data_loading=True)
        r_lazy = mdfreader.Mdf(out, no_data_loading=True)
        d_lazy = r_lazy.get_channel_data("SigMasked")
        assert isinstance(d_lazy, np.ma.MaskedArray)
        assert d_lazy.mask.tolist() == [False, True, False, False]






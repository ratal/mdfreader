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
                   '_write4_non_column', 'apply_invalid_bit', 'get_invalid_mask',
                   'apply_all_invalid_bit', 'get_channel_name4']:
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

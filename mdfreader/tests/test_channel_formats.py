"""Tests of the MDF4 data type -> numpy / struct format helpers (no sample file needed)."""
import warnings

import pytest

from mdfreader.channel import array_format4, data_type_format4, _array_format3, _data_type_format3


@pytest.mark.parametrize("data_type, nbytes, expected", [
    (0, 1, ('<', 'u1')), (0, 2, ('<', 'u2')), (0, 3, ('<', 'u4')), (0, 8, ('<', 'u8')),
    (1, 4, ('>', 'u4')),
    (2, 1, ('<', 'i1')), (2, 8, ('<', 'i8')),
    (3, 2, ('>', 'i2')),
    (4, 2, ('<', 'f2')), (4, 4, ('<', 'f4')), (4, 8, ('<', 'f8')),
    (5, 8, ('>', 'f8')),
    (6, 10, ('', 'S10')), (8, 10, ('<', 'S10')), (9, 10, ('>', 'S10')),
    (10, 5, ('', 'V5')),
    (13, 7, ('', '')),
    (15, 4, ('<', '(2,)f2')), (15, 8, ('<', 'c8')), (16, 16, ('>', 'c16')),
])
def test_array_format4(data_type, nbytes, expected):
    assert array_format4(data_type, nbytes) == expected


@pytest.mark.parametrize("data_type, nbytes, expected", [
    (0, 1, ('<', 'B')), (0, 3, ('<', 'I')), (0, 8, ('<', 'Q')),
    (1, 2, ('>', 'H')),
    (2, 1, ('<', 'b')), (2, 8, ('<', 'q')),
    (3, 4, ('>', 'i')),
    (4, 2, ('<', 'e')), (5, 4, ('>', 'f')), (4, 8, ('<', 'd')),
    (7, 6, ('', '6s')), (8, 6, ('<', '6s')), (9, 6, ('>', '6s')),
    (15, 8, ('<', '2f')), (16, 16, ('>', '2d')),
])
def test_data_type_format4(data_type, nbytes, expected):
    assert data_type_format4(data_type, nbytes) == expected


@pytest.mark.parametrize("func", [array_format4, data_type_format4])
@pytest.mark.parametrize("data_type, nbytes", [(2, 12), (4, 3), (15, 5), (99, 4)])
def test_unsupported_formats_warn_and_return_raw(func, data_type, nbytes):
    """Unsupported combinations must warn and fall back to a raw bytes format, not crash."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        endian, fmt = func(data_type, nbytes)
    assert caught
    assert str(nbytes) in fmt


@pytest.mark.parametrize("data_type, nbytes, byte_order, expected", [
    (0, 1, 0, ('<u1', 'u1')), (0, 2, 0, ('<u2', 'u2')), (0, 4, 1, ('>u4', 'u4')), (0, 8, 0, ('<u8', 'u8')),
    (1, 1, 0, ('<i1', 'i1')), (1, 4, 1, ('>i4', 'i4')), (1, 8, 0, ('<i8', 'i8')),
    (2, 4, 0, ('<f4', 'f4')), (2, 8, 1, ('>f8', 'f8')),
    (7, 12, 0, ('S12', 'S12')),
    (8, 6, 0, ('V6', 'V6')),
    (9, 4, 0, ('>u4', 'u4')),  # big endian regardless of byte_order
    (10, 4, 0, ('>i4', 'i4')),
    (11, 4, 0, ('>f4', 'f4')),
    (13, 4, 1, ('<u4', 'u4')),  # little endian regardless of byte_order
    (14, 4, 1, ('<i4', 'i4')),
    (15, 8, 1, ('<f8', 'f8')),
])
def test_array_format3(data_type, nbytes, byte_order, expected):
    assert _array_format3(data_type, nbytes, byte_order) == expected


@pytest.mark.parametrize("data_type, nbytes, byte_order, expected", [
    (0, 1, 0, '<B'), (0, 2, 1, '>H'), (0, 4, 0, '<I'), (0, 8, 1, '>Q'),
    (1, 1, 0, '<b'), (1, 2, 1, '>h'), (1, 4, 0, '<i'), (1, 8, 1, '>q'),
    (2, 4, 0, '<f'), (2, 8, 1, '>d'),
    (7, 8, 0, '8s'), (8, 16, 0, '16s'),
    (9, 4, 0, '>I'),
    (10, 2, 0, '>h'),
    (13, 4, 1, '<I'),
    (15, 8, 1, '<d'),
])
def test_data_type_format3(data_type, nbytes, byte_order, expected):
    assert _data_type_format3(data_type, nbytes, byte_order) == expected


@pytest.mark.parametrize("data_type, nbytes", [(1, 16), (2, 2), (99, 4)])
def test_unsupported_formats3_warn_and_return_raw(data_type, nbytes):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        dt = _data_type_format3(data_type, nbytes, 0)
        arr = _array_format3(data_type, nbytes, 0)
    assert caught
    assert str(nbytes) in dt
    assert str(nbytes) in arr[0]



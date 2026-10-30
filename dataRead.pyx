# cython: language_level=3, boundscheck=False
import numpy as np
cimport numpy as np
from sys import byteorder
from libc.stdint cimport int64_t, uint8_t, uint16_t, uint32_t, uint64_t
from libc.stdio cimport printf
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from cpython.bytes cimport PyBytes_AsString
from libc.string cimport memcpy
cimport cython

# pread() is POSIX-only.  On Windows we emulate it with _lseeki64 + _read.
# C-level conditional compilation avoids the deprecated Cython IF statement.
cdef extern from * nogil:
    """
    #ifdef _WIN32
      #include <io.h>
      #include <stdio.h>
      static Py_ssize_t c_pread(int fd, void *buf, size_t count, long long offset) {
          _lseeki64(fd, offset, SEEK_SET);
          return (Py_ssize_t)_read(fd, buf, (unsigned int)count);
      }
    #else
      #include <unistd.h>
      static Py_ssize_t c_pread(int fd, void *buf, size_t count, long long offset) {
          return (Py_ssize_t)pread(fd, buf, count, (off_t)offset);
      }
    #endif
    """
    Py_ssize_t c_pread(int fd, void *buf, size_t count, int64_t offset)

# ──────────────────────────────────────────────────────────────────────────────
# SymBufReader — bidirectional-buffered file reader
#
# MDF4 metadata blocks are linked by absolute file pointers that often point
# *backward* in the file (blocks are written newest-first). Python's default
# BufferedReader fills its buffer forward, so every backward seek past the
# buffer start forces a new read syscall. This class mirrors the Rust mdfr
# SymBufReader: it keeps a fixed C-level buffer filled *centred* on the current
# position, so backward seeks within BUFSIZE/2 are served from cache.
# ──────────────────────────────────────────────────────────────────────────────

DEF SYM_BUF_SIZE = 65536   # 64 KB — same default as Rust SymBufReader

cdef class SymBufReader:
    """Bidirectional-buffered wrapper around a Python file object.

    Drop-in replacement for any ``fid`` used in mdfreader metadata parsing:
    supports ``seek(pos[, whence])``, ``read([n])``, ``tell()``, ``fileno()``.
    """
    cdef object _fid                          # underlying Python file
    cdef unsigned char _buf[SYM_BUF_SIZE]     # C-level buffer (no GC pressure)
    cdef Py_ssize_t _buf_start                # file offset of _buf[0]
    cdef Py_ssize_t _buf_len                  # valid bytes currently in _buf
    cdef Py_ssize_t _pos                      # logical file position (cursor)

    def __init__(self, fid):
        self._fid = fid
        self._buf_start = 0
        self._buf_len = 0
        self._pos = 0

    cdef int _fill(self, Py_ssize_t pos) except -1:
        """Refill _buf centred on pos, reading from the underlying file."""
        cdef Py_ssize_t start = pos - (SYM_BUF_SIZE >> 1)
        cdef bytes raw
        cdef Py_ssize_t n
        if start < 0:
            start = 0
        self._fid.seek(start)
        raw = self._fid.read(SYM_BUF_SIZE)
        n = len(raw)
        memcpy(self._buf, <const unsigned char*>raw, n)
        self._buf_start = start
        self._buf_len = n
        return 0

    def seek(self, Py_ssize_t pos, int whence=0):
        """Update logical cursor; does NOT touch the underlying file."""
        if whence == 0:
            self._pos = pos
        elif whence == 1:
            self._pos += pos
        else:   # whence == 2 (from end)
            self._fid.seek(0, 2)
            self._pos = self._fid.tell() + pos
        return self._pos

    def tell(self):
        return self._pos

    def read(self, Py_ssize_t n=-1):
        """Return up to n bytes from the current position (or all if n<0)."""
        cdef Py_ssize_t pos = self._pos
        cdef Py_ssize_t buf_end, offset, available, end
        buf_end = self._buf_start + self._buf_len
        # Fast path: serve directly from buffer
        if self._buf_len > 0 and self._buf_start <= pos < buf_end:
            offset = pos - self._buf_start
            available = self._buf_len - offset
            if n < 0 or available >= n:
                end = offset + (n if n >= 0 else available)
                data = bytes(self._buf[offset:end])
                self._pos += len(data)
                return data
        # Buffer miss: refill centred on pos, then serve
        self._fill(pos)
        offset = pos - self._buf_start
        if offset >= self._buf_len:
            return b''
        available = self._buf_len - offset
        end = offset + (n if n >= 0 else available)
        data = bytes(self._buf[offset:end])
        self._pos += len(data)
        return data

    def fileno(self):
        return self._fid.fileno()


# ─── MDF4 fast metadata reader (pread + C structs) ────────────────────────────
#
# Goal: replace the Python struct.unpack + fid.seek/read hot path for the
# CN/CC/SI/TX block chain with a single Cython function that eliminates all
# Python overhead in the metadata-reading loop.
#
# Key techniques:
#   pread()      — POSIX atomic offset read; no seek, no GIL, no Python
#                  file-object dispatch (~575 k calls eliminated for a 36 k
#                  channel file).
#   memcpy/cast  — fields extracted directly from the raw buffer into C
#                  packed-struct variables; no struct.unpack tuple allocation.
#   <TX> scan    — fast bytes.find() scan for the common MD block pattern
#                  (<CNcomment><TX>…</TX></CNcomment>), replacing lxml for
#                  >95% of files.  lxml is used only for CDATA / namespaces.
#   SI cache     — SI blocks are cached by file offset in _si_cache so that
#                  shared sources (e.g., CAN bus) are read only once.
#
# Design constraints:
#   • CC val/ref for cc_type 3/7-11 (formula, text-table, tab conversions)
#     are left for the Python CCBlock.read_cc() path — they are rare and
#     involve variable-length link arrays not worth the complexity.
#   • Composition blocks (CA/CN/DS/CL/CV/CU) are post-processed in Python
#     because they are recursive and apply only to a minority of channels.
#   • _unique_channel_name() remains in Python — it requires full Info4
#     context (allChannelList, ChannelNamesByDG) and cannot be parallelised.
# ──────────────────────────────────────────────────────────────────────────────

_SI_TYPE_MAP = {0: 'OTHER', 1: 'ECU', 2: 'BUS', 3: 'I/O', 4: 'TOOL', 5: 'USER'}
_SI_BUS_MAP  = {0: 'NONE', 1: 'OTHER', 2: 'CAN', 3: 'LIN',
                4: 'MOST', 5: 'FLEXRAY', 6: 'K_LINE', 7: 'ETHERNET', 8: 'USB'}

# C struct layout matching the on-disk MDF4 little-endian packed format.
# All structs are declared ``packed`` so that Cython/GCC adds no padding.
# Field names match the ASAM MDF4 specification (section 4.x, Table n).
#
# Reading strategy:
#   1. pread(fd, &struct, sizeof(struct), file_offset)  — reads header
#   2. data section offset = 24 + link_count * 8         — past link section
#   3. pread(fd, &data_struct, sizeof(data_struct), offset + data_offset)

# _CNFixedHdr — 88 bytes: 24-byte common header + 8 standard CN link fields.
# Extra links (attachments, default-X) follow at offset 88 when link_count > 8.
cdef packed struct _CNFixedHdr:
    char     id[4]           # b'##CN'
    uint32_t reserved        # padding (0)
    uint64_t length          # total block length in bytes
    uint64_t link_count      # number of links in link section (≥8)
    # Standard links (offsets 24–87):
    uint64_t cn_cn_next      # next CN block in channel group (0 = last)
    uint64_t cn_composition  # composition block (CA/CN/DS/CL/CV/CU) or 0
    uint64_t cn_tx_name      # TX block with channel name
    uint64_t cn_si_source    # SI block for signal source (0 if none)
    uint64_t cn_cc_conversion # CC block for value conversion (0 = identity)
    uint64_t cn_data         # signal data block (SD/DZ/HL/DL/CG) or 0
    uint64_t cn_md_unit      # TX/MD block with physical unit (0 if none)
    uint64_t cn_md_comment   # TX/MD block with channel comment (0 if none)

# _CNData — 72 bytes: data section located at offset 24 + link_count*8.
cdef packed struct _CNData:
    uint8_t  cn_type           # 0=fixed-length, 1=VLSD, 2=master,
                               # 3=virtual master, 4=sync, 5=MLSD,
                               # 6=virtual data, 7=VLSC (MDF 4.3)
    uint8_t  cn_sync_type      # 0=none, 1=time, 2=angle, 3=distance, 4=index
    uint8_t  cn_data_type      # 0-15, see ASAM MDF4 spec §4.18.2
    uint8_t  cn_bit_offset     # bit start position within first byte (0-7)
    uint32_t cn_byte_offset    # byte start position within the record
    uint32_t cn_bit_count      # total bit width of channel in the record
    uint32_t cn_flags          # bitmask; bit 17 (0x20000) = CN_F_DATA_STREAM_MODE
    uint32_t cn_invalid_bit_pos # bit position of invalidation flag in the record
    uint8_t  cn_precision      # display decimal precision (0xFF = default)
    uint8_t  cn_reserved
    uint16_t cn_attachment_count # number of AT attachment references
    double   cn_val_range_min  # raw (record) value range minimum
    double   cn_val_range_max  # raw (record) value range maximum
    double   cn_limit_min      # physical limit minimum
    double   cn_limit_max      # physical limit maximum
    double   cn_limit_ext_min  # extended limit minimum
    double   cn_limit_ext_max  # extended limit maximum

# _CCFixedHdr — 56 bytes: 24-byte header + 4 standard CC link fields.
# Additional cc_ref links follow at offset 56 when link_count > 4.
cdef packed struct _CCFixedHdr:
    char     id[4]           # b'##CC'
    uint32_t reserved
    uint64_t length
    uint64_t link_count      # number of links (≥4)
    uint64_t cc_tx_name      # TX block with conversion name (0 if none)
    uint64_t cc_md_unit      # TX/MD block with unit (0 if none)
    uint64_t cc_md_comment   # TX/MD block with comment (0 if none)
    uint64_t cc_cc_inverse   # inverse CC block (0 if none)

# _CCData — 24 bytes: data section at offset 24 + link_count*8.
# cc_val doubles and cc_ref links follow after this section.
cdef packed struct _CCData:
    uint8_t  cc_type         # 0=identity, 1=linear, 2=rational, 3=formula,
                             # 4=tab-interp, 5=tab, 6=range→value,
                             # 7=value→text, 8=range→text, 9=text→value,
                             # 10=text→text, 11=bitfield-text
    uint8_t  cc_precision    # decimal places for display (0xFF = inherit from CN)
    uint16_t cc_flags        # bitmask
    uint16_t cc_ref_count    # number of cc_ref links beyond the 4 standard links
    uint16_t cc_val_count    # number of cc_val double entries
    double   cc_phy_range_min
    double   cc_phy_range_max

# _SIBlock — 56 bytes total: 24-byte header + 3 link fields + 8 data bytes.
cdef packed struct _SIBlock:
    char     id[4]           # b'##SI'
    uint32_t reserved
    uint64_t length
    uint64_t link_count      # always 3
    uint64_t si_tx_name      # TX block with source name
    uint64_t si_tx_path      # TX block with source path
    uint64_t si_md_comment   # MD/TX block with comment
    uint8_t  si_type         # 0=OTHER, 1=ECU, 2=BUS, 3=I/O, 4=TOOL, 5=USER
    uint8_t  si_bus_type     # 0=NONE, 1=OTHER, 2=CAN, 3=LIN, 4=MOST,
                             # 5=FLEXRAY, 6=K_LINE, 7=ETHERNET, 8=USB
    uint8_t  si_flags        # bit 0 = simulated source
    char     si_reserved[5]

# _TXHdr — 24-byte common header shared by TX (plain text) and MD (XML) blocks.
# The text/XML payload immediately follows at offset 24.
cdef packed struct _TXHdr:
    char     id[4]           # b'##TX' (plain text) or b'##MD' (XML)
    uint32_t reserved
    uint64_t length          # total block length in bytes (header + payload)
    uint64_t link_count      # always 0 for TX/MD


cdef inline Py_ssize_t _pread_or_mmap(const unsigned char* mmap_buf, size_t mmap_len, int fd,
                                       void* dest, size_t count, int64_t offset) nogil:
    if mmap_buf != NULL:
        if offset >= 0 and <size_t>offset + count <= mmap_len:
            memcpy(dest, mmap_buf + offset, count)
            return count
        return 0
    elif fd >= 0:
        return c_pread(fd, dest, count, offset)
    return 0


cdef str _fast_read_tx(const unsigned char* mmap_buf, size_t mmap_len, int fd, uint64_t pointer):
    """Read TX block text via mmap or pread. Returns '' if pointer is 0 or read fails."""
    cdef _TXHdr hdr
    cdef Py_ssize_t content_len, end, nread
    cdef unsigned char* buf
    cdef str result

    if pointer == 0:
        return ''
    if mmap_buf != NULL:
        if pointer + 24 > mmap_len:
            return ''
        memcpy(&hdr, mmap_buf + pointer, 24)
        if hdr.length <= 24 or pointer + hdr.length > mmap_len:
            return ''
        content_len = <Py_ssize_t>(hdr.length - 24)
        end = content_len
        while end > 0 and mmap_buf[pointer + 24 + end - 1] == 0:
            end -= 1
        return (<bytes>(mmap_buf + pointer + 24)[:end]).decode('UTF-8', 'ignore')

    with nogil:
        nread = c_pread(fd, &hdr, 24, pointer)
    if nread < 24 or hdr.length <= 24:
        return ''
    content_len = <Py_ssize_t>(hdr.length - 24)
    if content_len <= 0:
        return ''
    buf = <unsigned char*>PyMem_Malloc(content_len + 1)
    if buf == NULL:
        return ''
    try:
        with nogil:
            nread = c_pread(fd, buf, content_len, pointer + 24)
        if nread < 1:
            return ''
        end = nread
        while end > 0 and buf[end - 1] == 0:
            end -= 1
        buf[end] = 0
        result = (<bytes>buf[:end]).decode('UTF-8', 'ignore')
    finally:
        PyMem_Free(buf)
    return result


cdef str _fast_read_tx_or_md(const unsigned char* mmap_buf, size_t mmap_len, int fd, uint64_t pointer):
    """Read a TX or MD block and return its text content."""
    cdef _TXHdr hdr
    cdef Py_ssize_t content_len, nread, tx_start, tx_end, end
    cdef unsigned char* buf
    cdef bytes payload
    cdef str result

    if pointer == 0:
        return ''
    if mmap_buf != NULL:
        if pointer + 24 > mmap_len:
            return ''
        memcpy(&hdr, mmap_buf + pointer, 24)
        if hdr.length <= 24 or pointer + hdr.length > mmap_len:
            return ''
        content_len = <Py_ssize_t>(hdr.length - 24)
        end = content_len
        while end > 0 and mmap_buf[pointer + 24 + end - 1] == 0:
            end -= 1
        if hdr.id[2] == b'T' and hdr.id[3] == b'X':
            return (<bytes>(mmap_buf + pointer + 24)[:end]).decode('UTF-8', 'ignore')
        payload = <bytes>(mmap_buf + pointer + 24)[:end]
        tx_start = payload.find(b'<TX>')
        if tx_start >= 0:
            tx_start += 4
            tx_end = payload.find(b'</TX>', tx_start)
            if tx_end >= 0:
                return payload[tx_start:tx_end].decode('UTF-8', 'ignore')
        try:
            from lxml import objectify as _obj
            xml_tree = _obj.fromstring(payload)
            try:
                return str(xml_tree.TX)
            except AttributeError:
                pass
            try:
                return str(xml_tree.CNcomment.TX)
            except AttributeError:
                pass
            try:
                return str(xml_tree.CNunit.TX)
            except AttributeError:
                pass
        except Exception:
            pass
        return ''

    with nogil:
        nread = c_pread(fd, &hdr, 24, pointer)
    if nread < 24 or hdr.length <= 24:
        return ''
    content_len = <Py_ssize_t>(hdr.length - 24)
    if content_len <= 0:
        return ''
    buf = <unsigned char*>PyMem_Malloc(content_len + 1)
    if buf == NULL:
        return ''
    try:
        with nogil:
            nread = c_pread(fd, buf, content_len, pointer + 24)
        if nread < 1:
            return ''
        end = nread
        while end > 0 and buf[end - 1] == 0:
            end -= 1
        if hdr.id[2] == b'T' and hdr.id[3] == b'X':
            buf[end] = 0
            result = (<bytes>buf[:end]).decode('UTF-8', 'ignore')
            return result
        payload = bytes(buf[:end])
        tx_start = payload.find(b'<TX>')
        if tx_start >= 0:
            tx_start += 4
            tx_end = payload.find(b'</TX>', tx_start)
            if tx_end >= 0:
                return payload[tx_start:tx_end].decode('UTF-8', 'ignore')
        try:
            from lxml import objectify as _obj
            xml_tree = _obj.fromstring(payload)
            try:
                return str(xml_tree.TX)
            except AttributeError:
                pass
            try:
                return str(xml_tree.CNcomment.TX)
            except AttributeError:
                pass
            try:
                return str(xml_tree.CNunit.TX)
            except AttributeError:
                pass
        except Exception:
            pass
        return ''
    finally:
        PyMem_Free(buf)


cdef dict _fast_read_si(const unsigned char* mmap_buf, size_t mmap_len, int fd, uint64_t pointer):
    """Read an SI block and return a dict matching the SIBlock Python format."""
    cdef _SIBlock si
    cdef Py_ssize_t nread
    cdef dict result

    if pointer == 0:
        return None
    if mmap_buf != NULL:
        if pointer + 56 > mmap_len:
            return None
        memcpy(&si, mmap_buf + pointer, 56)
    else:
        with nogil:
            nread = c_pread(fd, &si, 56, pointer)
        if nread < 56:
            return None
    result = {
        'id': b'##SI',
        'length': si.length,
        'link_count': si.link_count,
        'si_tx_name': si.si_tx_name,
        'si_tx_path': si.si_tx_path,
        'si_md_comment': si.si_md_comment,
        'si_type': _SI_TYPE_MAP.get(si.si_type, 'OTHER'),
        'si_bus_type': _SI_BUS_MAP.get(si.si_bus_type, 'NONE'),
        'si_flags': si.si_flags,
        'source_name': {'Comment': _fast_read_tx(mmap_buf, mmap_len, fd, si.si_tx_name)},
        'source_path': {'Comment': _fast_read_tx(mmap_buf, mmap_len, fd, si.si_tx_path)},
    }
    return result


def read_cn_chain_fast(object fid, uint64_t first_pointer,
                       dict si_cache, int minimal, bint channel_name_list):
    """Read the CN linked list starting at first_pointer using mmap or pread()."""
    cdef const unsigned char* mmap_buf = NULL
    cdef size_t mmap_len = 0
    cdef const unsigned char[:] mview
    cdef int fd = -1
    cdef object local_mmap = None

    try:
        mview = fid
        mmap_buf = &mview[0]
        mmap_len = <size_t>mview.shape[0]
    except (TypeError, ValueError):
        mmap_buf = NULL
        mmap_len = 0

    if mmap_buf == NULL and hasattr(fid, 'fileno'):
        try:
            fd = fid.fileno()
            import mmap as _py_mmap
            local_mmap = _py_mmap.mmap(fd, 0, access=_py_mmap.ACCESS_READ)
            mview = local_mmap
            mmap_buf = &mview[0]
            mmap_len = <size_t>mview.shape[0]
        except Exception:
            local_mmap = None
            mmap_buf = NULL
            mmap_len = 0

    cdef uint64_t pointer = first_pointer
    cdef _CNFixedHdr cn_hdr
    cdef _CNData cn_dat
    cdef _CCFixedHdr cc_hdr
    cdef _CCData cc_dat
    cdef Py_ssize_t nread, data_offset, cc_data_offset, n_extra, n_bytes
    cdef uint64_t cn_key_uint, si_ptr, cc_ptr
    cdef int64_t cn_key_neg
    cdef object cn_key
    cdef list results = []
    cdef dict cn_dict, cc_dict, si_dict
    cdef str cn_name, unit_str, desc_str
    cdef unsigned char extra_buf[512]
    cdef double* dbl_ptr
    cdef unsigned char* cc_val_buf
    cdef list cc_val_list
    cdef uint32_t i

    try:
        while pointer != 0:
            if mmap_buf != NULL:
                if pointer + 88 > mmap_len:
                    break
                memcpy(&cn_hdr, mmap_buf + pointer, 88)
            else:
                with nogil:
                    nread = c_pread(fd, &cn_hdr, 88, pointer)
                if nread < 88:
                    break

            data_offset = 24 + <Py_ssize_t>(cn_hdr.link_count) * 8
            if mmap_buf != NULL:
                if pointer + data_offset + 72 > mmap_len:
                    break
                memcpy(&cn_dat, mmap_buf + pointer + data_offset, 72)
            else:
                with nogil:
                    nread = c_pread(fd, &cn_dat, 72, pointer + data_offset)
                if nread < 72:
                    break

            if cn_dat.cn_flags & 0x20000:   # CN_F_DATA_STREAM_MODE
                cn_key_neg = -<int64_t>pointer
                cn_key = cn_key_neg
            else:
                cn_key_uint = (<uint64_t>cn_dat.cn_byte_offset) * 8 + cn_dat.cn_bit_offset
                cn_key = cn_key_uint

            cn_name = _fast_read_tx(mmap_buf, mmap_len, fd, cn_hdr.cn_tx_name) if cn_hdr.cn_tx_name else ''

            cn_dict = {
                'pointer': pointer,
                'id': b'##CN',
                'length': cn_hdr.length,
                'link_count': cn_hdr.link_count,
                'cn_cn_next': cn_hdr.cn_cn_next,
                'cn_composition': cn_hdr.cn_composition,
                'cn_tx_name': cn_hdr.cn_tx_name,
                'cn_si_source': cn_hdr.cn_si_source,
                'cn_cc_conversion': cn_hdr.cn_cc_conversion,
                'cn_data': cn_hdr.cn_data,
                'cn_md_unit': cn_hdr.cn_md_unit,
                'cn_md_comment': cn_hdr.cn_md_comment,
                'cn_type': cn_dat.cn_type,
                'cn_sync_type': cn_dat.cn_sync_type,
                'cn_data_type': cn_dat.cn_data_type,
                'cn_bit_offset': cn_dat.cn_bit_offset,
                'cn_byte_offset': cn_dat.cn_byte_offset,
                'cn_bit_count': cn_dat.cn_bit_count,
                'cn_flags': cn_dat.cn_flags,
                'cn_invalid_bit_pos': cn_dat.cn_invalid_bit_pos,
                'cn_precision': cn_dat.cn_precision,
                'cn_reserved': 0,
                'cn_attachment_count': cn_dat.cn_attachment_count,
                'cn_val_range_min': cn_dat.cn_val_range_min,
                'cn_val_range_max': cn_dat.cn_val_range_max,
                'cn_default_x': None,
                'name': cn_name,
            }

            if cn_hdr.link_count > 8:
                n_extra = (<Py_ssize_t>(cn_hdr.link_count) - 8) * 8
                if mmap_buf != NULL:
                    if pointer + 88 + n_extra <= mmap_len:
                        extra_links = []
                        for i in range(0, n_extra, 8):
                            lnk = 0
                            memcpy(&lnk, mmap_buf + pointer + 88 + i, 8)
                            extra_links.append(lnk)
                        if cn_dat.cn_attachment_count > 0:
                            cn_dict['cn_at_reference'] = extra_links[:cn_dat.cn_attachment_count]
                        if cn_hdr.link_count > 8 + cn_dat.cn_attachment_count:
                            cn_dict['cn_default_x'] = extra_links[cn_dat.cn_attachment_count:]
                elif n_extra <= 512:
                    with nogil:
                        nread = c_pread(fd, extra_buf, n_extra, pointer + 88)
                    if nread == n_extra:
                        extra_links = []
                        for i in range(0, n_extra, 8):
                            lnk = 0
                            memcpy(&lnk, extra_buf + i, 8)
                            extra_links.append(lnk)
                        if cn_dat.cn_attachment_count > 0:
                            cn_dict['cn_at_reference'] = extra_links[:cn_dat.cn_attachment_count]
                        if cn_hdr.link_count > 8 + cn_dat.cn_attachment_count:
                            cn_dict['cn_default_x'] = extra_links[cn_dat.cn_attachment_count:]

            if not channel_name_list:
                if cn_hdr.cn_md_unit:
                    unit_str = _fast_read_tx_or_md(mmap_buf, mmap_len, fd, cn_hdr.cn_md_unit)
                    if unit_str:
                        cn_dict['unit'] = unit_str
                    elif cn_dat.cn_sync_type == 1:
                        cn_dict['unit'] = 's'
                    elif cn_dat.cn_sync_type == 2:
                        cn_dict['unit'] = 'rad'
                    elif cn_dat.cn_sync_type == 3:
                        cn_dict['unit'] = 'm'
                elif cn_dat.cn_sync_type == 1:
                    cn_dict['unit'] = 's'
                elif cn_dat.cn_sync_type == 2:
                    cn_dict['unit'] = 'rad'
                elif cn_dat.cn_sync_type == 3:
                    cn_dict['unit'] = 'm'

                if cn_hdr.cn_md_comment:
                    desc_str = _fast_read_tx_or_md(mmap_buf, mmap_len, fd, cn_hdr.cn_md_comment)
                    if desc_str:
                        cn_dict['Comment'] = {'description': desc_str}

            cc_ptr = cn_hdr.cn_cc_conversion
            cc_dict = {'cc_type': 0}
            if cc_ptr != 0:
                if mmap_buf != NULL:
                    if cc_ptr + 56 <= mmap_len:
                        memcpy(&cc_hdr, mmap_buf + cc_ptr, 56)
                        nread = 56
                    else:
                        nread = 0
                else:
                    with nogil:
                        nread = c_pread(fd, &cc_hdr, 56, cc_ptr)
                if nread == 56:
                    cc_data_offset = 24 + <Py_ssize_t>(cc_hdr.link_count) * 8
                    if mmap_buf != NULL:
                        if cc_ptr + cc_data_offset + 24 <= mmap_len:
                            memcpy(&cc_dat, mmap_buf + cc_ptr + cc_data_offset, 24)
                            nread = 24
                        else:
                            nread = 0
                    else:
                        with nogil:
                            nread = c_pread(fd, &cc_dat, 24, cc_ptr + cc_data_offset)
                    if nread == 24:
                        cc_dict = {
                            'pointer': cc_ptr,
                            'id': b'##CC',
                            'length': cc_hdr.length,
                            'link_count': cc_hdr.link_count,
                            'cc_tx_name': cc_hdr.cc_tx_name,
                            'cc_md_unit': cc_hdr.cc_md_unit,
                            'cc_md_comment': cc_hdr.cc_md_comment,
                            'cc_cc_inverse': cc_hdr.cc_cc_inverse,
                            'cc_type': cc_dat.cc_type,
                            'cc_precision': cc_dat.cc_precision,
                            'cc_flags': cc_dat.cc_flags,
                            'cc_ref_count': cc_dat.cc_ref_count,
                            'cc_val_count': cc_dat.cc_val_count,
                            'cc_phy_range_min': cc_dat.cc_phy_range_min,
                            'cc_phy_range_max': cc_dat.cc_phy_range_max,
                        }
                        if cc_dat.cc_val_count > 0 and cc_dat.cc_type not in (3, 7, 8, 9, 10, 11):
                            n_bytes = cc_dat.cc_val_count * 8
                            if mmap_buf != NULL:
                                if cc_ptr + cc_data_offset + 24 + n_bytes <= mmap_len:
                                    dbl_ptr = <double*>(mmap_buf + cc_ptr + cc_data_offset + 24)
                                    cc_val_list = []
                                    for i in range(cc_dat.cc_val_count):
                                        cc_val_list.append(dbl_ptr[i])
                                    cc_dict['cc_val'] = cc_val_list
                            else:
                                cc_val_buf = <unsigned char*>PyMem_Malloc(n_bytes)
                                if cc_val_buf != NULL:
                                    try:
                                        with nogil:
                                            nread = c_pread(fd, cc_val_buf, n_bytes, cc_ptr + cc_data_offset + 24)
                                        if nread == n_bytes:
                                            dbl_ptr = <double*>cc_val_buf
                                            cc_val_list = []
                                            for i in range(cc_dat.cc_val_count):
                                                cc_val_list.append(dbl_ptr[i])
                                            cc_dict['cc_val'] = cc_val_list
                                    finally:
                                        PyMem_Free(cc_val_buf)
                        if cc_dat.cc_type in (3, 7, 8, 9, 10, 11):
                            cc_dict['_needs_completion'] = True

            if not minimal and not channel_name_list and cn_hdr.cn_si_source != 0:
                si_ptr = cn_hdr.cn_si_source
                if si_ptr not in si_cache:
                    si_cache[si_ptr] = _fast_read_si(mmap_buf, mmap_len, fd, si_ptr)
                si_dict = si_cache.get(si_ptr)
                if si_dict is not None:
                    cn_dict['SI'] = si_dict

            results.append((cn_key, cn_dict, cc_dict))
            pointer = cn_hdr.cn_cn_next

    finally:
        if local_mmap is not None:
            try:
                local_mmap.close()
            except Exception:
                pass

    return results


# ─── MDF3 fast metadata reader ────────────────────────────────────────────────

cdef packed struct _CNBlock3:
    char     id[2]            # b'CN'
    uint16_t length           # 228
    uint32_t cn_cn_next       # pointerToNextCNBlock
    uint32_t cn_cc_conversion # pointerToConversionFormula
    uint32_t cn_ce_source     # pointerToCEBlock
    uint32_t cn_cd_dependency # pointerToCDBlock
    uint32_t cn_tx_comment    # pointerToChannelCommentBlock
    uint16_t channel_type     # channelType (0=data, 1=master)
    char     signal_name[32]  # signalName
    char     signal_desc[128] # signalDescription
    uint16_t first_bit        # numberOfTheFirstBits
    uint16_t bit_count        # numberOfBits
    uint16_t signal_data_type # signalDataType
    uint16_t value_range_known# valueRangeKnown
    double   val_range_min    # valueRangeMinimum
    double   val_range_max    # valueRangeMaximum
    double   rate_sample      # rateVariableSampled
    uint32_t cn_tx_asam_name  # pointerToASAMNameBlock
    uint32_t cn_tx_si_ident   # pointerToSignalIdentifierBlock
    uint16_t byte_offset      # ByteOffset

cdef packed struct _CCBlock3:
    char     id[2]            # b'CC'
    uint16_t length           # 46
    uint16_t value_range_known# valueRangeKnown
    double   val_range_min    # valueRangeMinimum
    double   val_range_max    # valueRangeMaximum
    char     unit[20]         # physicalUnit
    uint16_t cc_type          # cc_type
    uint16_t val_pair_count   # numberOfValuePairs


cdef str _fast_read_tx3(const unsigned char* mmap_buf, size_t mmap_len, int fd, uint32_t pointer):
    cdef char block_type[2]
    cdef uint16_t block_size = 0
    cdef Py_ssize_t content_len, end, nread
    cdef unsigned char* buf
    if pointer == 0:
        return ''
    if mmap_buf != NULL:
        if pointer + 4 > mmap_len:
            return ''
        memcpy(block_type, mmap_buf + pointer, 2)
        memcpy(&block_size, mmap_buf + pointer + 2, 2)
        if block_size < 4 or pointer + block_size > mmap_len:
            return ''
        content_len = block_size - 4
        end = content_len
        while end > 0 and mmap_buf[pointer + 4 + end - 1] == 0:
            end -= 1
        return (<bytes>(mmap_buf + pointer + 4)[:end]).decode('latin1', 'replace')
    else:
        with nogil:
            nread = c_pread(fd, block_type, 2, pointer)
        if nread < 2:
            return ''
        with nogil:
            nread = c_pread(fd, &block_size, 2, pointer + 2)
        if nread < 2 or block_size < 4:
            return ''
        content_len = block_size - 4
        buf = <unsigned char*>PyMem_Malloc(content_len + 1)
        if buf == NULL:
            return ''
        try:
            with nogil:
                nread = c_pread(fd, buf, content_len, pointer + 4)
            if nread < 1:
                return ''
            end = nread
            while end > 0 and buf[end - 1] == 0:
                end -= 1
            buf[end] = 0
            return (<bytes>buf[:end]).decode('latin1', 'replace')
        finally:
            PyMem_Free(buf)


cdef dict _fast_read_ce3(const unsigned char* mmap_buf, size_t mmap_len, int fd, uint32_t pointer):
    cdef uint16_t ext_type = 0
    cdef unsigned char ce_buf[128]
    cdef Py_ssize_t nread, end
    cdef dict result = {}
    cdef str ecu_id, desc, msg, sender
    if pointer == 0:
        return None
    if mmap_buf != NULL:
        if pointer + 6 > mmap_len:
            return None
        memcpy(&ext_type, mmap_buf + pointer + 4, 2)
        if ext_type == 2:
            if pointer + 6 + 116 <= mmap_len:
                memcpy(ce_buf, mmap_buf + pointer + 6, 116)
                end = 80
                while end > 0 and ce_buf[4 + end - 1] == 0:
                    end -= 1
                desc = (<bytes>ce_buf[4:4 + end]).decode('latin1', 'replace')
                end = 32
                while end > 0 and ce_buf[84 + end - 1] == 0:
                    end -= 1
                ecu_id = (<bytes>ce_buf[84:84 + end]).decode('latin1', 'replace')
                result['tail'] = ecu_id
                result['ECU_id'] = ecu_id
                result['description'] = desc
                return result
        elif ext_type == 19:
            if pointer + 6 + 78 <= mmap_len:
                memcpy(ce_buf, mmap_buf + pointer + 6, 78)
                end = 36
                while end > 0 and ce_buf[6 + end - 1] == 0:
                    end -= 1
                msg = (<bytes>ce_buf[6:6 + end]).decode('latin1', 'replace')
                end = 36
                while end > 0 and ce_buf[42 + end - 1] == 0:
                    end -= 1
                sender = (<bytes>ce_buf[42:42 + end]).decode('latin1', 'replace')
                result['tail'] = msg
                result['message'] = msg
                result['sender'] = sender
                return result
    else:
        with nogil:
            nread = c_pread(fd, &ext_type, 2, pointer + 4)
        if nread < 2:
            return None
        if ext_type == 2:
            with nogil:
                nread = c_pread(fd, ce_buf, 116, pointer + 6)
            if nread == 116:
                end = 80
                while end > 0 and ce_buf[4 + end - 1] == 0:
                    end -= 1
                desc = (<bytes>ce_buf[4:4 + end]).decode('latin1', 'replace')
                end = 32
                while end > 0 and ce_buf[84 + end - 1] == 0:
                    end -= 1
                ecu_id = (<bytes>ce_buf[84:84 + end]).decode('latin1', 'replace')
                result['tail'] = ecu_id
                result['ECU_id'] = ecu_id
                result['description'] = desc
                return result
        elif ext_type == 19:
            with nogil:
                nread = c_pread(fd, ce_buf, 78, pointer + 6)
            if nread == 78:
                end = 36
                while end > 0 and ce_buf[6 + end - 1] == 0:
                    end -= 1
                msg = (<bytes>ce_buf[6:6 + end]).decode('latin1', 'replace')
                end = 36
                while end > 0 and ce_buf[42 + end - 1] == 0:
                    end -= 1
                sender = (<bytes>ce_buf[42:42 + end]).decode('latin1', 'replace')
                result['tail'] = msg
                result['message'] = msg
                result['sender'] = sender
                return result
    return None


def read_cn_chain_fast3(object fid, uint32_t first_pointer, dict ce_cache):
    """Read MDF3 CN block chain starting at first_pointer using mmap or pread().
    Returns list of (cn_dict, cc_dict).
    """
    cdef const unsigned char* mmap_buf = NULL
    cdef size_t mmap_len = 0
    cdef const unsigned char[:] mview
    cdef int fd = -1
    cdef object local_mmap = None

    try:
        mview = fid
        mmap_buf = &mview[0]
        mmap_len = <size_t>mview.shape[0]
    except (TypeError, ValueError):
        mmap_buf = NULL
        mmap_len = 0

    if mmap_buf == NULL and hasattr(fid, 'fileno'):
        try:
            fd = fid.fileno()
            import mmap as _py_mmap
            local_mmap = _py_mmap.mmap(fd, 0, access=_py_mmap.ACCESS_READ)
            mview = local_mmap
            mmap_buf = &mview[0]
            mmap_len = <size_t>mview.shape[0]
        except Exception:
            local_mmap = None
            mmap_buf = NULL
            mmap_len = 0

    cdef uint32_t pointer = first_pointer
    cdef _CNBlock3 cn_dat
    cdef _CCBlock3 cc_dat
    cdef Py_ssize_t nread, end
    cdef list results = []
    cdef dict cn_dict, cc_dict, ce_dict
    cdef str sig_name, sig_desc, long_name, comment_str, si_str, unit_str
    cdef double p1, p2
    cdef uint32_t cc_ptr

    try:
        while pointer != 0:
            if mmap_buf != NULL:
                if pointer + 228 > mmap_len:
                    break
                memcpy(&cn_dat, mmap_buf + pointer, 228)
            else:
                with nogil:
                    nread = c_pread(fd, &cn_dat, 228, pointer)
                if nread < 228:
                    break

            end = 32
            while end > 0 and cn_dat.signal_name[end - 1] == 0:
                end -= 1
            sig_name = (<bytes>cn_dat.signal_name[:end]).decode('latin1', 'replace')

            end = 128
            while end > 0 and cn_dat.signal_desc[end - 1] == 0:
                end -= 1
            sig_desc = (<bytes>cn_dat.signal_desc[:end]).decode('latin1', 'replace')

            if cn_dat.cn_tx_asam_name > 0:
                long_name = _fast_read_tx3(mmap_buf, mmap_len, fd, cn_dat.cn_tx_asam_name)
            else:
                long_name = ''

            if cn_dat.cn_tx_comment > 0:
                comment_str = _fast_read_tx3(mmap_buf, mmap_len, fd, cn_dat.cn_tx_comment)
            else:
                comment_str = ''

            if cn_dat.cn_tx_si_ident > 0:
                si_str = _fast_read_tx3(mmap_buf, mmap_len, fd, cn_dat.cn_tx_si_ident)
            else:
                si_str = ''

            if cn_dat.cn_ce_source > 0:
                if cn_dat.cn_ce_source in ce_cache:
                    ce_dict = ce_cache[cn_dat.cn_ce_source]
                else:
                    ce_dict = _fast_read_ce3(mmap_buf, mmap_len, fd, cn_dat.cn_ce_source)
                    if ce_dict is not None:
                        ce_cache[cn_dat.cn_ce_source] = ce_dict

            cn_dict = {
                'BlockType': b'CN',
                'BlockSize': 228,
                'pointerToNextCNBlock': cn_dat.cn_cn_next,
                'pointerToConversionFormula': cn_dat.cn_cc_conversion,
                'pointerToCEBlock': cn_dat.cn_ce_source,
                'pointerToCDBlock': cn_dat.cn_cd_dependency,
                'pointerToChannelCommentBlock': cn_dat.cn_tx_comment,
                'channelType': cn_dat.channel_type,
                'signalName': sig_name,
                'signalDescription': sig_desc,
                'numberOfTheFirstBits': cn_dat.first_bit,
                'numberOfBits': cn_dat.bit_count,
                'signalDataType': cn_dat.signal_data_type,
                'valueRangeKnown': cn_dat.value_range_known,
                'valueRangeMinimum': cn_dat.val_range_min,
                'valueRangeMaximum': cn_dat.val_range_max,
                'rateVariableSampled': cn_dat.rate_sample,
                'pointerToASAMNameBlock': cn_dat.cn_tx_asam_name,
                'pointerToSignalIdentifierBlock': cn_dat.cn_tx_si_ident,
                'ByteOffset': cn_dat.byte_offset,
                'orig_name': sig_name,
                'asam_name': long_name,
                'ChannelCommentBlock': comment_str,
                'SignalIdentifierBlock': si_str,
            }

            cc_ptr = cn_dat.cn_cc_conversion
            if cc_ptr == 0:
                cc_dict = {'cc_type': 65535, 'BlockType': b'CC', 'BlockSize': 46}
            else:
                if mmap_buf != NULL:
                    if cc_ptr + 46 <= mmap_len:
                        memcpy(&cc_dat, mmap_buf + cc_ptr, 46)
                        nread = 46
                    else:
                        nread = 0
                else:
                    with nogil:
                        nread = c_pread(fd, &cc_dat, 46, cc_ptr)
                if nread == 46:
                    end = 20
                    while end > 0 and cc_dat.unit[end - 1] == 0:
                        end -= 1
                    unit_str = (<bytes>cc_dat.unit[:end]).decode('latin1', 'replace')

                    cc_dict = {
                        'BlockType': b'CC',
                        'BlockSize': 46,
                        'valueRangeKnown': cc_dat.value_range_known,
                        'valueRangeMinimum': cc_dat.val_range_min,
                        'valueRangeMaximum': cc_dat.val_range_max,
                        'physicalUnit': unit_str,
                        'cc_type': cc_dat.cc_type,
                        'numberOfValuePairs': cc_dat.val_pair_count,
                        'conversion': {},
                    }
                    if cc_dat.cc_type == 0:
                        if mmap_buf != NULL:
                            if cc_ptr + 46 + 16 <= mmap_len:
                                memcpy(&p1, mmap_buf + cc_ptr + 46, 8)
                                memcpy(&p2, mmap_buf + cc_ptr + 46 + 8, 8)
                                cc_dict['conversion'] = {'P1': p1, 'P2': p2}
                        else:
                            with nogil:
                                nread = c_pread(fd, &p1, 8, cc_ptr + 46)
                                nread = c_pread(fd, &p2, 8, cc_ptr + 46 + 8)
                            cc_dict['conversion'] = {'P1': p1, 'P2': p2}
                    elif cc_dat.cc_type != 65535:
                        cc_dict['_needs_completion'] = True
                else:
                    cc_dict = {'cc_type': 65535, 'BlockType': b'CC', 'BlockSize': 46}

            results.append((cn_dict, cc_dict))
            pointer = cn_dat.cn_cn_next

    finally:
        if local_mmap is not None:
            try:
                local_mmap.close()
            except Exception:
                pass

    return results


@cython.boundscheck(False)
@cython.wraparound(False)
def sorted_data_read(bytes tmp, unsigned short bit_count,
        unsigned short signal_data_type, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned char bit_offset,
        unsigned long pos_byte_beg, unsigned long n_bytes, array):
    """dataRead function to read in cython a channel from a byte stream

    Parameters
    ------------
    tmp : bytes
        byte stream
    bit_count : unsigned short
        number of bit taken by the channel in the record
    signal_data_type : unsigned short
        int to describe data type
    record_format : string
        basic numpy dtype description of data type, used to create
        returned numpy ndarray
    number_of_records : unsigned long long
        number of records in byte stream
    record_byte_size : unsigned long
        number of bytes taken by one record repeated in byte stream
    bit_offset : unsigned char
        bit offset of data in C aligned bytes
    pos_byte_beg : unsigned long
        beginning byte position of channel in record
    n_bytes : unsigned long
        bytes length of channel in record
    array : boolean
        reads an array, not a vector

    Returns
    -------
    ndarray of type record_format with number_of_records records.
    Byte order is swapped if necessary to match machine byte order before bits offset and masking
    """
    cdef const char* bit_stream = PyBytes_AsString(tmp)
    if not array:
        if 'V' in record_format or 'S' in record_format or record_format is None:
            return read_byte(bit_stream, record_format, number_of_records,
                                record_byte_size, pos_byte_beg, n_bytes, bit_count, bit_offset)
        elif signal_data_type in (4, 5) and n_bytes == 4:  # float
            if (byteorder == 'little' and signal_data_type == 4) or \
                    (byteorder == 'big' and signal_data_type == 5):
                return read_float(bit_stream, record_format, number_of_records,
                                     record_byte_size, pos_byte_beg, 0)
            else: #  swap bytes
                return read_float(bit_stream, record_format, number_of_records,
                                     record_byte_size, pos_byte_beg, 1)
        elif signal_data_type in (4, 5) and n_bytes == 8:  # double
            if (byteorder == 'little' and signal_data_type == 4) or \
                    (byteorder == 'big' and signal_data_type == 5):
                return read_double(bit_stream, record_format, number_of_records,
                                      record_byte_size, pos_byte_beg, 0)
            else: #  swap bytes
                return read_double(bit_stream, record_format, number_of_records,
                                      record_byte_size, pos_byte_beg, 1)
        elif signal_data_type in (4, 5) and n_bytes == 2:  # half precision
            if (byteorder == 'little' and signal_data_type == 4) or \
                    (byteorder == 'big' and signal_data_type == 5):
                return read_half(bit_stream, record_format, number_of_records,
                                      record_byte_size, pos_byte_beg, 0)
            else: #  swap bytes
                return read_half(bit_stream, record_format, number_of_records,
                                      record_byte_size, pos_byte_beg, 1)
        elif signal_data_type in (0, 1, 13) and n_bytes == 1:  # unsigned char
            return read_unsigned_char(bit_stream, record_format, number_of_records,
                                 record_byte_size, pos_byte_beg, bit_count, bit_offset)
        elif signal_data_type in (2, 3) and n_bytes == 1:  # signed char
            return read_signed_char(bit_stream, record_format, number_of_records,
                                record_byte_size, pos_byte_beg, bit_count, bit_offset)
        elif signal_data_type in (0, 1, 13, 14) and n_bytes <= 2:  # unsigned short
            if (byteorder == 'little' and signal_data_type == 0) or \
                    (byteorder == 'big' and signal_data_type == 1):
                return read_unsigned_short(bit_stream, record_format, number_of_records,
                                      record_byte_size, pos_byte_beg, bit_count, bit_offset, 0)
            else: #  swap bytes
                return read_unsigned_short(bit_stream, record_format, number_of_records,
                                      record_byte_size, pos_byte_beg, bit_count, bit_offset, 1)
        elif signal_data_type in (2, 3) and n_bytes <= 2:  # signed short
            if (byteorder == 'little' and signal_data_type == 2) or \
                    (byteorder == 'big' and signal_data_type == 3):
                return read_signed_short(bit_stream, record_format, number_of_records,
                                     record_byte_size, pos_byte_beg, bit_count, bit_offset, 0)
            else: #  swap bytes
                return read_signed_short(bit_stream, record_format, number_of_records,
                                     record_byte_size, pos_byte_beg, bit_count, bit_offset, 1)
        elif signal_data_type in (0, 1, 14) and n_bytes <= 4:  # unsigned int
            if (byteorder == 'little' and signal_data_type == 0) or \
                    (byteorder == 'big' and signal_data_type == 1):
                return read_unsigned_int(bit_stream, record_format, number_of_records,
                                    record_byte_size, pos_byte_beg, bit_count, bit_offset, n_bytes, 0)
            else: #  swap bytes
                return read_unsigned_int(bit_stream, record_format, number_of_records,
                                    record_byte_size, pos_byte_beg, bit_count, bit_offset, n_bytes, 1)
        elif signal_data_type in (2, 3) and n_bytes <= 4:  # signed int
            if (byteorder == 'little' and signal_data_type == 2) or \
                    (byteorder == 'big' and signal_data_type == 3):
                return read_signed_int(bit_stream, record_format, number_of_records,
                                   record_byte_size, pos_byte_beg, bit_count, bit_offset, n_bytes, 0)
            else: #  swap bytes
                return read_signed_int(bit_stream, record_format, number_of_records,
                                   record_byte_size, pos_byte_beg, bit_count, bit_offset, n_bytes, 1)
        elif signal_data_type in (0, 1) and n_bytes <= 8:  # unsigned long long
            if (byteorder == 'little' and signal_data_type == 0) or \
                    (byteorder == 'big' and signal_data_type == 1):
                return read_unsigned_longlong(bit_stream, record_format, number_of_records,
                                         record_byte_size, pos_byte_beg, bit_count, bit_offset, n_bytes, 0)
            else: #  swap bytes
                return read_unsigned_longlong(bit_stream, record_format, number_of_records,
                                         record_byte_size, pos_byte_beg, bit_count, bit_offset, n_bytes, 1)
        elif signal_data_type in (2, 3) and n_bytes <= 8:  # signed long long
            if (byteorder == 'little' and signal_data_type == 2) or \
                    (byteorder == 'big' and signal_data_type == 3):
                return read_signed_longlong(bit_stream, record_format, number_of_records,
                                        record_byte_size, pos_byte_beg, bit_count, bit_offset, n_bytes, 0)
            else: #  swap bytes
                return read_signed_longlong(bit_stream, record_format, number_of_records,
                                        record_byte_size, pos_byte_beg, bit_count, bit_offset, n_bytes, 1)
        elif signal_data_type in (15, 16):  # complex
            if (byteorder == 'little' and signal_data_type == 15) or \
                    (byteorder == 'big' and signal_data_type == 16):
                swap_flag = 0
            else: #  swap bytes
                swap_flag = 1
            if n_bytes == 16:
                return read_cdouble(bit_stream, record_format, number_of_records,
                                      record_byte_size, pos_byte_beg, 0)
            elif n_bytes == 8:
                return read_cfloat(bit_stream, record_format, number_of_records,
                                      record_byte_size, pos_byte_beg, 0)
            elif n_bytes == 4:
                return read_chalf(bit_stream, record_format, number_of_records,
                                      record_byte_size, pos_byte_beg, 0)
        elif n_bytes <= 4:
            # VLSD/VLSC channels: record stores a uint pointer/size (signal_data_type 6-12)
            return read_unsigned_int(bit_stream, record_format, number_of_records,
                                     record_byte_size, pos_byte_beg, n_bytes * 8, bit_offset, n_bytes,
                                     0 if byteorder == 'little' else 1)
        elif n_bytes <= 8:
            return read_unsigned_longlong(bit_stream, record_format, number_of_records,
                                          record_byte_size, pos_byte_beg, n_bytes * 8, bit_offset, n_bytes,
                                          0 if byteorder == 'little' else 1)
        else:
            return read_byte(bit_stream, record_format, number_of_records,
                                record_byte_size, pos_byte_beg, n_bytes, bit_count, bit_offset)
    else: # array
        if (byteorder == 'little' and signal_data_type in (0, 2, 4)) or \
                    (byteorder == 'big' and signal_data_type in (1, 3, 5)):
            return read_array(bit_stream, record_format, number_of_records,
                                 record_byte_size, pos_byte_beg, n_bytes, bit_count, bit_offset, 0)
        else: #  swap bytes
            return read_array(bit_stream, record_format, number_of_records,
                                 record_byte_size, pos_byte_beg, n_bytes, bit_count, bit_offset, 1)

cdef inline read_half(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg, unsigned char swap):
    cdef uint16_t[:] buf = np.empty(number_of_records, dtype=np.uint16)
    cdef unsigned long long i
    cdef uint16_t temp_uint16 = 0  # using uint16 because float16_t is not existing
    for i in range(number_of_records):
        memcpy(&temp_uint16, &bit_stream[pos_byte_beg + record_byte_size * i], 2)
        buf[i] = temp_uint16
    if swap == 0:
        return np.asarray(buf).view(dtype=np.float16)
    else:
        return np.asarray(buf).view(dtype=np.float16).byteswap()

cdef inline read_chalf(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg, unsigned char swap):
    # complex32 = real(f16) + imag(f16): return as (n_records, 2) float16 array
    cdef np.ndarray[np.uint16_t, ndim=2] buf = np.empty((number_of_records, 2), dtype=np.uint16)
    cdef unsigned long long i
    for i in range(number_of_records):
        memcpy(&buf[i, 0], &bit_stream[pos_byte_beg + record_byte_size * i], 2)
        memcpy(&buf[i, 1], &bit_stream[pos_byte_beg + record_byte_size * i + 2], 2)
    if swap == 0:
        return buf.view(dtype=np.float16)
    else:
        return buf.view(dtype=np.float16).byteswap()

cdef inline read_float(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg, unsigned char swap):
    cdef np.ndarray[np.float32_t] buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef float temp_float = 0
    for i in range(number_of_records):
        memcpy(&temp_float, &bit_stream[pos_byte_beg + record_byte_size * i], 4)
        buf[i] = temp_float
    if swap == 0:
        return buf
    else:
        return buf.byteswap()

cdef inline read_cfloat(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg, unsigned char swap):
    cdef np.ndarray[np.complex64_t] buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef float complex temp_cfloat = 0
    for i in range(number_of_records):
        memcpy(&temp_cfloat, &bit_stream[pos_byte_beg + record_byte_size * i], 8)
        buf[i] = temp_cfloat
    if swap == 0:
        return buf
    else:
        return buf.byteswap()

cdef inline read_double(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg, unsigned char swap):
    cdef np.ndarray[np.float64_t] buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef double temp_double = 0
    for i in range(number_of_records):
        memcpy(&temp_double, &bit_stream[pos_byte_beg + record_byte_size * i], 8)
        buf[i] = temp_double
    if swap == 0:
        return buf
    else:
        return buf.byteswap()

cdef inline read_cdouble(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg, unsigned char swap):
    cdef np.ndarray[np.complex128_t] buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef double complex temp_cdouble = 0
    for i in range(number_of_records):
        memcpy(&temp_cdouble, &bit_stream[pos_byte_beg + record_byte_size * i], 16)
        buf[i] = temp_cdouble
    if swap == 0:
        return buf
    else:
        return buf.byteswap()

cdef inline read_unsigned_char(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg,
        unsigned long bit_count, unsigned char bit_offset):
    cdef np.ndarray[np.uint8_t] buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef unsigned char mask = ((1 << bit_count) - 1)
    cdef unsigned char temp1byte = 0
    if bit_count == 8:
        for i in range(number_of_records):
            memcpy(&temp1byte, &bit_stream[pos_byte_beg + record_byte_size * i], 1)
            buf[i] = temp1byte
    else:
        for i in range(number_of_records):
            memcpy(&temp1byte, &bit_stream[pos_byte_beg + record_byte_size * i], 1)
            # right shift 
            if bit_offset > 0:
                temp1byte = temp1byte >> bit_offset
            # mask left part
            temp1byte &= mask
            buf[i] = temp1byte
    return buf

cdef inline read_signed_char(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg,
        unsigned long bit_count, unsigned char bit_offset):
    cdef np.ndarray[np.int8_t] buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef unsigned char mask = ((1 << bit_count) - 1)
    cdef char temp1byte = 0
    cdef unsigned char sign_bit = 0
    cdef unsigned char sign_bit_mask = (1 << (bit_count-1))
    cdef unsigned char sign_extend = ((1 << (8 - bit_count)) - 1) << bit_count
    if bit_count == 8:
        for i in range(number_of_records):
            memcpy(&temp1byte, &bit_stream[pos_byte_beg + record_byte_size * i], 1)
            buf[i] = temp1byte
    else:
        for i in range(number_of_records):
            memcpy(&temp1byte, &bit_stream[pos_byte_beg + record_byte_size * i], 1)
            # right shift 
            if bit_offset > 0:
                temp1byte = temp1byte >> bit_offset
            # mask left part
            temp1byte &= mask
            sign_bit = temp1byte & sign_bit_mask
            if sign_bit: #  negative value, sign extend
                temp1byte |= sign_extend
            buf[i] = temp1byte
    return buf

cdef inline read_unsigned_short(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg,
        unsigned long bit_count, unsigned char bit_offset, unsigned char swap):
    cdef np.ndarray[np.uint16_t] buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef unsigned short mask = ((1 << bit_count) - 1)
    cdef unsigned short temp2byte = 0
    cdef unsigned char temp[2]
    if bit_count == 16:
        for i in range(number_of_records):
            memcpy(&temp2byte, &bit_stream[pos_byte_beg + record_byte_size * i], 2)
            buf[i] = temp2byte
        if swap == 0:
            return buf
        else:
            return buf.byteswap()
    else:
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp2byte, &bit_stream[pos_byte_beg + record_byte_size * i], 2)
                # right shift
                if bit_offset > 0:
                    temp2byte = temp2byte >> bit_offset
                # mask left part
                if bit_count < 16:
                    temp2byte &= mask
                buf[i] = temp2byte
        else:
            for i in range(number_of_records):
                memcpy(&temp, &bit_stream[pos_byte_beg + record_byte_size * i], 2)
                temp2byte = temp[0]<<8 | temp[1]  #  swap bytes
                # right shift
                if bit_offset > 0:
                    temp2byte = temp2byte >> bit_offset
                # mask left part
                if bit_count < 16:
                    temp2byte &= mask
                buf[i] = temp2byte
        return buf

cdef inline read_signed_short(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg,
        unsigned long bit_count, unsigned char bit_offset, unsigned char swap):
    cdef np.ndarray[np.int16_t] buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef unsigned short mask = ((1 << bit_count) - 1)
    cdef short temp2byte = 0
    cdef unsigned short sign_bit = 0
    cdef unsigned short sign_bit_mask = (1 << (bit_count-1))
    cdef unsigned short sign_extend = ((1 << (16 - bit_count)) - 1) << bit_count
    cdef unsigned char temp[2]
    if bit_count == 16:
        for i in range(number_of_records):
            memcpy(&temp2byte, &bit_stream[pos_byte_beg + record_byte_size * i], 2)
            buf[i] = temp2byte
        if swap == 0:
            return buf
        else:
            return buf.byteswap()
    else:
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp2byte, &bit_stream[pos_byte_beg + record_byte_size * i], 2)
                # right shift
                if bit_offset > 0:
                    temp2byte = temp2byte >> bit_offset
                # mask left part
                temp2byte &= mask
                sign_bit = temp2byte & sign_bit_mask
                if sign_bit: #  negative value, sign extend
                    temp2byte |= sign_extend
                buf[i] = temp2byte
        else:
            for i in range(number_of_records):
                memcpy(&temp, &bit_stream[pos_byte_beg + record_byte_size * i], 2)
                temp2byte = temp[0]<<8 | temp[1]  #  swap bytes
                # right shift
                if bit_offset > 0:
                    temp2byte = temp2byte >> bit_offset
                # mask left part
                temp2byte &= mask
                sign_bit = temp2byte & sign_bit_mask
                if sign_bit: #  negative value, sign extend
                    temp2byte |= sign_extend
                buf[i] = temp2byte
        return buf

cdef inline read_unsigned_int(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg,
        unsigned long bit_count, unsigned char bit_offset, unsigned long n_bytes, unsigned char swap):
    cdef np.ndarray[np.uint32_t] buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef unsigned int mask = ((1 << bit_count) - 1)
    cdef unsigned int temp4byte = 0
    cdef unsigned char temp4[4]
    cdef unsigned char temp3[3]
    if bit_count == 32:
        for i in range(number_of_records):
            memcpy(&temp4byte, &bit_stream[pos_byte_beg + record_byte_size * i], 4)
            buf[i] = temp4byte
        if swap == 0:
            return buf
        else:
            return buf.byteswap()
    elif n_bytes == 4:
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp4byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                # right shift
                if bit_offset > 0:
                    temp4byte = temp4byte >> bit_offset
                # mask left part
                if bit_count < 32:
                    temp4byte &= mask
                buf[i] = temp4byte
        else:
            for i in range(number_of_records):
                memcpy(&temp4, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                temp4byte = temp4[0]<<24 | temp4[1]<<16 | temp4[2]<<8 | temp4[3]  #  swap bytes
                # right shift
                if bit_offset > 0:
                    temp4byte = temp4byte >> bit_offset
                # mask left part
                if bit_count < 32:
                    temp4byte &= mask
                buf[i] = temp4byte
        return buf
    else:  # on 3 bytes
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp4byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                # right shift 
                if bit_offset > 0:
                    temp4byte = temp4byte >> bit_offset
                # mask left part
                if bit_count < 24:
                    temp4byte &= mask
                buf[i] = temp4byte
        else:
            for i in range(number_of_records):
                memcpy(&temp3, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                temp4byte = temp3[0]<<16 | temp3[1]<<8 | temp3[2]  #  swap bytes
                # right shift 
                if bit_offset > 0:
                    temp4byte = temp4byte >> bit_offset
                # mask left part
                if bit_count < 24:
                    temp4byte &= mask
                buf[i] = temp4byte
        return buf


cdef inline read_signed_int(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg,
        unsigned long bit_count, unsigned char bit_offset, unsigned long n_bytes, unsigned char swap):
    cdef np.ndarray[np.int32_t] buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef unsigned int mask = ((1 << bit_count) - 1)
    cdef int temp4byte = 0
    cdef unsigned int sign_bit = 0
    cdef unsigned int sign_bit_mask = (1 << (bit_count-1))
    cdef unsigned int sign_extend = ((1 << (32 - bit_count)) - 1) << bit_count
    cdef unsigned char temp4[4]
    cdef unsigned char temp3[3]
    if bit_count == 32:
        for i in range(number_of_records):
            memcpy(&temp4byte, &bit_stream[pos_byte_beg + record_byte_size * i], 4)
            buf[i] = temp4byte
        if swap == 0:
            return buf
        else:
            return buf.byteswap()
    elif n_bytes == 4:
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp4byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                # right shift
                if bit_offset > 0:
                    temp4byte = temp4byte >> bit_offset
                # mask left part
                if bit_count < 32:
                    temp4byte &= mask
                sign_bit = temp4byte & sign_bit_mask # assumes return in little endian, to be reviewed
                if sign_bit: #  negative value, sign extend
                    temp4byte |= sign_extend
                buf[i] = temp4byte
        else:
            for i in range(number_of_records):
                memcpy(&temp4, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                temp4byte = temp4[0]<<24 | temp4[1]<<16 | temp4[2]<<8 | temp4[3]  #  swap bytes
                # right shift
                if bit_offset > 0:
                    temp4byte = temp4byte >> bit_offset
                # mask left part
                if bit_count < 32:
                    temp4byte &= mask
                sign_bit = temp4byte & sign_bit_mask # assumes return in little endian, to be reviewed
                if sign_bit: #  negative value, sign extend
                    temp4byte |= sign_extend
                buf[i] = temp4byte
        return buf
    else:  # on 3 bytes
        if swap == 0:
            for i in range(number_of_records):
                temp4byte = 0  # must zero high byte: memcpy only writes n_bytes=3, and sign_extend may have set it to 0xFF
                memcpy(&temp4byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                # right shift
                if bit_offset > 0:
                    temp4byte = temp4byte >> bit_offset
                # mask left part
                if bit_count < 24:
                    temp4byte &= mask
                sign_bit = temp4byte & sign_bit_mask # assumes return in little endian, to be reviewed
                if sign_bit: #  negative value, sign extend
                    temp4byte |= sign_extend
                buf[i] = temp4byte
        else:
            for i in range(number_of_records):
                memcpy(&temp3, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                temp4byte = temp3[0]<<16 | temp3[1]<<8 | temp3[2]  #  swap bytes
                # right shift 
                if bit_offset > 0:
                    temp4byte = temp4byte >> bit_offset
                # mask left part
                if bit_count < 24:
                    temp4byte &= mask
                sign_bit = temp4byte & sign_bit_mask # assumes return in little endian, to be reviewed
                if sign_bit: #  negative value, sign extend
                    temp4byte |= sign_extend
                buf[i] = temp4byte
        return buf

cdef inline read_unsigned_longlong(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg,
        unsigned long bit_count, unsigned char bit_offset, unsigned long n_bytes, unsigned char swap):
    cdef np.ndarray[np.uint64_t] buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long long mask = ((1 << bit_count) - 1)
    cdef unsigned long long temp8byte = 0
    cdef unsigned char temp8[8]
    cdef unsigned char temp7[7]
    cdef unsigned char temp6[6]
    cdef unsigned char temp5[5]
    if bit_count == 64:
        for i in range(number_of_records):
            memcpy(&temp8byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
            buf[i] = temp8byte
        if swap == 0:
            return buf
        else:
            return buf.byteswap()
    elif n_bytes == 8:
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp8byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 64:
                    temp8byte &= mask
                buf[i] = temp8byte
        else:
            for i in range(number_of_records):
                memcpy(&temp8, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                temp8byte = <uint64_t>temp8[0]<<56 | <uint64_t>temp8[1]<<48 | \
                            <uint64_t>temp8[2]<<40 | <uint64_t>temp8[3]<<32 | \
                            temp8[4]<<24 | temp8[5]<<16 | temp8[6]<<8 | temp8[7] #  swap bytes
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 64:
                    temp8byte &= mask
                buf[i] = temp8byte
    elif n_bytes == 7:
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp8byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 56:
                    temp8byte &= mask
                buf[i] = temp8byte
        else:
            for i in range(number_of_records):
                memcpy(&temp7, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                temp8byte = <uint64_t>temp7[0]<<48 | <uint64_t>temp7[1]<<40 | <uint64_t>temp7[2]<<32 | \
                            temp7[3]<<24 | temp7[4]<<16 | temp7[5]<<8 | temp7[6] #  swap bytes
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 56:
                    temp8byte &= mask
                buf[i] = temp8byte
    elif n_bytes == 6:
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp8byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 48:
                    temp8byte &= mask
                buf[i] = temp8byte
        else:
            for i in range(number_of_records):
                memcpy(&temp6, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                temp8byte = <uint64_t>temp6[0]<<40 | <uint64_t>temp6[1]<<32 | temp6[2]<<24 | \
                            temp6[3]<<16 | temp6[4]<<8 | temp6[5] #  swap bytes
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 48:
                    temp8byte &= mask
                buf[i] = temp8byte
    elif n_bytes == 5:
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp8byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 32:
                    temp8byte &= mask
                buf[i] = temp8byte
        else:
            for i in range(number_of_records):
                memcpy(&temp5, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                temp8byte = <uint64_t>temp5[0]<<32 | temp5[1]<<24 | \
                            temp5[2]<<16 | temp5[3]<<8 | temp5[4] #  swap bytes
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 32:
                    temp8byte &= mask
                buf[i] = temp8byte
    return buf

cdef inline read_signed_longlong(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg,
        unsigned long bit_count, unsigned char bit_offset, unsigned long n_bytes, unsigned char swap):
    cdef np.ndarray[np.int64_t] buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long long mask = ((1 << bit_count) - 1)
    cdef long long temp8byte = 0
    cdef unsigned long sign_bit = 0
    cdef unsigned long long sign_bit_mask = (1 << (bit_count-1))
    cdef unsigned long long sign_extend = ((1 << (64 - bit_count)) - 1) << bit_count
    cdef unsigned char temp8[8]
    cdef unsigned char temp7[7]
    cdef unsigned char temp6[6]
    cdef unsigned char temp5[5]
    if bit_count == 64:
        for i in range(number_of_records):
            memcpy(&temp8byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
            buf[i] = temp8byte
        if swap == 0:
            return buf
        else:
            return buf.byteswap()
    elif n_bytes == 8:
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp8byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 64:
                    temp8byte &= mask
                sign_bit = temp8byte & sign_bit_mask
                if sign_bit: #  negative value, sign extend
                    temp8byte |= sign_extend
                buf[i] = temp8byte
        else:
            for i in range(number_of_records):
                memcpy(&temp8, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                temp8byte = <uint64_t>temp8[0]<<56 | <uint64_t>temp8[1]<<48 | \
                            <uint64_t>temp8[2]<<40 | <uint64_t>temp8[3]<<32 | \
                            temp8[4]<<24 | temp8[5]<<16 | temp8[6]<<8 | temp8[7] #  swap bytes
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 64:
                    temp8byte &= mask
                sign_bit = temp8byte & sign_bit_mask
                if sign_bit: #  negative value, sign extend
                    temp8byte |= sign_extend
                buf[i] = temp8byte
    elif n_bytes == 7:
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp8byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 56:
                    temp8byte &= mask
                sign_bit = temp8byte & sign_bit_mask
                if sign_bit: #  negative value, sign extend
                    temp8byte |= sign_extend
                buf[i] = temp8byte
        else:
            for i in range(number_of_records):
                memcpy(&temp7, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                temp8byte = <uint64_t>temp7[0]<<48 | <uint64_t>temp7[1]<<40 | <uint64_t>temp7[2]<<32 | \
                            temp7[3]<<24 | temp7[4]<<16 | temp7[5]<<8 | temp7[6] #  swap bytes
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 56:
                    temp8byte &= mask
                sign_bit = temp8byte & sign_bit_mask
                if sign_bit: #  negative value, sign extend
                    temp8byte |= sign_extend
                buf[i] = temp8byte
    elif n_bytes == 6:
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp8byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 48:
                    temp8byte &= mask
                sign_bit = temp8byte & sign_bit_mask
                if sign_bit: #  negative value, sign extend
                    temp8byte |= sign_extend
                buf[i] = temp8byte
        else:
            for i in range(number_of_records):
                memcpy(&temp6, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                temp8byte = <uint64_t>temp6[0]<<40 | <uint64_t>temp6[1]<<32 | temp6[2]<<24 | \
                            temp6[3]<<16 | temp6[4]<<8 | temp6[5] #  swap bytes
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 48:
                    temp8byte &= mask
                sign_bit = temp8byte & sign_bit_mask
                if sign_bit: #  negative value, sign extend
                    temp8byte |= sign_extend
                buf[i] = temp8byte
    elif n_bytes == 5:
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp8byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 40:
                    temp8byte &= mask
                sign_bit = temp8byte & sign_bit_mask
                if sign_bit: #  negative value, sign extend
                    temp8byte |= sign_extend
                buf[i] = temp8byte
        else:
            for i in range(number_of_records):
                memcpy(&temp5, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                temp8byte = <uint64_t>temp5[0]<<32 | temp5[1]<<24 | \
                            temp5[2]<<16 | temp5[3]<<8 | temp5[4] #  swap bytes
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 40:
                    temp8byte &= mask
                sign_bit = temp8byte & sign_bit_mask
                if sign_bit: #  negative value, sign extend
                    temp8byte |= sign_extend
                buf[i] = temp8byte
    return buf

cdef inline read_byte(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg, unsigned long n_bytes,
        unsigned long bit_count, unsigned char bit_offset):
    cdef np.ndarray buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long pos_byte_end = pos_byte_beg + n_bytes
    for i in range(number_of_records):
            buf[i] = bytes(bit_stream[pos_byte_beg + record_byte_size * i:\
                    pos_byte_end + record_byte_size * i])
    return buf

cdef inline read_array(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg, unsigned long n_bytes,
        unsigned long bit_count, unsigned char bit_offset, unsigned char swap):
    cdef np.ndarray buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long pos_byte_end = pos_byte_beg + n_bytes
    for i in range(number_of_records):
        buf[i] = np.frombuffer(bit_stream[pos_byte_beg + record_byte_size * i:\
            pos_byte_end + record_byte_size * i], dtype=record_format)
    if swap == 0:
        return buf
    else:
        return buf.byteswap()

def unsorted_data_read4(record, info, bytes tmp,
                        const unsigned short record_id_size,
                        const unsigned long long data_block_length):
    """ reads only the channels using offset functions, channel by channel within unsorted data

    Parameters
    ------------
    record : class
        record class
    info: class
        info class
    tmp : bytes
        byte stream
    record_id_size : unsigned short
        record id length
    data_block_length : unsigned long long
        length of data block minus header

    Returns
    --------
    buf : array
        data array

    """
    cdef const char* bit_stream = PyBytes_AsString(tmp)
    cdef unsigned long long position = 0
    cdef unsigned char record_id_char = 0
    cdef unsigned short record_id_short = 0
    cdef unsigned long record_id_long = 0
    cdef unsigned long long record_id_long_long = 0
    # initialise data structure
    # key is channel name
    cdef dict buf = {}
    cdef dict VLSD = {}
    cdef dict pos_byte_beg = {}
    cdef dict pos_byte_end = {}
    cdef dict c_format_structure = {}
    cdef dict byte_length = {}
    cdef dict numpy_format = {}
    # key is record id
    cdef dict index = {}
    cdef dict CGrecordLength = {}
    cdef dict VLSD_flag = {}
    cdef dict VLSD_CG_name = {}
    cdef dict VLSD_CG_signal_data_type = {}
    cdef dict channel_name_set = {}
    for record_id in record:
        if record[record_id]['record'].Flags & 0b100001:  # VLSD (bit 0) or VLSC compact (bit 5)
            VLSD_flag[record_id] = True
            VLSD[record[record_id]['record'].VLSD_CG[record_id]['channelName']] = []
            VLSD_CG_name[record_id] = record[record_id]['record'].VLSD_CG[record_id]['channelName']
            VLSD_CG_signal_data_type[record_id] = record[record_id]['record'].VLSD_CG[record_id]['channel'].signal_data_type(info)
        else:
            VLSD_flag[record_id] = False
            for Channel in record[record_id]['record'].values():
                #if not Channel.VLSD_CG_Flag:
                buf[Channel.name] = np.empty((record[record_id]['record'].numberOfRecords,),
                                             dtype='V{}'.format(Channel.nBytes_aligned), order ='C')
                numpy_format[Channel.name] = Channel.data_format(info)
                pos_byte_beg[Channel.name] = record_id_size + Channel.byteOffset
                pos_byte_end[Channel.name] = pos_byte_beg[Channel.name] + Channel.nBytes_aligned
            index[record_id] = 0
            CGrecordLength[record_id] = record[record_id]['record'].CGrecordLength
            channel_name_set[record_id] = record[record_id]['record'].channelNames
    # read data
    if record_id_size == 1:
        while position < data_block_length:
            memcpy(&record_id_char, &bit_stream[position], 1)
            position, buf, VLSD, index = unsorted_read4(bit_stream, tmp, record_id_char, 1,
                                                        position, buf, VLSD, pos_byte_beg, pos_byte_end,
                                                        c_format_structure, index, CGrecordLength,
                                                        VLSD_flag, VLSD_CG_name, VLSD_CG_signal_data_type, channel_name_set)
    elif record_id_size == 2:
        while position < data_block_length:
            memcpy(&record_id_short, &bit_stream[position], 2)
            position, buf, VLSD, index = unsorted_read4(bit_stream, tmp, record_id_short, 2,
                                                        position, buf, VLSD, pos_byte_beg, pos_byte_end,
                                                        c_format_structure, index, CGrecordLength,
                                                        VLSD_flag, VLSD_CG_name, VLSD_CG_signal_data_type, channel_name_set)
    elif record_id_size == 3:
        while position < data_block_length:
            memcpy(&record_id_long, &bit_stream[position], 4)
            position, buf, VLSD, index = unsorted_read4(bit_stream, tmp, record_id_long, 4,
                                                        position, buf, VLSD, pos_byte_beg, pos_byte_end,
                                                        c_format_structure, index, CGrecordLength,
                                                        VLSD_flag, VLSD_CG_name, VLSD_CG_signal_data_type, channel_name_set)
    elif record_id_size == 4:
        while position < data_block_length:
            memcpy(&record_id_long_long, &bit_stream[position], 8)
            position, buf, VLSD, index = unsorted_read4(bit_stream, tmp, record_id_long_long, 8,
                                                        position, buf, VLSD, pos_byte_beg, pos_byte_end,
                                                        c_format_structure, index, CGrecordLength,
                                                        VLSD_flag, VLSD_CG_name, VLSD_CG_signal_data_type, channel_name_set)
    # changing from bytes type to desired type
    if buf:
        for name in buf.keys():
            buf[name] = buf[name].view(dtype=numpy_format[name])
    # convert list to array for VLSD only
    if VLSD:
        for channel_name in VLSD:
            VLSD[channel_name] = np.array(VLSD[channel_name])
        buf.update(VLSD)
    return buf

cdef inline unsorted_read4(const char* bit_stream, bytes tmp, record_id,
                           unsigned short record_id_size, unsigned long long position,
                           buf, VLSD, pos_byte_beg, pos_byte_end, c_format_structure,
                           index, CGrecordLength, VLSD_flag, VLSD_CG_name, VLSD_CG_signal_data_type,
                           channel_name_set):
    cdef unsigned long VLSDLen = 0
    if not VLSD_flag[record_id]:  # not VLSD CG)
        for channel_name in channel_name_set[record_id]:  # list of channel classes
            buf[channel_name][index[record_id]] = \
                   tmp[position + pos_byte_beg[channel_name]:position + pos_byte_end[channel_name]]
        index[record_id] += 1
        position += CGrecordLength[record_id]
    else:  # VLSD CG
        position += <unsigned long long> record_id_size
        memcpy(&VLSDLen, &bit_stream[position], 4)  # VLSD length
        position += 4
        signal_data_type = VLSD_CG_signal_data_type[record_id]
        temp = bytes(bit_stream[position:position + VLSDLen - 1])  # default: raw bytes
        if signal_data_type == 6:
            temp = temp.decode('ISO8859')
        elif signal_data_type == 7:
            temp = temp.decode('utf-8')
        elif signal_data_type == 8:
            temp = temp.decode('<utf-16')
        elif signal_data_type == 9:
            temp = temp.decode('>utf-16')
        VLSD[VLSD_CG_name[record_id]].append(temp)
        position += <unsigned long long> VLSDLen
    return position, buf, VLSD, index


def sd_data_read(unsigned short signal_data_type, bytes sd_block,
                 unsigned long long sd_block_length, unsigned long long n_records):
    """ Reads vlsd channel from its SD Block bytes

    Parameters
    ----------------
    signal_data_type : int

    sd_block : bytes
        SD Block bytes

    sd_block_length: int
        SD Block data length (header not included)

    n_records: int
        number of records

    Returns
    -----------
    array
    """
    if n_records == 0:
        return None
    cdef const char* bit_stream = PyBytes_AsString(sd_block)
    cdef unsigned long long buf_len = len(sd_block)
    cdef unsigned long max_len = 0
    cdef unsigned long vlsd_len = 0
    cdef unsigned long *VLSDLen = <unsigned long *> PyMem_Malloc(n_records * sizeof(unsigned long))
    cdef unsigned long long *pointer = <unsigned long long *> PyMem_Malloc(n_records * sizeof(unsigned long long))
    cdef unsigned long long rec = 0
    if not VLSDLen or not pointer:
        raise MemoryError()
    try:
        pointer[0] = 0
        VLSDLen[0] = 0
        for rec from 0 <= rec < n_records - 1 by 1:
            if pointer[rec] + 4 > buf_len:
                printf('VLSD sd_data_read: pointer out of bounds at record %llu\n', rec)
                return None
            memcpy(&vlsd_len, &bit_stream[pointer[rec]], 4)
            VLSDLen[rec] = vlsd_len
            pointer[rec + 1] = VLSDLen[rec] + 4 + pointer[rec]
            if VLSDLen[rec] > max_len:
                    max_len = VLSDLen[rec]
        if pointer[rec] + 4 > buf_len:
            printf('VLSD sd_data_read: pointer out of bounds at last record %llu\n', rec)
            return None
        memcpy(&vlsd_len, &bit_stream[pointer[rec]], 4)
        VLSDLen[rec] = vlsd_len
        if VLSDLen[rec] > max_len:
            max_len = VLSDLen[rec]
        if max_len != 0:
            if signal_data_type < 10 or signal_data_type == 17:
                if signal_data_type == 6:
                    channel_format = 'ISO8859'
                elif signal_data_type == 7:
                    channel_format = 'utf-8'
                elif signal_data_type == 8:
                    channel_format = '<utf-16'
                elif signal_data_type == 9:
                    channel_format = '>utf-16'
                elif signal_data_type == 17:
                    channel_format = 'bom'  # BOM-per-value detection
                else:
                    channel_format = 'utf-8'
                    printf('signal_data_type should have fixed length')
                return equalize_string_length(bit_stream, pointer, VLSDLen, max_len, n_records, channel_format)
            else:  # byte arrays or mime types
                return equalize_byte_length(bit_stream, pointer, VLSDLen, max_len, n_records)
        else:
            printf('VLSD channel could not be properly read\n')
            return None
    finally:
        PyMem_Free(pointer)
        PyMem_Free(VLSDLen)

cdef inline equalize_byte_length(const char* bit_stream, unsigned long long *pointer, unsigned long *VLSDLen,
                                 unsigned long max_len, unsigned long long n_records):
    cdef np.ndarray output = np.zeros((n_records, ), dtype='V{}'.format(max_len))
    cdef unsigned long rec = 0
    for rec from 0 <= rec < n_records by 1:  # resize string to same length, numpy constrain
        output[rec] = bytearray(bit_stream[pointer[rec]+4:pointer[rec]+4+VLSDLen[rec]]).rjust(max_len,  b'\x00')
    return output

cdef inline equalize_string_length(const char* bit_stream, unsigned long long *pointer, unsigned long *VLSDLen,
                                   unsigned long max_len, unsigned long long n_records, channel_format):
    cdef np.ndarray output = np.zeros((n_records, ), dtype='U{}'.format(max_len))
    cdef unsigned long rec = 0
    cdef bytes raw
    if channel_format == 'bom':
        # signal_data_type 17: detect BOM per value
        for rec from 0 <= rec < n_records by 1:
            raw = bytes(bit_stream[pointer[rec]+4:pointer[rec]+4+VLSDLen[rec]])
            if len(raw) >= 3 and raw[0] == 0xEF and raw[1] == 0xBB and raw[2] == 0xBF:
                output[rec] = raw[3:].decode('utf-8', errors='replace').rstrip('\x00')
            elif len(raw) >= 2 and raw[0] == 0xFF and raw[1] == 0xFE:
                output[rec] = raw[2:].decode('utf-16-le', errors='replace').rstrip('\x00')
            elif len(raw) >= 2 and raw[0] == 0xFE and raw[1] == 0xFF:
                output[rec] = raw[2:].decode('utf-16-be', errors='replace').rstrip('\x00')
            elif len(raw) > 0:
                output[rec] = raw.decode('utf-8', errors='replace').rstrip('\x00')
    else:
        for rec from 0 <= rec < n_records by 1:  # resize string to same length, numpy constrain
            output[rec] = bit_stream[pointer[rec]+4:pointer[rec]+4+VLSDLen[rec]].decode(channel_format).rstrip('\x00')
    return output


def vd_data_read(unsigned short signal_data_type, bytes vd_block,
                 object offsets_array, object sizes_array,
                 unsigned long long n_records):
    """Read VLSC channel data from raw VD block bytes.

    Parameters
    ----------------
    signal_data_type : int
        signal data type (6=ISO-8859-1, 7=UTF-8, 8=UTF-16-LE, 9=UTF-16-BE,
        17=BOM, 10+=byte array)
    vd_block : bytes
        raw VD block data (no per-value length prefix, unlike SD blocks)
    offsets_array : numpy uint64 array
        byte offset of each value within vd_block
    sizes_array : numpy uint64 array
        byte size of each value
    n_records : int
        number of records

    Returns
    -------
    numpy array of decoded strings (dtype U...) or byte arrays (dtype V...)
    """
    cdef const char* bit_stream = PyBytes_AsString(vd_block)
    cdef unsigned long long i
    cdef unsigned long long offset, size
    cdef unsigned long max_len = 0

    # Find max_len from sizes
    for i in range(n_records):
        size = <unsigned long long> sizes_array[i]
        if size > max_len:
            max_len = <unsigned long> size

    if max_len == 0:
        return None

    if signal_data_type < 10 or signal_data_type == 17:
        return vd_equalize_string(bit_stream, offsets_array, sizes_array,
                                  max_len, n_records, signal_data_type)
    else:
        return vd_equalize_bytes(bit_stream, offsets_array, sizes_array,
                                 max_len, n_records)


cdef inline vd_equalize_string(const char* bit_stream, object offsets_array, object sizes_array,
                                unsigned long max_len, unsigned long long n_records,
                                unsigned short signal_data_type):
    """Decode string values from VD block using offset/size pairs."""
    cdef np.ndarray output = np.zeros((n_records,), dtype='U{}'.format(max_len))
    cdef unsigned long long i
    cdef unsigned long long offset
    cdef unsigned long size
    cdef bytes raw
    if signal_data_type == 6:
        for i in range(n_records):
            size = <unsigned long> sizes_array[i]
            if size > 0:
                offset = <unsigned long long> offsets_array[i]
                output[i] = bit_stream[offset:offset+size].decode('ISO-8859-1', errors='replace').rstrip('\x00')
    elif signal_data_type == 7:
        for i in range(n_records):
            size = <unsigned long> sizes_array[i]
            if size > 0:
                offset = <unsigned long long> offsets_array[i]
                output[i] = bit_stream[offset:offset+size].decode('utf-8', errors='replace').rstrip('\x00')
    elif signal_data_type == 8:
        for i in range(n_records):
            size = <unsigned long> sizes_array[i]
            if size > 0:
                offset = <unsigned long long> offsets_array[i]
                output[i] = bit_stream[offset:offset+size].decode('utf-16-le', errors='replace').rstrip('\x00')
    elif signal_data_type == 9:
        for i in range(n_records):
            size = <unsigned long> sizes_array[i]
            if size > 0:
                offset = <unsigned long long> offsets_array[i]
                output[i] = bit_stream[offset:offset+size].decode('utf-16-be', errors='replace').rstrip('\x00')
    elif signal_data_type == 17:
        # BOM-per-value: detect encoding from BOM bytes
        for i in range(n_records):
            size = <unsigned long> sizes_array[i]
            if size > 0:
                offset = <unsigned long long> offsets_array[i]
                raw = bytes(bit_stream[offset:offset+size])
                if size >= 3 and raw[0] == 0xEF and raw[1] == 0xBB and raw[2] == 0xBF:
                    output[i] = raw[3:].decode('utf-8', errors='replace').rstrip('\x00')
                elif size >= 2 and raw[0] == 0xFF and raw[1] == 0xFE:
                    output[i] = raw[2:].decode('utf-16-le', errors='replace').rstrip('\x00')
                elif size >= 2 and raw[0] == 0xFE and raw[1] == 0xFF:
                    output[i] = raw[2:].decode('utf-16-be', errors='replace').rstrip('\x00')
                else:
                    output[i] = raw.decode('utf-8', errors='replace').rstrip('\x00')
    return output


cdef inline vd_equalize_bytes(const char* bit_stream, object offsets_array, object sizes_array,
                               unsigned long max_len, unsigned long long n_records):
    """Return byte-array values from VD block using offset/size pairs."""
    cdef np.ndarray output = np.zeros((n_records,), dtype='V{}'.format(max_len))
    cdef unsigned long long i
    cdef unsigned long long offset
    cdef unsigned long size
    for i in range(n_records):
        size = <unsigned long> sizes_array[i]
        if size > 0:
            offset = <unsigned long long> offsets_array[i]
            output[i] = bytearray(bit_stream[offset:offset+size]).rjust(max_len, b'\x00')
    return output
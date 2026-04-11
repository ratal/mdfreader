.. _performance:

Performance
===========

When the ``dataRead`` Cython extension is compiled, mdfreader applies three
complementary optimisations that together give a **3–5× speedup** on large
MDF4 files compared to the pure-Python path.

Benchmark (184 MB, 36 000-channel MDF4 file)
---------------------------------------------

.. list-table::
   :widths: 45 20 20
   :header-rows: 1

   * - Scenario
     - Wall time
     - Speedup
   * - Pure Python (no Cython)
     - ~1.9 s
     - 1×
   * - v4.2 with Cython (bit-exact reader only)
     - ~1.9 s
     - 1×
   * - v4.3 with all optimisations
     - **~0.6 s**
     - **~3×**

Optimisation 1 — Fast CN/CC/SI/TX metadata reader
--------------------------------------------------

**Function:** ``read_cn_chain_fast()`` in ``dataRead.pyx``

**Called from:** :meth:`~mdfreader.mdfinfo4.Info4.read_cn_blocks`

MDF4 channel metadata is organised as a singly-linked list of CN blocks.
For a 36 000-channel file this involves:

* ~575 000 ``fid.seek()`` + ``fid.read()`` Python call pairs
* ~36 000 ``struct.unpack()`` calls per channel
* ~36 000 ``lxml.objectify`` XML parses for unit/description MD blocks
* ~36 000 SI-block reads (often duplicated — same source for many channels)

``read_cn_chain_fast()`` replaces all of the above with a single C-level loop:

``pread()``
    POSIX atomic offset-based read (``<unistd.h>``).  No ``seek()``, no GIL
    acquired during I/O, no Python file-object method dispatch.  Each block
    is read in one syscall directly into a C stack buffer.

C packed structs
    All MDF4 block fields are extracted with ``memcpy`` into typed
    ``cdef packed struct`` variables matching the on-disk little-endian
    layout.  No ``struct.unpack()`` tuple allocation, no intermediate
    ``bytes`` object.

Fast ``<TX>`` scan
    Unit and description text lives in MD blocks (XML) of the form
    ``<CNcomment><TX>text</TX></CNcomment>``.  A plain ``bytes.find()`` scan
    for ``<TX>`` / ``</TX>`` handles >95 % of real files without invoking
    the XML parser.  ``lxml`` is only called for CDATA sections or
    namespace-qualified elements.

SI-block cache
    Source Information blocks are shared across many channels (e.g. all
    channels from one CAN bus share the same SI block).  Results are stored
    in ``Info4._si_cache`` (keyed by file offset) so each unique source is
    read and decoded only once.

Unchanged paths (Python)
    * ``_unique_channel_name()`` — requires full ``Info4`` context.
    * Composition blocks (CA/CN/DS/CL/CV/CU) — recursive, rare.
    * CC val/ref for complex types (cc\_type 3, 7–11) — variable-length
      link arrays; handled by :class:`~mdfreader.mdfinfo4.CCBlock`.

Optimisation 2 — SymBufReader (bidirectional file buffer)
---------------------------------------------------------

**Class:** ``SymBufReader`` in ``dataRead.pyx``

**Called from:** :meth:`~mdfreader.mdfinfo4.Info4.__init__`

MDF4 metadata blocks are written in reverse order (newest first) and linked
by backward-pointing file offsets.  Python's default ``BufferedReader``
fills its buffer *forward* from the current position, so every backward seek
past the buffer boundary causes a kernel ``read()`` call.

``SymBufReader`` keeps a **64 KB buffer centred on the current position**
(mirroring the Rust `mdfr` library design).  Seeks within ±32 KB of the
last fill position are served from the C-level ``unsigned char`` buffer
without any kernel interaction.  The buffer is stored as a C array (not a
Python ``bytes`` object) to avoid allocation and GC pressure.

.. code-block:: python

   # Transparent drop-in for fid:
   from dataRead import SymBufReader
   reader = SymBufReader(fid)
   reader.seek(offset)
   data = reader.read(88)   # served from buffer if within ±32 KB

Optimisation 3 — Single-call record read (sorted data)
------------------------------------------------------

**Function:** :meth:`~mdfreader.mdf4reader.Mdf4Record.read_all_channels_sorted_record`

For sorted channel groups (the common case), all records for the group are
read in **one** ``readinto()`` call into a pre-allocated flat ``uint8``
buffer.  The buffer is then reinterpreted as a structured NumPy record array
via ``.view()``.  This eliminates the per-chunk Python loop that previously
issued thousands of ``fid.read()`` calls for large data blocks.

.. code-block:: python

   raw = numpy.empty(total_size, dtype='u1')
   fid.readinto(raw)                          # single syscall
   return raw.view(record_dtype).view(recarray)[:n_records]

Building the Cython extension
-----------------------------

The extension is built automatically by ``setup.py`` when Cython and NumPy
are available::

   pip install cython numpy
   python setup.py build_ext --inplace

If compilation fails, mdfreader falls back to the pure-Python path (with
``bitarray`` for bit-exact channel reading).  A ``ImportWarning`` is issued
at import time to alert you.

To verify the fast path is active:

.. code-block:: python

   from mdfreader.mdfinfo4 import _CN_CHAIN_FAST, _SYMBUF_AVAILABLE
   print('fast CN reader:', _CN_CHAIN_FAST)
   print('SymBufReader:', _SYMBUF_AVAILABLE)

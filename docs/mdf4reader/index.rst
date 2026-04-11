mdf4reader — MDF4 sample-data reader
=====================================

This module reads the *sample data* from MDF 4.x files.
Metadata is parsed separately by :mod:`mdfreader.mdfinfo4`.

Key classes
-----------

:class:`~mdfreader.mdf4reader.Mdf4`
    Version-4 reader.  Inherits from :class:`~mdfreader.mdf.MdfSkeleton`
    and is mixed into :class:`~mdfreader.mdfreader.Mdf`.

:class:`~mdfreader.mdf4reader.Mdf4Record`
    Per-channel-group record layout.  Built once during reading; holds
    dtype, byte offsets and channel references.  The
    :meth:`~mdfreader.mdf4reader.Mdf4Record.read_all_channels_sorted_record`
    method reads all records in a single ``readinto()`` call.

Data block types
----------------

.. list-table::
   :widths: 15 85
   :header-rows: 1

   * - Block
     - Description
   * - DT
     - Plain data (sorted records concatenated)
   * - DZ
     - Deflate/transposed-deflate or LZ4/ZStandard compressed data
   * - HL
     - Header list grouping DL blocks
   * - DL
     - Data list: ordered sequence of DT/DZ blocks for one channel group
   * - SD
     - Signal data for VLSD channels (variable-length strings/arrays)
   * - DV/DI
     - Data Values / Data Invalidation (unsorted records with record IDs)

Conversion types
----------------

Raw-to-physical conversion (``cc_type`` in the CC block) is applied after
reading when ``convert_after_read=True`` (default):

.. list-table::
   :widths: 10 20 70
   :header-rows: 1

   * - Type
     - Name
     - Description
   * - 0
     - Identity
     - No conversion; raw value is the physical value
   * - 1
     - Linear
     - ``phys = a0 + a1 * raw``
   * - 2
     - Rational
     - Six-coefficient rational function
   * - 3
     - Formula
     - Algebraic formula (requires ``sympy``)
   * - 4
     - Tab (interpolated)
     - Piecewise-linear table with interpolation
   * - 5
     - Tab (no interpolation)
     - Step table
   * - 6
     - Range → value
     - Maps value ranges to scalar outputs
   * - 7
     - Value → text
     - Maps discrete values to strings
   * - 8
     - Range → text
     - Maps value ranges to strings
   * - 9–11
     - Text ↔ value/text
     - Text input or output conversions

.. automodule:: mdfreader.mdf4reader
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

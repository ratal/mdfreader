channel — per-channel record layout
=====================================

This module provides :class:`~mdfreader.channel.Channel4` and
:class:`~mdfreader.channel.Channel3`, which encapsulate the record-level
layout of a single channel within its channel group.

These objects are created once during reading and cached in the
:class:`~mdfreader.mdf4reader.Mdf4Record` /
:class:`~mdfreader.mdf3reader.Mdf3Record` structures.  They expose helper
methods used by the data readers:

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Method
     - Returns
   * - ``unit(info)``
     - Physical unit string from CC or CN block
   * - ``desc(info)``
     - Channel description from CN comment block
   * - ``conversion(info)``
     - CC block dict for raw→physical conversion
   * - ``calc_bytes(info)``
     - Aligned byte width of channel in the record
   * - ``set(info)``
     - Populate all attributes from the ``Info4`` dict

.. automodule:: mdfreader.channel
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

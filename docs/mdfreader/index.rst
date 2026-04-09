mdfreader — top-level API
=========================

This module provides the two public entry points:

:class:`~mdfreader.mdfreader.Mdf`
    The main class for reading, writing, converting and exporting MDF files.
    Inherits from :class:`~mdfreader.mdf4reader.Mdf4` and
    :class:`~mdfreader.mdf3reader.Mdf3`, which in turn inherit from
    :class:`~mdfreader.mdf.MdfSkeleton`.

    Typical usage::

       import mdfreader
       mdf = mdfreader.Mdf('file.mf4')
       data = mdf.get_channel_data('speed')
       mdf.resample(0.01)
       mdf.export_to_csv('output.csv')

:class:`~mdfreader.mdfreader.MdfInfo`
    Read-only parser that extracts the block structure of an MDF file
    without loading sample data.  Useful for listing channels, inspecting
    metadata, or checking file integrity before a full read.

    Typical usage::

       info = mdfreader.MdfInfo('file.mf4')
       channels = info.list_channels()

.. automodule:: mdfreader.mdfreader
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

mdf3reader — MDF3 sample-data reader
=====================================

This module reads sample data from MDF 3.x files.
Metadata is parsed by :mod:`mdfreader.mdfinfo3`.

MDF3 differences from MDF4
---------------------------

* No compressed data blocks (DZ/DL/HL are MDF4 features)
* Channel conversion types differ — cc\_type 0 is linear (not identity)
* No VLSD or VLSC string channels
* No Source Information (SI) blocks

.. automodule:: mdfreader.mdf3reader
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

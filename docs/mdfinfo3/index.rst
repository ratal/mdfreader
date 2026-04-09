mdfinfo3 — MDF3 metadata parser
================================

This module parses the block structure of MDF 3.x files.
The result is stored in an :class:`~mdfreader.mdfinfo3.Info3` nested dict,
analogous to :class:`~mdfreader.mdfinfo4.Info4` for MDF4.

MDF3 uses different block names than MDF4:

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Key
     - Description
   * - ``info['HDBlock']``
     - Header block
   * - ``info['DGBlock'][dg]``
     - Data group block
   * - ``info['CGBlock'][dg][cg]``
     - Channel group block
   * - ``info['CNBlock'][dg][cg][cn]``
     - Channel block
   * - ``info['CCBlock'][dg][cg][cn]``
     - Channel conversion block

.. automodule:: mdfreader.mdfinfo3
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

mdfinfo4 — MDF4 metadata parser
================================

This module parses the *block structure* of MDF 4.x files and stores the
result in the nested :class:`~mdfreader.mdfinfo4.Info4` dict.  It does not
read sample data (that is handled by :mod:`mdfreader.mdf4reader`).

Reading paths
-------------

Two paths exist for the CN/CC/SI/TX metadata hot loop.  The fast path is
selected automatically when the ``dataRead`` Cython extension is available.
See :doc:`../performance` for a full description.

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Path
     - Description
   * - **Fast (Cython)**
     - ``read_cn_chain_fast()`` — ``pread()`` + C packed structs, ``<TX>``
       bytes scan, SI cache.  Active when ``_CN_CHAIN_FAST`` is ``True``.
   * - **Fallback (Python)**
     - :meth:`~mdfreader.mdfinfo4.Info4.read_cn_block` — ``struct.unpack``
       via :class:`~mdfreader.mdfinfo4.CNBlock` /
       :class:`~mdfreader.mdfinfo4.CCBlock` / :class:`~mdfreader.mdfinfo4.SIBlock`.

Info4 dict structure
--------------------

.. code-block:: python

   info = mdfreader.MdfInfo('file.mf4')

   info['HD']                         # header block
   info['DG'][dg]                     # data group block
   info['CG'][dg][cg]                 # channel group block
   info['CN'][dg][cg][cn]             # channel block
   info['CC'][dg][cg][cn]             # conversion block
   info['AT']                         # attachment blocks
   info['VLSD'][dg][cg]               # list of VLSD channel keys
   info['VLSD_CG'][record_id]         # VLSD CG → (cg, cn)
   info['MLSD'][dg][cg]               # MLSD channel cross-reference
   info['masters']                    # master channel registry
   info['allChannelList']             # set of all channel names
   info['ChannelNamesByDG'][dg]       # set of channel names per data group

Block parsers
-------------

.. automodule:: mdfreader.mdfinfo4
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

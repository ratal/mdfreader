mdf — base class and channel dict
==================================

:mod:`mdfreader.mdf` provides :class:`~mdfreader.mdf.MdfSkeleton`, the base
class inherited by all version-specific readers.  It defines the channel dict
structure, all export methods, resampling, and plotting.

Channel dict layout
-------------------

After reading, the :class:`~mdfreader.mdf.MdfSkeleton` object (a ``dict``
subclass) stores one entry per channel:

.. code-block:: python

   mdf['channel_name'] = {
       'data':        numpy_array,
       'unit':        'm/s',
       'master':      'time_master',
       'masterType':  1,        # 1=time, 2=angle, 3=distance, 4=index
       'description': 'Vehicle speed',
       # present only when convert_after_read=False:
       'conversion':  {'type': 1, 'parameters': {'cc_val': [0.0, 1.0]}},
       'id':          ((dg, cg, cn), (name, source, path), (grp, gs, gp)),
   }

Key constants (importable from :mod:`mdfreader.mdf`):

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Constant
     - Dict key it names
   * - ``dataField``
     - ``'data'``
   * - ``unitField``
     - ``'unit'``
   * - ``masterField``
     - ``'master'``
   * - ``masterTypeField``
     - ``'masterType'``
   * - ``descriptionField``
     - ``'description'``
   * - ``conversionField``
     - ``'conversion'``
   * - ``idField``
     - ``'id'``

.. automodule:: mdfreader.mdf
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

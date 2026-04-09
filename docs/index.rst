mdfreader |version| documentation
===================================

**mdfreader** is a Python library for reading and writing MDF files
(Measured Data Format, versions 3.x and 4.x), the standard binary format
used in automotive ECU data logging (INCA, CANape, CANoe).

.. code-block:: python

   import mdfreader

   # Read all channels
   mdf = mdfreader.Mdf('recording.mf4')
   time  = mdf.get_channel_data('t')
   speed = mdf.get_channel_data('vehicle_speed')

   # Inspect the file structure without loading data
   info = mdfreader.MdfInfo('recording.mf4')

   # List channels grouped by master (raster)
   print(mdf.masterChannelList)

   # Resample and export
   mdf.resample(sampling=0.01)
   mdf.export_to_csv('output.csv')

Quick-start
-----------

Install from PyPI (binary wheel includes the Cython extension)::

   pip install mdfreader

Build from source (recommended for development)::

   pip install cython numpy
   python setup.py build_ext --inplace
   pip install -e .

.. seealso::

   :ref:`genindex`, :ref:`modindex`

Architecture
------------

mdfreader is split into two layers per MDF version:

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Module
     - Responsibility
   * - :mod:`mdfreader.mdfreader`
     - Top-level :class:`~mdfreader.mdfreader.Mdf` and
       :class:`~mdfreader.mdfreader.MdfInfo` classes; dispatch by file version
   * - :mod:`mdfreader.mdf`
     - :class:`~mdfreader.mdf.MdfSkeleton` base class; channel dict structure,
       export methods, resampling, plotting
   * - :mod:`mdfreader.mdfinfo4`
     - MDF4 metadata parser (:class:`~mdfreader.mdfinfo4.Info4`); reads all
       block types (HD/DG/CG/CN/CC/SI/TX/MD/AT/EV) into a nested dict
   * - :mod:`mdfreader.mdf4reader`
     - MDF4 sample-data reader; sorted/unsorted records, DT/DZ/DL/HL blocks,
       VLSD/VLSC strings, channel conversion
   * - :mod:`mdfreader.mdfinfo3`
     - MDF3 metadata parser
   * - :mod:`mdfreader.mdf3reader`
     - MDF3 sample-data reader
   * - :mod:`mdfreader.channel`
     - :class:`~mdfreader.channel.Channel4` /
       :class:`~mdfreader.channel.Channel3` — per-channel record layout helper
   * - ``dataRead`` (Cython)
     - Low-level fast helpers: :func:`read_cn_chain_fast`,
       :class:`SymBufReader`, :func:`sorted_data_read`, bit-exact decoders.
       See :doc:`performance`.

Channel dict structure
----------------------

After reading, each channel in the :class:`~mdfreader.mdfreader.Mdf` dict
has the following keys:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Key
     - Description
   * - ``data``
     - :class:`numpy.ndarray` of channel samples
   * - ``unit``
     - Physical unit string (e.g. ``'m/s'``)
   * - ``master``
     - Name of the master channel for this raster
   * - ``masterType``
     - Master type: ``1``\ =time, ``2``\ =angle, ``3``\ =distance, ``4``\ =index
   * - ``description``
     - Free-text channel description
   * - ``conversion``
     - Raw→physical conversion dict (present when ``convert_after_read=False``)
   * - ``id``
     - ``(dg, cg, cn)`` index tuple into the :class:`~mdfreader.mdfinfo4.Info4` dict

Contents
--------

.. toctree::
   :maxdepth: 2

   performance
   mdfreader/index
   mdf/index
   mdf4reader/index
   mdfinfo4/index
   mdf3reader/index
   mdfinfo3/index
   channel/index

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

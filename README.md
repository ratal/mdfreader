**MDFREADER**
**************

Abstract:
=========
This Module imports MDF files (Measured Data Format V3.x and V4.x), typically from INCA (ETAS), CANape or CANOe. It is widely used in automotive industry to record data from ECUs. The main module mdfreader.py inherits from 2 modules (One pair for each MDF version X) : The first one to read the file's blocks descriptions (mdfinfoX) and the second (mdfXreader) to read the raw data from the file. It can optionally run multithreaded. It was built in mind to process efficently big amount of data in a batch, generally data from endurance evaluation.

The structure of the mdf object inheriting from python dict
===========================================================
for each channel: mdf[channelName] below keys exist
* data: numpy array
* unit: unitName
* master : vector name corresponding to master channel
* masterType : type of master channel (time, angle, distance, etc.)
* description : physical meaning of channel
* conversion: (exist when reading with convertAfterRead=False) dictionary describing how to convert raw data into meaningful/physical data

mdf object main attribute: masterChannelList, a dict containing one list of channel names per datagroup


Mdfreader module methods:
=========================
* resample channels to one sampling frequency
* merge files
* plot one or a list of channels

It is also possible to export mdf data into:
* CSV file (excel dialect by default)
* NetCDF file for a compatibility with Uniplot for instance (needs netcdf4, Scientific.IO)
* HDF5 (needs h5py)
* Excel 95 to 2003 (needs xlwt, really slooow, be careful about data size)
* Excel 2007/2010 (needs openpyxl, slow if not resampled data)
* Matlab .mat (needs scipy.io)
* MDF simplified file. It allows you to modify data, units, description and save it again in the same major mdf version (3.x->3.3, 4.x-> 4.1)
* Pandas dataframe(s) (only in command line, not in mdfconverter). One dataframe per raster.

Compatibility:
==============
This code is compatible for both python 2.6-7 and python 3.2+
Evaluated for Windows and Linux platforms (x86 and AMD64)

Requirements:
=============
Mdfreader is mostly relying on numpy/scipy/matplotlib.

Reading channels defined by a formula will require sympy.

Cython is required to compile dataRead module for reading quickly exotic data (not byte aligned or containing hidden bytes). However, if cython compilation fails, bitarray becomes required (slower, pure python).

Export requirements (optional): scipy, csv, h5py, xlwt(3), openpyxl, pandas

Mdfconverter user interface requires PyQt

Installation:
=============
pip package existing:
```shell
pip install mdfreader
```
or with only source from github from instance
```python
python setup.py develop
```

Graphical interface: mdfconverter (PyQt4&PyQt5)
==================================
User interface in PyQt4 or PyQt5 to convert batch of files is part of package. You can launch it with command 'mdfconverter'. By right clicking a channel in the interface list, you can plot it. You can also drag-drop channels between columns to tune import list. Channel list from a .lab text file can be imported. You can optionally merge several files into one and even resample all of them.

Others:
=======
In the case of big files and lack of memory, you can optionally:
* Read only a channel list (slightly slower, argument channelList = ['channel', 'list'])
* Keep raw data as stored in mdf without data type conversion (mdfreader argument convertAfterRead=False). Data will then be converted on the fly by the other functions (plot, exportTo..., getChannelData, etc.) but raw data type will remain as in mdf file along with conversion information.

For great data visualization, dataPlugin for Veusz (from 1.16, http://home.gna.org/veusz/) is also existing ; please follow instructions from Veusz documentation and plugin file's header.

Warning:
========
MDF 4.x specification is much complex compared to 3.x and its implementation is young and not 100% complete. Chances of bug are higher with version 4.x compared to 3.x

Command example in ipython:
===========================
```python
    import mdfreader
    # loads whole mdf file content in yop mdf object
    yop=mdfreader.mdf('NameOfFile')
    # alternatively, for max speed and smallest memory footprint
    yop=mdfreader.mdf('NameOfFile',channelList=['channel', 'list'],convertAfterRead=False)
    # to get file mdf versoin
    yop.MDFVersionNumber
    # to only get file struture, you can instead create a mdfinfo instance
    info=mdfreader.mdfinfo()
    info.listChannels('NameOfFile') # returns the list of channels
    info.readinfo('NameOfFile') # more complete file structure object
    # to list channels names after reading
    yop.keys()
    # to list channels names grouped by raster, below dict mdf attirbute contains pairs (key=masterChannelName : value=listOfChannelNamesForThisMaster)
    yop.masterChannelList
    # quick plot of channel(s)
    yop.plot('channelName') or yop.plot({'channel1','channel2'})
    # file manipulations
    yop.resample(0.1) or yop.resample(channelName='master3')
    yop.exportoCSV(sampling=0.01)
    yop.exportNetCDF()
    yop.exporttoHDF5()
    yop.exporttoMatlab()
    # converts data groups into pandas dataframes
    yop.convertToPandas()
    # drops all the channels except the one in argument
    yop.keepChannels({'channel1','channel2','channel3'})
    # can write mdf file after modifications (by default, same version of orignal file)
    yop.write()
    # to get/show raw data from channel after read
    yop.getChannelData('channelName') # returns channel numpy array
```
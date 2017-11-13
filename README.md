**MDFREADER**
**************

Abstract:
=========
This Module imports MDF files (Measured Data Format V3.x and V4.x), typically from INCA (ETAS), CANape or CANOe. It is widely used in automotive industry to record data from ECUs. The main module mdfreader.py inherits from 2 modules (One pair for each MDF version X) : The first one to read the file's blocks descriptions (mdfinfoX) and the second (mdfXreader) to read the raw data from the file. It can optionally run multithreaded. It was built in mind to process efficently big amount of data in a batch, endurance evaluation files for data mining.

The structure of the mdf object inheriting from python dict
===========================================================
for each channel: mdf[channelName] below keys exist
* data: numpy array
* unit: unit name
* master : master channel name of channelName
* masterType : type of master channel (time, angle, distance, etc.)
* description : description of channel
* conversion: (exist when reading with convertAfterRead=False) dictionary describing how to convert raw data into meaningful/physical data

mdf object main attribute: masterChannelList, a dict containing a list of channel names per datagroup


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
* MDF file. It allows you to create or modify data, units, description and save it again.
* Pandas dataframe(s) (only in command line, not in mdfconverter). One dataframe per raster.

Compatibility:
==============
This code is compatible for both python 2.7 and python 3.4+
Evaluated for Windows and Linux platforms (x86 and AMD64)

Requirements:
=============
Mdfreader is mostly relying on numpy/scipy/matplotlib.

Reading channels defined by a formula will require sympy.

Cython is required to compile dataRead module for reading quickly exotic data (not byte aligned or containing hidden bytes) or only a list of channels. However, if cython compilation fails, bitarray becomes required (slower, pure python and maybe not so robust as not so much tested).

Export requirements (optional): scipy, csv, h5py, xlwt(3), openpyxl, pandas

Mdfconverter graphical user interface requires PyQt (versions 4 or 5)

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

Graphical interface: mdfconverter (PyQt4 and PyQt5)
==================================
User interface in PyQt4 or PyQt5 to convert batch of files is part of package. You can launch it with command 'mdfconverter' from shell. By right clicking a channel in the interface list, you can plot it. You can also drag-drop channels between columns to tune import list. Channel list from a .lab text file can be imported. You can optionally merge several files into one and even resample all of them.

Others:
=======
In the case of big files or lack of memory, you can optionally:
* Read only a channel list (argument channelList = ['channel', 'list'])
* Keep raw data as stored in mdf without data type conversion (argument convertAfterRead=False). Data will then be converted on the fly by the other functions (plot, exportTo..., getChannelData, etc.) but raw data type will remain as in mdf file along with conversion information.
* Compress data in memory with blosc or bcolz with argument compression. If integer or boolean is given, it will use by default bcolz with integer compression level. If 'blosc' is given, default compression level is 9.
* Create a mdf dict with its metadata but without raw data (argument noDataLoading=True). Data will be loaded on demand by mdfreader methods (in general by getChannelData method)

For great data visualization, dataPlugin for Veusz (from 1.16, http://home.gna.org/veusz/) is also existing ; please follow instructions from Veusz documentation and plugin file's header.

Command example in ipython:
===========================
```python
    import mdfreader
    # loads whole mdf file content in yop mdf object
    yop=mdfreader.mdf('NameOfFile')
    # alternatively, for max speed and smaller memory footprint, read only few channels
    yop=mdfreader.mdf('NameOfFile',channelList=['channel1', 'channel2'],convertAfterRead=False)
    # also possible to keep data compressed for small memory footprint
    yop=mdfreader.mdf('NameOfFile',compression=True)
    # for interactive file exploration, possible to read the file but not its data to save memory
    yop=mdfreader.mdf('NameOfFile',noDataLoading=True) # channel data will be loaded from file if needed
    yop.getChannel('channelName') # to yield one channel and keep its content in mdf object
    # to get file mdf version
    yop.MDFVersionNumber
    # to get file structure or attachments, you can create a mdfinfo instance
    info=mdfreader.mdfinfo()
    info.listChannels('NameOfFile') # returns only the list of channels
    info.readinfo('NameOfFile') # complete file structure object
    yop.info # same class is stored in mdfreader class
    # to list channels names after reading
    yop.keys()
    # to list channels names grouped by raster, below dict mdf attribute contains
    # pairs (key=masterChannelName : value=listOfChannelNamesForThisMaster)
    yop.masterChannelList
    # quick plot of channel(s)
    yop.plot('channelName') or yop.plot({'channel1','channel2'})
    # file manipulations
    yop.resample(0.1) or yop.resample(masterChannel='master3')
    yop.cut(begin=10, end=15)  # keep only data between begin and end
    yop.exporToCSV(sampling=0.01)
    yop.exportNetCDF()
    yop.exportToHDF5()
    yop.exportToMatlab()
    # converts data groups into pandas dataframes
    yop.convertToPandas()
    # drops all the channels except the one in argument
    yop.keepChannels({'channel1','channel2','channel3'})
    # merge 2 files
    yop2=mdfreader.mdf('NameOfFile_2')
    yop=mergeMDF(yop2)
    # can write mdf file after modifications
    yop.write()  # write in same version as original file
    yop.write4()  # write mdf version 4 file
    yop.write3()  # write mdf version 3 file
    # to get/show raw data from channel after read
    yop.getChannelData('channelName') # returns channel numpy array
```
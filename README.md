**MDFREADER**
**************

Abstract:
=========
This Module imports MDF files (Measured Data Format V3.x and V4.x), typically from INCA (ETAS), CANape or CANOe. It is widely used in automotive industry to record data from ECUs. The main module mdfreader.py inherits from 2 modules (One pair for each MDF version X) : The first one to read the file's blocks descriptions (mdfinfoX) and the second (mdfXreader) to read the raw data from the file. It can optionally run multithreaded. It was built in mind to process efficiently big amount of data in a batch, endurance evaluation files for data mining.

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
* plot one channel, several channels on one graph (list) or several channels on subplots (list of lists)

It is also possible to export mdf data into:
* CSV file (excel dialect by default)
* NetCDF file for a compatibility with Uniplot for instance (needs netcdf4, Scientific.IO)
* HDF5 (needs h5py)
* Excel 95 to 2003 (needs xlwt, extremely slooow, be careful about data size)
* Excel 2007/2010 (needs openpyxl, can be also slow with big files)
* Matlab .mat (needs hdf5storage)
* MDF file. It allows you to create, convert or modify data, units, description and save it again.
* Pandas dataframe(s) (only in command line, not in mdfconverter). One dataframe per raster.

Compatibility:
==============
This code is compatible for python 3.4+
Evaluated for Windows and Linux platforms (x86 and AMD64)

Requirements:
=============
Mdfreader is mostly relying on numpy/scipy/matplotlib and lxml for parsing the metadata in mdf version 4.x files

Reading channels defined by a formula will require sympy.

Cython is strongly advised and allows to compile dataRead module for reading quickly exotic data (not byte aligned or containing hidden bytes) or only a list of channels. However, if cython compilation fails, bitarray becomes required (slower, pure python and maybe not so robust as not so much tested).

Export requirements (optional): scipy, csv, h5py, hdf5storage, xlwt(3), openpyxl, pandas

Blosc for data compression (optional)

Mdfconverter graphical user interface requires PyQt (versions 4 or 5)

Installation:
=============
pip package existing:
```shell
pip install mdfreader
```
or from source cloned from github from instance
```shell
python setup.py develop
```

Graphical interface: mdfconverter (PyQt4 and PyQt5)
==================================
User interface in PyQt4 or PyQt5 to convert batch of files is part of package. You can launch it with command 'mdfconverter' from shell. By right clicking a channel in the interface list, you can plot it. You can also drag-drop channels between columns to tune import list. Channel list from a .lab text file can be imported. You can optionally merge several files into one and even resample all of them.

Others:
=======
In the case of big files or lack of memory, you can optionally:
* Read only a channel list (argument channel_list = ['channel', 'list'], you can get the file channel list without loading data with mdfinfo)
* Keep raw data as stored in mdf without data type conversion (argument convert_after_read=False). Data will then be converted on the fly by the other functions (plot, export_to..., get_channel_data, etc.) but raw data type will remain as in mdf file along with conversion information.
* Compress data in memory with blosc with argument compression. Default compression level is 9.
* Create a mdf dict with its metadata but without data (argument no_data_loading=True). Data will be read from file on demand by mdfreader methods (in general by get_channel_data method)

For great data visualization, dataPlugin for Veusz (from 1.16, http://home.gna.org/veusz/) is also existing ; please follow instructions from Veusz documentation and plugin file's header.

Command example in ipython:
===========================
```python
    import mdfreader
    # loads whole mdf file content in yop mdf object.
    yop=mdfreader.Mdf('NameOfFile')
    # you can print file content in ipython with a simple:
    yop
    # alternatively, for max speed and smaller memory footprint, read only few channels
    yop=mdfreader.Mdf('NameOfFile', channel_list=['channel1', 'channel2'], convert_after_read=False)
    # also possible to keep data compressed for small memory footprint, using Blosc module
    yop=mdfreader.Mdf('NameOfFile', compression=True)
    # for interactive file exploration, possible to read the file but not its data to save memory
    yop=mdfreader.Mdf('NameOfFile', no_data_loading=True) # channel data will be loaded from file if needed
    # parsing xml metadata from mdf4.x for many channels can take more than just reading data.
    # You can reduce to minimum metadata reading with below argument (no source information, attachment, etc.) 
    yop=mdfreader.Mdf('NameOfFile', metadata=0)  # 0: full, 2: minimal
    # only for mdf4.x, you can search for the mdf key of a channel name that can have been recorded by different sources
    yop.get_channel_name4('channelName', 'source path or name')  # returns list of mdf keys
    # to yield one channel and keep its content in mdf object
    yop.get_channel('channelName')
    # to yield one channel numpy array
    yop.get_channel_data('channelName')
    # to get file mdf version
    yop.MDFVersionNumber
    # to get file structure or attachments, you can create a mdfinfo instance
    info=mdfreader.MdfInfo()
    info.list_channels('NameOfFile') # returns only the list of channels
    info.read_info('NameOfFile') # complete file structure object
    yop.info # same class is stored in mdfreader class
    # to list channels names after reading
    yop.keys()
    # to list channels names grouped by raster, below dict mdf attribute contains
    # pairs (key=masterChannelName : value=listOfChannelNamesForThisMaster)
    yop.masterChannelList
    # quick plot or subplot (with lists) of channel(s)
    yop.plot(['channel1',['channel2','channel3']])
    # file manipulations
    yop.resample(0.1)
    # or
    yop.resample(master_channel='master3')
    # keep only data between begin and end
    yop.cut(begin=10, end=15)
    # export to other file formats :
    yop.export_to_csv(sampling=0.01)
    yop.export_to_NetCDF()
    yop.export_to_hdf5()
    yop.export_to_matlab()
    yop.export_to_xlsx()
    yop.export_to_parquet()
    # return pandas dataframe from master channel name
    yop.return_pandas_dataframe('master_channel_name')
    # converts data groups into pandas dataframes and keeps it in mdf object
    yop.convert_to_pandas()
    # drops all the channels except the one in argument
    yop.keep_channels({'channel1','channel2','channel3'})
    # merge 2 files
    yop2=mdfreader.Mdf('NameOfFile_2')
    yop.merge_mdf(yop2)
    # can write mdf file after modifications or creation from scratch
    # write4 and write3 also allow to convert file versions
    yop.write('NewNameOfFile')  # write in same version as original file after modifications
    yop.write4('NameOfFile', compression=True)  # write mdf version 4.1 file, data compressed
    yop.write3()  # write mdf version 3 file
    yop.attachments  # to get attachments, embedded or paths to files 
```
<a href="https://scan.coverity.com/projects/ratal-mdfreader">
  <img alt="Coverity Scan Build Status"
       src="https://scan.coverity.com/projects/21368/badge.svg"/>
</a>
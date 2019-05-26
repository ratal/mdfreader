# -*- coding: utf-8 -*-
""" Measured Data Format blocks parser for version 3.x

Created on Thu Dec 9 12:57:28 2014

:Author: `Aymeric Rateau <https://github.com/ratal/mdfreader>`__


Dependencies
-------------------
- Python >2.6, >3.2 <http://www.python.org>
- Numpy >1.14 <http://numpy.scipy.org>


Attributes
--------------
PythonVersion : float
    Python version currently running, needed for compatibility of both python 2.6+ and 3.2+


mdfinfo3
--------------------------
"""
from __future__ import absolute_import  # for consistency between python 2 and 3
from __future__ import print_function
from sys import version_info
from warnings import warn
from numpy import sort, zeros
from struct import unpack, Struct
from .mdf import dataField, descriptionField, unitField, masterField, masterTypeField, idField
PythonVersion = version_info
PythonVersion = PythonVersion[0]

cn_struct = Struct('<2sH5IH32s128s4H3d2IH')
tx_struct = Struct('<2sH')
cc_struct = Struct('<2s2H2d20s2H')
dg_struct = Struct('<2sH4I2H')
cg_struct = Struct('<2sH3I3HI')
ce_struct = Struct('<2s2H')


class Info3(dict):
    __slots__ = ['fileName', 'fid', 'filterChannelNames']
    """ mdf file info class version 3.x
    MDFINFO is a class information about an MDF (Measure Data Format) file
    Based on following specification http://powertrainnvh.com/nvh/MDFspecificationv03.pdf

    Attributes
    --------------
    filterChannelNames : bool, optional
        flag to filter long channel names including module names separated by a '.'
    fileName : str
        name of file
    fid :
        file identifier

    Notes
    --------
    mdfinfo(FILENAME) contains a dict of structures, for
    each data group, containing key information about all channels in each
    group. FILENAME is a string that specifies the name of the MDF file.
    General dictionary structure is the following

    - mdfinfo['HDBlock'] header block
    - mdfinfo['DGBlock'][dataGroup] Data Group block
    - mdfinfo['CGBlock'][dataGroup][channelGroup] Channel Group block
    - mdfinfo['CNBlock'][dataGroup][channelGroup][channel] Channel block including text blocks for comment and identifier
    - mdfinfo['CCBlock'][dataGroup][channelGroup][channel] Channel conversion information
    """

    def __init__(self, file_name=None, fid=None, filter_channel_names=False, minimal=0):
        """ info3 class constructor

        Parameters
        ----------------
        file_name : str, optional
            name of file
        fid : float, optional
            file identifier
        filter_channel_names : bool, optional
            flag to filter long channel names including module names separated by a '.'
        minimal : int
            0 will load every metadata
            1 will load DG, CG, CN and CC
            2 will load only DG

        Notes
        --------
        If fileName is given it will read file blocks directly by calling method readinfo3
        """
        self['IDBlock'] = dict()  # Identifier Block
        self['HDBlock'] = dict()  # Header Block
        self['DGBlock'] = dict()  # Data Group Block
        self['CGBlock'] = dict()  # Channel Group Block
        self['CNBlock'] = dict()  # Channel Block
        self['CCBlock'] = dict()  # Conversion block
        self['allChannelList'] = set()  # all channels
        self['ChannelNamesByDG'] = dict()
        self.filterChannelNames = filter_channel_names
        self.fileName = file_name
        self.fid = None
        if file_name is not None and fid is None:
            try:
                self.fid = open(self.fileName, 'rb')
            except IOError:
                raise IOError('Can not find file ' + self.fileName)
            self.read_info3(self.fid, minimal)
        elif file_name is None and fid is not None:
            self.read_info3(fid, minimal)

    def read_info3(self, fid, minimal=0):
        """ read all file blocks except data

        Parameters
        ----------------
        fid : float
            file identifier
        minimal : int
            0 will load every metadata
            1 will load DG, CG, CN and CC
            2 will load only DG
        """
        # reads IDBlock
        fid.seek(24)
        (self['IDBlock']['ByteOrder'],
         self['IDBlock']['FloatingPointFormat'],
         self['IDBlock']['id_ver'],
         self['IDBlock']['CodePageNumber']) = unpack('<4H', fid.read(8))

        # Read header block (HDBlock) information
        # Read Header block info into structure, HD pointer at 64 as mentioned in specification
        self['HDBlock'] = read_hd_block(fid, 64, self['IDBlock']['id_ver'])

        # Read text block (TXBlock) information
        self['HDBlock']['TXBlock'] = read_tx_block(fid, self['HDBlock']['pointerToTXBlock'])

        # Read text block (PRBlock) information
        self['HDBlock']['PRBlock'] = read_tx_block(fid, self['HDBlock']['pointerToPRBlock'])

        # Read Data Group blocks (DGBlock) information
        # Get pointer to first Data Group block
        dg_pointer = self['HDBlock']['pointerToFirstDGBlock']
        for dg in range(self['HDBlock']['numberOfDataGroups']):

            # Read data Data Group block info into structure
            self['DGBlock'][dg] = read_dg_block(fid, dg_pointer)
            # Get pointer to next Data Group block
            dg_pointer = self['DGBlock'][dg]['pointerToNextDGBlock']
            self['ChannelNamesByDG'][dg] = set()
            if minimal < 2:
                self.read_cg_block(fid, dg, minimal)

        # Close the file
        fid.close()

    def read_cg_block(self, fid, dg, minimal=0):
        """ read all CG blocks and relying CN & CC

        Parameters
        ----------------
        fid : float
            file identifier
        dg : int
            data group number
        minimal : int
            0 will load every metadata
            1 will load DG, CG, CN and CC
            2 will load only DG
        """
        # Read Channel Group block (CGBlock) information - offset set already

        # Read data Channel Group block info into structure
        cg_pointer = self['DGBlock'][dg]['pointerToNextCGBlock']
        self['CGBlock'][dg] = dict()
        self['CNBlock'][dg] = dict()
        self['CCBlock'][dg] = dict()
        for cg in range(self['DGBlock'][dg]['numberOfChannelGroups']):
            self['CNBlock'][dg][cg] = dict()
            self['CCBlock'][dg][cg] = dict()
            self['CGBlock'][dg][cg] = read_cg_block(fid, cg_pointer)
            cg_pointer = self['CGBlock'][dg][cg]['pointerToNextCGBlock']

            # Read data Text block info into structure
            self['CGBlock'][dg][cg]['TXBlock'] = \
                read_tx_block(fid, self['CGBlock'][dg][cg]['pointerToChannelGroupCommentText'])

            # Get pointer to next first Channel block
            cn_pointer = self['CGBlock'][dg][cg]['pointerToFirstCNBlock']

            # For each Channel
            for channel in range(self['CGBlock'][dg][cg]['numberOfChannels']):

                # Read data Channel block info into structure

                self['CNBlock'][dg][cg][channel] = read_cn_block(fid, cn_pointer)
                cn_pointer = self['CNBlock'][dg][cg][channel]['pointerToNextCNBlock']

                # Read Channel text blocks (TXBlock)

                # Clean signal name
                short_signal_name = self['CNBlock'][dg][cg][channel]['signalName']  # short name of signal
                # keep original non unique channel name
                self['CNBlock'][dg][cg][channel]['orig_name'] = short_signal_name

                if self['CNBlock'][dg][cg][channel]['pointerToASAMNameBlock'] > 0:
                    # long name of signal
                    long_signal_name = \
                        read_tx_block(fid, self['CNBlock'][dg][cg][channel]['pointerToASAMNameBlock'])
                    if len(long_signal_name) > len(short_signal_name):  # long name should be used
                        signal_name = long_signal_name
                    else:
                        signal_name = short_signal_name
                else:
                    signal_name = short_signal_name
                signal_name = signal_name.split('\\')
                if len(signal_name) > 1:
                    self['CNBlock'][dg][cg][channel]['deviceName'] = signal_name[1]
                signal_name = signal_name[0]
                if self.filterChannelNames:
                    signal_name = signal_name.split('.')[-1]  # filters channels modules

                if signal_name in self['ChannelNamesByDG'][dg]:  # for unsorted data
                    pointer = self['CNBlock'][dg][cg][channel]['pointerToCEBlock']
                    if pointer:
                        temp = read_ce_block(fid, pointer)
                    else:
                        temp = dict()
                        temp['tail'] = channel
                    self['CNBlock'][dg][cg][channel]['signalName'] = \
                        '{0}_{1}_{2}_{3}'.format(signal_name, dg, cg, temp['tail'])
                elif signal_name in self['allChannelList']:
                    # doublon name or master channel
                    pointer = self['CNBlock'][dg][cg][channel]['pointerToCEBlock']
                    if pointer:
                        temp = read_ce_block(fid, pointer)
                    else:
                        temp = dict()
                        temp['tail'] = dg
                    if '{0}_{1}'.format(signal_name, temp['tail']) not in self['allChannelList']:
                        self['CNBlock'][dg][cg][channel]['signalName'] = '{0}_{1}'.format(signal_name, temp['tail'])
                    else:
                        self['CNBlock'][dg][cg][channel]['signalName'] = '{0}_{1}_{2}'.format(signal_name,
                                                                                              dg, temp['tail'])
                else:
                    self['CNBlock'][dg][cg][channel]['signalName'] = signal_name
                self['ChannelNamesByDG'][dg].add(self['CNBlock'][dg][cg][channel]['signalName'])
                self['allChannelList'].add(self['CNBlock'][dg][cg][channel]['signalName'])

                # Read channel description
                self['CNBlock'][dg][cg][channel]['ChannelCommentBlock'] = \
                    read_tx_block(fid, self['CNBlock'][dg][cg][channel]['pointerToChannelCommentBlock'])

                # Read data Text block info into structure
                self['CNBlock'][dg][cg][channel]['SignalIdentifierBlock'] = \
                    read_tx_block(fid, self['CNBlock'][dg][cg][channel]['pointerToSignalIdentifierBlock'])

                # Read Channel Conversion block (CCBlock)
                if self['CNBlock'][dg][cg][channel]['pointerToConversionFormula'] == 0:
                    # If no conversion formula, set to 1:1
                    self['CCBlock'][dg][cg][channel] = {}
                    self['CCBlock'][dg][cg][channel]['cc_type'] = 65535
                else:  # Otherwise get conversion formula, parameters and physical units
                    # Read data Channel Conversion block info into structure
                    self['CCBlock'][dg][cg][channel] = \
                        read_cc_block(fid, self['CNBlock'][dg][cg][channel]['pointerToConversionFormula'])

            # reorder channel blocks and related blocks based on first bit position
            # this reorder is meant to improve performance while parsing records using core.records.fromfile
            # as it will not use cn_byte_offset
            # first, calculate new mapping/order
            n_channel = len(self['CNBlock'][dg][cg])
            Map = zeros(shape=len(self['CNBlock'][dg][cg]), dtype=[('index', 'u4'), ('first_bit', 'u4')])
            for cn in range(n_channel):
                Map[cn] = (cn, self['CNBlock'][dg][cg][cn]['numberOfTheFirstBits'])
            ordered_map = sort(Map, order='first_bit')

            to_change_index = Map == ordered_map
            for cn in range(n_channel):
                if not to_change_index[cn]:
                    # offset all indexes of indexes to be moved
                    self['CNBlock'][dg][cg][cn + n_channel] = self['CNBlock'][dg][cg].pop(cn)
                    self['CCBlock'][dg][cg][cn + n_channel] = self['CCBlock'][dg][cg].pop(cn)
            for cn in range(n_channel):
                if not to_change_index[cn]:
                    # change to ordered index
                    self['CNBlock'][dg][cg][cn] = self['CNBlock'][dg][cg].pop(ordered_map[cn][0] + n_channel)
                    self['CCBlock'][dg][cg][cn] = self['CCBlock'][dg][cg].pop(ordered_map[cn][0] + n_channel)

    def clean_dg_info(self, dg):
        """ delete CN,CC and CG blocks related to data group

        Parameters
        ----------------
        dg : int
            data group number
        """
        try:
            self['CNBlock'][dg] = {}
        except KeyError:
            pass
        try:
            self['CCBlock'][dg] = {}
        except KeyError:
            pass
        try:
            self['CGBlock'][dg] = {}
        except KeyError:
            pass

    def list_channels3(self, file_name=None, fid=None):
        """ reads data, channel group and channel blocks to list channel names

        Attributes
        --------------
        file_name : str
            file name

        Returns
        -----------
        list of channel names
        """
        # Read MDF file and extract its complete structure
        if file_name is not None:
            self.fileName = file_name
        # Open file
        if fid is None and file_name is not None:
            fid = open(self.fileName, 'rb')
        channel_name_list = []

        # Read header block (HDBlock) information
        # Set file pointer to start of HDBlock
        hd_pointer = 64
        # Read Header block info into structure
        self['HDBlock'] = read_hd_block(fid, hd_pointer)

        # Read Data Group blocks (DGBlock) information
        # Get pointer to first Data Group block
        dg_pointer = self['HDBlock']['pointerToFirstDGBlock']
        for dataGroup in range(self['HDBlock']['numberOfDataGroups']):

            # Read data Data Group block info into structure
            self['DGBlock'][dataGroup] = read_dg_block(fid, dg_pointer)
            # Get pointer to next Data Group block
            dg_pointer = self['DGBlock'][dataGroup]['pointerToNextDGBlock']

            # Read Channel Group block (CGBlock) information - offset set already

            # Read data Channel Group block info into structure
            cg_pointer = self['DGBlock'][dataGroup]['pointerToNextCGBlock']
            self['CGBlock'][dataGroup] = {}
            self['CNBlock'][dataGroup] = {}
            self['CCBlock'][dataGroup] = {}
            for channelGroup in range(self['DGBlock'][dataGroup]['numberOfChannelGroups']):
                self['CNBlock'][dataGroup][channelGroup] = {}
                self['CCBlock'][dataGroup][channelGroup] = {}
                self['CGBlock'][dataGroup][channelGroup] = read_cg_block(fid, cg_pointer)
                cg_pointer = self['CGBlock'][dataGroup][channelGroup]['pointerToNextCGBlock']

                # Get pointer to next first Channel block
                cn_pointer = self['CGBlock'][dataGroup][channelGroup]['pointerToFirstCNBlock']

                names_set = set()
                # For each Channel
                for channel in range(self['CGBlock'][dataGroup][channelGroup]['numberOfChannels']):

                    # Read Channel block (CNBlock) information
                    # self.numberOfChannels += 1
                    # Read data Channel block info into structure
                    self['CNBlock'][dataGroup][channelGroup][channel] = read_cn_block(fid, cn_pointer)
                    cn_pointer = self['CNBlock'][dataGroup][channelGroup][channel]['pointerToNextCNBlock']

                    # Read Channel text blocks (TXBlock)

                    # Clean signal name
                    short_signal_name = self['CNBlock'][dataGroup][channelGroup][channel]['signalName']
                    cn_tx_pointer = self['CNBlock'][dataGroup][channelGroup][channel]['pointerToASAMNameBlock']
                    if cn_tx_pointer > 0:
                        long_signal_name = read_tx_block(fid, cn_tx_pointer)  # long name of signal
                        if len(long_signal_name) > len(short_signal_name):  # long name should be used
                            signal_name = long_signal_name
                        else:
                            signal_name = short_signal_name
                    else:
                        signal_name = short_signal_name
                    signal_name = signal_name.split('\\')
                    signal_name = signal_name[0]
                    if self.filterChannelNames:
                        signal_name = signal_name.split('.')[-1]

                    if signal_name in names_set:
                        pointer = self['CNBlock'][dataGroup][channelGroup][channel]['pointerToCEBlock']
                        if pointer:
                            temp = read_ce_block(fid, pointer)
                        else:
                            temp = dict()
                            temp['tail'] = channel
                        self['CNBlock'][dataGroup][channelGroup][channel]['signalName'] = \
                            '{0}_{1}'.format(signal_name, temp['tail'])
                        warn('WARNING added number to duplicate signal name: {}'.
                             format(self['CNBlock'][dataGroup][channelGroup][channel]['signalName']))
                    else:
                        self['CNBlock'][dataGroup][channelGroup][channel]['signalName'] = signal_name
                        names_set.add(signal_name)

                    channel_name_list.append(signal_name)

        # CLose the file
        fid.close()
        return channel_name_list


def read_hd_block(fid, pointer, version=0):
    """ header block reading """
    if pointer != 0 and pointer is not None:
        temp = dict()
        fid.seek(pointer)
        if version < 320:
            (temp['BlockType'],
             temp['BlockSize'],
             temp['pointerToFirstDGBlock'],
             temp['pointerToTXBlock'],
             temp['pointerToPRBlock'],
             temp['numberOfDataGroups'],
             temp['Date'],
             temp['Time'],
             temp['Author'],
             temp['Organization'],
             temp['ProjectName'],
             temp['Subject']) = unpack('<2sH3IH10s8s32s32s32s32s', fid.read(164))
        else:
            (temp['BlockType'],
             temp['BlockSize'],
             temp['pointerToFirstDGBlock'],
             temp['pointerToTXBlock'],
             temp['pointerToPRBlock'],
             temp['numberOfDataGroups'],
             temp['Date'],
             temp['Time'],
             temp['Author'],
             temp['Organization'],
             temp['ProjectName'],
             temp['Subject'],
             temp['TimeStamp'],
             temp['UTCTimeOffset'],
             temp['TimeQualityClass'],
             temp['TimeIdentification']) = unpack('<2sH3IH10s8s32s32s32s32sQ2H32s', fid.read(208))
            temp['TimeIdentification'] = temp['TimeIdentification'].rstrip(b'\x00').decode('latin1', 'replace')
        temp['Date'] = temp['Date'].rstrip(b'\x00').decode('latin1', 'replace')
        temp['Time'] = temp['Time'].rstrip(b'\x00').decode('latin1', 'replace')
        temp['Author'] = temp['Author'].rstrip(b'\x00').decode('latin1', 'replace')
        temp['Organization'] = temp['Organization'].rstrip(b'\x00').decode('latin1', 'replace')
        temp['ProjectName'] = temp['ProjectName'].rstrip(b'\x00').decode('latin1', 'replace')
        temp['Subject'] = temp['Subject'].rstrip(b'\x00').decode('latin1', 'replace')
        return temp
    else:
        return None


def read_dg_block(fid, pointer):
    """ data group block reading """
    if pointer != 0 and pointer is not None:
        temp = dict()
        fid.seek(pointer)
        (temp['BlockType'],
         temp['BlockSize'],
         temp['pointerToNextDGBlock'],
         temp['pointerToNextCGBlock'],
         temp['reserved'],
         temp['pointerToDataRecords'],
         temp['numberOfChannelGroups'],
         temp['numberOfRecordIDs']) = dg_struct.unpack(fid.read(24))
        return temp
    else:
        return None


def read_cg_block(fid, pointer):
    """ channel block reading """
    if pointer != 0 and pointer is not None:
        temp = dict()
        fid.seek(pointer)
        (temp['BlockType'],
         temp['BlockSize'],
         temp['pointerToNextCGBlock'],
         temp['pointerToFirstCNBlock'],
         temp['pointerToChannelGroupCommentText'],
         temp['recordID'],
         temp['numberOfChannels'],
         temp['dataRecordSize'],
         temp['numberOfRecords']) = cg_struct.unpack(fid.read(26))
        return temp
    else:
        return None


def read_cn_block(fid, pointer):
    """ channel block reading """
    if pointer != 0 and pointer is not None:
        temp = dict()
        fid.seek(pointer)
        (temp['BlockType'],
         temp['BlockSize'],
         temp['pointerToNextCNBlock'],
         temp['pointerToConversionFormula'],
         temp['pointerToCEBlock'],
         temp['pointerToCDBlock'],
         temp['pointerToChannelCommentBlock'],
         temp['channelType'],
         temp['signalName'],
         temp['signalDescription'],
         temp['numberOfTheFirstBits'],
         temp['numberOfBits'],
         temp['signalDataType'],
         temp['valueRangeKnown'],
         temp['valueRangeMinimum'],
         temp['valueRangeMaximum'],
         temp['rateVariableSampled'],
         temp['pointerToASAMNameBlock'],
         temp['pointerToSignalIdentifierBlock'],
         temp['ByteOffset']) = cn_struct.unpack(fid.read(228))

        temp['signalName'] = temp['signalName'].rstrip(b'\x00').decode('latin1', 'replace')
        temp['signalDescription'] = temp['signalDescription'].rstrip(b'\x00').decode('latin1', 'replace')
        return temp
    else:
        return None


def read_cc_block(fid, pointer):
    """ channel conversion block reading """
    if pointer != 0 and pointer is not None:
        temp = dict()
        fid.seek(pointer)
        (temp['BlockType'],
         temp['BlockSize'],
         temp['valueRangeKnown'],
         temp['valueRangeMinimum'],
         temp['valueRangeMaximum'],
         temp['physicalUnit'],
         temp['cc_type'],
         temp['numberOfValuePairs']) = cc_struct.unpack(fid.read(46))
        temp['physicalUnit'] = temp['physicalUnit'].rstrip(b'\x00').decode('latin1', 'replace')
        temp['conversion'] = dict()

        if temp['cc_type'] == 0:  # Parametric, Linear: Physical = Integer*P2 + P1
            (temp['conversion']['P1'],
             temp['conversion']['P2']) = unpack('2d', fid.read(16))

        elif temp['cc_type'] in (1, 2):  # Table look up with or without interpolation
            for pair in range(temp['numberOfValuePairs']):
                temp['conversion'][pair] = dict()
                (temp['conversion'][pair]['int'],
                 temp['conversion'][pair]['phys']) = unpack('2d', fid.read(16))

        elif temp['cc_type'] in (6, 9):  # Polynomial or rational
            (temp['conversion']['P1'],
             temp['conversion']['P2'],
             temp['conversion']['P3'],
             temp['conversion']['P4'],
             temp['conversion']['P5'],
             temp['conversion']['P6']) = unpack('6d', fid.read(48))

        elif temp['cc_type'] in (7, 8):  # Exponential or logarithmic
            (temp['conversion']['P1'],
             temp['conversion']['P2'],
             temp['conversion']['P3'],
             temp['conversion']['P4'],
             temp['conversion']['P5'],
             temp['conversion']['P6'],
             temp['conversion']['P7']) = unpack('7d', fid.read(56))

        elif temp['cc_type'] == 10:  # Text Formula
            temp['conversion']['textFormula'] = \
                unpack('256s', fid.read(256))[0].rstrip(b'\x00').decode('latin1', 'replace')

        elif temp['cc_type'] == 11:  # ASAM-MCD2 text table
            for pair in range(temp['numberOfValuePairs']):
                temp['conversion'][pair] = dict()
                (temp['conversion'][pair]['int'],
                 temp['conversion'][pair]['text']) = unpack('d32s', fid.read(40))
                temp['conversion'][pair]['text'] = \
                    temp['conversion'][pair]['text'].rstrip(b'\x00').decode('latin1', 'replace')

        elif temp['cc_type'] == 12:  # Text range table
            for pair in range(temp['numberOfValuePairs']):
                temp['conversion'][pair] = dict()
                (temp['conversion'][pair]['lowerRange'],
                 temp['conversion'][pair]['upperRange'],
                 temp['conversion'][pair]['pointerToTXBlock']) = unpack('2dI', fid.read(20))
            for pair in range(temp['numberOfValuePairs']):  # get text range
                # Read corresponding text
                try:
                    temp['conversion'][pair]['Textrange'] = read_tx_block(fid, temp['conversion'][pair]['pointerToTXBlock'])
                except:
                    temp['conversion'][pair]['Textrange'] = ""

        elif temp['cc_type'] == 65535:  # No conversion int=phys
            pass
        else:
            # Give warning that conversion formula is not being
            # made
            warn('Conversion Formula type (cc_type={})not supported.'.format(temp['cc_type']))

        return temp
    else:
        return None


def read_tx_block(fid, pointer):
    """ reads text block """
    if pointer != 0 and pointer is not None:
        fid.seek(pointer)
        (block_type,
         block_size) = tx_struct.unpack(fid.read(4))
        text = unpack('{}s'.format(block_size - 4), fid.read(block_size - 4))[0]

        return text.rstrip(b'\x00').decode('latin1', 'replace')


def read_ce_block(fid, pointer):
    """ reads source block """
    if pointer != 0 and pointer is not None:
        fid.seek(pointer)
        temp = dict()
        (block_type,
         block_size,
         temp['extension_type']) = ce_struct.unpack(fid.read(6))
        if temp['extension_type'] == 2:
            (temp['n_module'],
             temp['address'],
             temp['description'],
             temp['ECU_id']) = unpack('HI80s32s', fid.read(120))
            temp['ECU_id'] = temp['ECU_id'].rstrip(b'\x00').decode('latin1', 'replace')
            temp['description'] = temp['description'].rstrip(b'\x00').decode('latin1', 'replace')
            temp['tail'] = temp['ECU_id']
        elif temp['extension_type'] == 19:
            (temp['CAN_id'],
             temp['CAN_channel_id'],
             temp['message'],
             temp['sender']) = unpack('2I36s36s', fid.read(80))
            temp['message'] = temp['message'].rstrip(b'\x00').decode('latin1', 'replace')
            temp['sender'] = temp['sender'].rstrip(b'\x00').decode('latin1', 'replace')
            temp['tail'] = temp['message']
        else:
            return None
        return temp


def _generate_dummy_mdf3(info, channel_list):
    """ computes MasterChannelList and mdf dummy dict from an info object

    Parameters
    ----------------
    info : info object
        information structure of file

    channel_list : list of str
        list of channel names

    Returns
    -----------
    a dict which keys are master channels in files with values a list of related channels of the raster
    """
    master_channel_list = {}
    all_channel_list = set()
    mdf_dict = {}
    for dg in info['DGBlock']:
        master = ''
        master_type = 0
        channel_names_by_dg = set()
        for cg in info['CGBlock'][dg]:
            channel_name_list = []
            for cn in info['CNBlock'][dg][cg]:
                name = info['CNBlock'][dg][cg][cn]['signalName']
                if name in channel_names_by_dg:
                    name = '{0}_{1}_{2}_{3}'.format(name, dg, cg, cn)
                elif name in all_channel_list:
                    name = ''.join([name, '_{}'.format(dg)])
                if channel_list is None or name in channel_list:
                    channel_name_list.append(name)
                    all_channel_list.add(name)
                    channel_names_by_dg.add(name)
                    # create mdf channel
                    mdf_dict[name] = {}
                    mdf_dict[name][dataField] = None
                    if 'signalDescription' in info['CNBlock'][dg][cg][cn]:
                        mdf_dict[name][descriptionField] = \
                                info['CNBlock'][dg][cg][cn]['signalDescription']
                    else:
                        mdf_dict[name][descriptionField] = ''
                    if 'physicalUnit' in info['CCBlock'][dg][cg][cn]:
                        mdf_dict[name][unitField] = info['CCBlock'][dg][cg][cn]['physicalUnit']
                    else:
                        mdf_dict[name][unitField] = ''
                    mdf_dict[name][masterField] = 0  # default is time
                    mdf_dict[name][masterTypeField] = None
                    mdf_dict[name][idField] = (dg, cg, cn)
                if info['CNBlock'][dg][cg][cn]['channelType']:  # master channel of cg
                    master = name
                    master_type = info['CNBlock'][dg][cg][cn]['channelType']
            for chan in channel_name_list:
                mdf_dict[chan][masterField] = master
                mdf_dict[chan][masterTypeField] = master_type
        master_channel_list[master] = channel_name_list
    return (master_channel_list, mdf_dict)

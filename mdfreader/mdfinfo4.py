# -*- coding: utf-8 -*-
""" Measured Data Format blocks paser for version 4.x

Platform and python version
----------------------------------------
With Unix and Windows for python 2.6+ and 3.2+

Created on Sun Dec 15 12:57:28 2013

:Author: `Aymeric Rateau <https://github.com/ratal/mdfreader>`__

Dependencies
-------------------
- Python >2.6, >3.2 <http://www.python.org>
- Numpy >1.6 <http://numpy.scipy.org>

Attributes
--------------
PythonVersion : float
    Python version currently running, needed for compatibility of both
    python 2.6+ and 3.2+

mdfinfo4 module
--------------------------
"""
from __future__ import print_function

from struct import calcsize, unpack, pack, Struct
from os.path import dirname, abspath
from os import remove
from sys import version_info, stderr, path
from numpy import sort, zeros
import time
from xml.etree.ElementTree import Element, SubElement, \
    tostring, register_namespace
from lxml import objectify
_root = dirname(abspath(__file__))
path.append(_root)
from mdf import _open_MDF, dataField, descriptionField, unitField, \
    masterField, masterTypeField, _convertName

PythonVersion = version_info
PythonVersion = PythonVersion[0]

# datatypes
LINK = '<Q'
REAL = '<d'
BOOL = '<h'
UINT8 = '<B'
BYTE = '<c'
INT16 = '<h'
UINT16 = '<H'
UINT32 = '<I'
INT32 = '<i'
UINT64 = '<Q'
INT64 = '<q'
CHAR = '<c'

HeaderStruct = Struct('<4sI2Q')
DZStruct = Struct('2s2BI2Q')
HLStruct = Struct('<QHB5s')
SIStruct = Struct('<4sI5Q3B5s')
DGStruct = Struct('<4sI6QB7s')
CGStruct = Struct('<4sI10Q2H3I')
CNStruct = Struct('<4B4IBcH6d')
CCStruct1 = Struct('<4sI6Q')
CCStruct2 = Struct('<2B3H2d')
SRStruct = Struct('<4sI5Qd2B6s')


def _loadHeader(fid, pointer):
    """ reads block's header and put in class dict

    Parameters
    ----------------
    fid : float
        file identifier
    pointer : int
        position of block in file
    """
    # All blocks have the same header
    if pointer != 0 and pointer is not None:
        fid.seek(pointer)
        temp = dict()
        (temp['id'],
         reserved,
         temp['length'],
         temp['link_count']) = HeaderStruct.unpack(fid.read(24))
        temp['pointer'] = pointer
        return temp
    else:
        return None


def _mdfblockread(fid, Type, count):
    """ converts a byte array of length count to a given data Type

    Parameters
    ----------------
    Type : str
        C format data type
    count : int
        number of elements to sequentially read

    Returns
    -----------
    array of values of 'Type' parameter
    """
    value = fid.read(calcsize(Type) * count)
    if value:
        if count == 1:
            return unpack(Type, value)[0]
        else:
            if '<' in Type or '>' in Type:
                endian = Type[0]
                Type = Type.strip('<>')
                return unpack(endian+count*Type, value)
            else:
                return unpack(count*Type, value)
    else:
        return None


def _writeHeader(fid, Id, block_length, link_count):
    """ Writes header of a block

    Parameters
    ----------------
    fid : float
        file identifier
    Id : str
        4 character id of block, for instance '##HD'
    block_length : int
        total block length
    link_count : int
        number of links in the block

    Returns
    -------
    (block_length_pointer, link_count_pointer)
    """
    # make sure beginning of block starts with a multiple of 8 position
    current_position = fid.tell()
    remain = current_position % 8
    if not remain == 0:
        current_position = current_position - remain + 8
        fid.seek(current_position)
    head = (Id, 0, block_length, link_count)
    fid.write(HeaderStruct.pack(*head))
    return current_position


def _writePointer(fid, pointer, value):
    """ Writes a value at pointer position and comes back to original position

    Parameters
    ----------------
    fid : float
        file identifier
    pointer : float
        pointer where to write value
    value : int
        value to write (LINK)
    """
    currentPosition = fid.tell()
    fid.seek(pointer)
    fid.write(pack(LINK, value))
    return fid.seek(currentPosition)


class IDBlock(dict):

    """ reads or writes ID Block
    """
    def __init__(self, fid=None):
        if fid is not None:
            self.read(fid)

    def read(self, fid):
        """ reads IDBlock
        """
        fid.seek(0)
        (self['id_file'],
         self['id_vers'],
         self['id_prog'],
         id_reserved1,
         self['id_ver'],
         id_reserved2,
         self['id_unfi_flags'],
         self['id_custom_unfi_flags']) = unpack('<8s8s8sIH30s2H',
                                                fid.read(64))
        # treatment of unfinalised file
        if self['id_ver'] > 410 and 'UnFin' in self['id_file']:
            print('  ! unfinalised file', file=stderr)
            if self['id_unfi_flags'] & 1:
                print('Update of cycle counters for CG/CA blocks required',
                      file=stderr)
            if self['id_unfi_flags'] & (1 << 1):
                print('Update of cycle counters for SR blocks required',
                      file=stderr)
            if self['id_unfi_flags'] & (1 << 2):
                print('Update of length for last DT block required',
                      file=stderr)
            if self['id_unfi_flags'] & (1 << 3):
                print('Update of length for last RD block required',
                      file=stderr)
            if self['id_unfi_flags'] & (1 << 4):
                print('Update of last DL block in each chained list of \
                      DL blocks required', file=stderr)
            if self['id_unfi_flags'] & (1 << 5):
                print('Update of cg_data_bytes and cg_inval_bytes in VLSD \
                      CG block required', file=stderr)
            if self['id_unfi_flags'] & (1 << 6):
                print('Update of offset values for VLSD channel required \
                      in case a VLSD CG block is used', file=stderr)

    def write(self, fid):
        """ Writes IDBlock
        """
        # MDF versionTxt tool reserved version_int
        head = (b'MDF     ', b'4.11    ', b'MDFreadr', b'\0' * 4, 411,
                b'\0' * 30, 0, 0)
        fid.write(pack('<8s8s8s4sH30s2H', *head))


class HDBlock(dict):

    """ reads Header block and save in class dict
    """

    def __init__(self, fid=None, pointer=64):
        if fid is not None:
            self.read(fid)

    def read(self, fid=None, pointer=64):
        fid.seek(pointer)
        (self['id'],
         reserved,
         self['length'],
         self['link_count'],
         self['hd_dg_first'],
         self['hd_fh_first'],
         self['hd_ch_first'],
         self['hd_at_first'],
         self['hd_ev_first'],
         self['hd_md_comment'],
         self['hd_start_time_ns'],
         self['hd_tz_offset_min'],
         self['hd_dst_offset_min'],
         self['hd_time_flags'],
         self['hd_time_class'],
         self['hd_flags'],
         hd_reserved,
         self['hd_start_angle_rad'],
         self['hd_start_distance']) = unpack('<4sI9Q2h4B2Q', fid.read(104))
        if self['hd_md_comment']:  # if comments exist
            self['Comment'] = {}
            self['Comment'].update(CommentBlock(fid, self['hd_md_comment'], 'HD'))

    def write(self, fid):
        # write block header
        currentPosition = _writeHeader(fid, b'##HD', 104, 6)
        # link section
        currentPosition += 24
        pointers = dict()
        pointers['HD'] = dict()
        pointers['HD']['DG'] = currentPosition
        pointers['HD']['FH'] = currentPosition + 8 * 1
        pointers['HD']['CH'] = currentPosition + 8 * 2
        pointers['HD']['AT'] = currentPosition + 8 * 3
        pointers['HD']['EV'] = currentPosition + 8 * 4
        pointers['HD']['MD'] = currentPosition + 8 * 5
        # (first Data group pointer, first file history block pointer,
        # pointer to hierarchy Block file, pointer to attachment Block,
        # pointer to event Block, pointer to comment Block,
        # time in ns, timezone offset in min, time daylight offest in min,
        # time flags, time class, hd flags, reserved, start angle in radians
        # start distance in meters)
        dataBytes = (0, 0, 0, 0, 0, 0,
                     int(time.time()*1E9),
                     int(time.timezone/60),
                     int(time.daylight*60),
                     2, 0, 0, b'\0', 0, 0)
        fid.write(pack('<7Q2h3Bs2d', *dataBytes))
        return pointers


class FHBlock(dict):

    """ reads File History block and save in class dict
    """

    def __init__(self, fid=None, pointer=None):
        if fid is not None:
            self.read(fid, pointer)

    def read(self, fid, pointer):
        fid.seek(pointer)
        (self['id'],
         reserved,
         self['length'],
         self['link_count'],
         self['fh_fh_next'],
         self['fh_md_comment'],
         self['fh_time_ns'],
         self['fh_tz_offset_min'],
         self['fh_dst_offset_min'],
         self['fh_time_flags'],
         fh_reserved) = unpack('<4sI5Q2hB3s', fid.read(56))
        if self['fh_md_comment']:  # comments exist
            self['Comment'] = {}
            self['Comment'].update(CommentBlock(fid, self['fh_md_comment'], 'FH'))

    def write(self, fid):
        # write block header
        currentPosition = _writeHeader(fid, b'##FH', 56, 2)
        # link section
        currentPosition += 24
        pointers = {}
        pointers['FH'] = dict()
        pointers['FH']['FH'] = currentPosition
        pointers['FH']['MD'] = currentPosition + 8
        # (No next FH, comment block pointer
        # time in ns, timezone offset in min, time daylight offest in min,
        # time flags, reserved)
        dataBytes = (0, 0,
                     int(time.time()*1E9),
                     int(time.timezone/60),
                     int(time.daylight*60),
                     2, b'\0' * 3)
        fid.write(pack('<3Q2hB3s', *dataBytes))
        return pointers


class CHBlock(dict):

    """ reads Channel Hierarchy block and saves in class dict
    """

    def __init__(self, fid, pointer):
        # block header
        self.update(_loadHeader(fid, pointer))
        # Channel hierarchy block
        (self['ch_ch_next'],
         self['ch_ch_first'],
         self['ch_tx_name'],
         self['ch_md_comment']) = unpack('<4Q', fid.read(32))
        nLinks = self['link_count'] - 4
        self['ch_element'] = unpack('<{}Q'.format(nLinks),
                                    fid.read(nLinks * 8))
        (self['ch_element_count'],
         self['ch_type'],
         ch_reserved) = unpack('<IB3s', fid.read(8))
        if self['ch_md_comment']:  # comments exist
            self['Comment'] = {}
            self['Comment'].update(CommentBlock(fid, self['ch_md_comment']))
        if self['ch_tx_name']:  # text block containing name of hierarchy level
            self['ch_name_level'] = {}
            self['ch_name_level'].update(CommentBlock(fid, self['ch_tx_name']))


class CommentBlock(dict):
    __slots__ = []
    """ reads or writes Comment block and saves in class dict
    """

    def __init__(self, fid=None, pointer=None, MDType=None):
        if fid is not None:
            self.read(fid, pointer, MDType)

    def read(self, fid, pointer, MDType=None):
        """ reads Comment block and saves in class dict

        Notes
        --------
        Can read xml (MD metadata) or text (TX) comments from several
        kind of blocks
        """

        if pointer > 0:
            fid.seek(pointer)
            # block header
            (self['id'],
             reserved,
             self['length'],
             self['link_count']) = HeaderStruct.unpack(fid.read(24))
            if self['id'] in ('##MD', b'##MD'):
                # Metadata block
                # removes normal 0 at end
                # self['Comment'] = None
                try:
                    xml_tree = objectify.fromstring(fid.read(self['length'] - 24).rstrip(b'\x00'))
                except:
                    print('xml metadata malformed', file=stderr)
                    xml_tree = None
                # specific action per comment block type,
                # #extracts specific tags from xml
                if MDType == 'CN':  # channel comment
                    try:
                        self['description'] = xml_tree.TX.text
                    except AttributeError:
                        print('Could not parse CN block TX tag')
                    try:
                        self['names'] = xml_tree.names.text
                    except AttributeError:
                        pass  # optional
                    try:
                        self['linker_name'] = xml_tree.linker_name.text
                    except AttributeError:
                        pass  # optional
                    try:
                        self['linker_address'] = xml_tree.linker_address.text
                    except AttributeError:
                        pass  # optional
                    try:
                        self['address'] = xml_tree.address.text
                    except AttributeError:
                        pass  # optional
                    try:
                        self['axis_monotony'] = xml_tree.axis_monotony.text
                    except AttributeError:
                        pass  # optional
                    try:
                        self['raster'] = xml_tree.raster.text
                    except AttributeError:
                        pass  # optional
                    try:
                        self['formula'] = xml_tree.formula.text
                    except AttributeError:
                        pass  # optional
                elif MDType == 'unit':  # channel comment
                    try:
                        self['unit'] = xml_tree.TX
                    except AttributeError:
                        print('Could not parse unit TX tag')
                elif MDType == 'HD':  # header comment
                    try:
                        self['TX'] = xml_tree.TX.text
                    except AttributeError:
                        print('Could not parse HD block TX tag')
                    try:
                        self['time_source'] = xml_tree.time_source.text
                    except AttributeError:
                        pass  # optional
                    try:
                        tmp = xml_tree.common_properties
                        for t in range(tmp.countchildren()):
                            self[tmp.e[t].attrib.values()[0]] = tmp.e[t].text
                    except AttributeError:
                        pass  # optional
                elif MDType == 'FH':  # File History comment
                    try:
                        self['TX'] = xml_tree.TX.text
                    except AttributeError:
                        print('Could not parse FH block TX tag')
                    try:
                        self['tool_id'] = xml_tree.tool_id.text
                    except AttributeError:
                        print('Could not parse FH block tool_id tag')
                    try:
                        self['tool_vendor'] = xml_tree.tool_vendor.text
                    except AttributeError:
                        print('Could not parse extract HD block tool_vendor tag')
                    try:
                        self['tool_version'] = xml_tree.tool_version.text
                    except AttributeError:
                        print('Could not parse extract HD block tool_vendor tag')
                    try:
                        self['user_name'] = xml_tree.user_name.text
                    except AttributeError:
                        pass  # optional
                elif MDType == 'SI':
                    try:
                        self['TX'] = xml_tree.TX.text
                    except AttributeError:
                        print('Could not parse SI block TX tag')
                    try:
                        self['names'] = xml_tree.names.text
                    except AttributeError:
                        pass  # optional
                    try:
                        self['path'] = xml_tree.path.text
                    except AttributeError:
                        pass  # optional
                    try:
                        self['bus'] = xml_tree.bus.text
                    except AttributeError:
                        pass  # optional
                    try:
                        self['protocol'] = xml_tree.protocol.text
                    except AttributeError:
                        pass  # optional
                elif MDType == 'CC':
                    try:
                        self['TX'] = xml_tree.TX.text
                    except AttributeError:
                        print('Could not parse CC block TX tag')
                    try:
                        self['names'] = xml_tree.names.text
                    except AttributeError:
                        pass  # optional
                    try:
                        self['COMPU_METHOD'] = xml_tree.COMPU_METHOD.text
                    except AttributeError:
                        pass  # optional
                    try:
                        self['formula'] = xml_tree.formula.text
                    except AttributeError:
                        pass  # optional
                else:
                    if MDType is not None:
                        print('No recognized MDType', file=stderr)
                        print(MDType, file=stderr)
            elif self['id'] in ('##TX', b'##TX'):
                if MDType == 'CN':  # channel comment
                    self['name'] = fid.read(self['length'] - 24).rstrip(b'\x00').decode('UTF-8', 'ignore')
                else:
                    self['Comment'] = fid.read(self['length'] - 24).rstrip(b'\x00').decode('UTF-8', 'ignore')


    def write(self, fid, data, MDType):

        if MDType == 'TX':
            data = b''.join([data.encode('utf-8', 'replace'), b'\0'])
            # make sure block is multiple of 8
            data_lentgth = len(data)
            remain = data_lentgth % 8
            if not remain == 0:
                data = b''.join([data, b'\0' * (8 - (remain % 8))])
            block_start = _writeHeader(fid, b'##TX', 24 + len(data), 0)
            fid.write(data)
        else:
            register_namespace('', 'http://www.asam.net/mdf/v4')
            if MDType == 'HD':
                root = Element('HDcomment')
                root.set('xmlns', 'http://www.asam.net/mdf/v4')
                TX = SubElement(root, 'TX')
                TX.text = data['comment']
                common_properties = SubElement(root, 'common_properties')
                e = SubElement(common_properties, 'e',
                               attrib={'name': 'subject'})
                e.text = data['subject']
                e = SubElement(common_properties, 'e',
                               attrib={'name': 'project'})
                e.text = data['project']
                e = SubElement(common_properties, 'e',
                               attrib={'name': 'department'})
                e.text = data['organisation']
                e = SubElement(common_properties, 'e',
                               attrib={'name': 'author'})
                e.text = data['author']
            elif MDType == 'CN':
                pass
            elif MDType == 'FH':
                root = Element('FHcomment')
                root.set('xmlns', 'http://www.asam.net/mdf/v4')
                TX = SubElement(root, 'TX')
                TX.text = data['comment']
                tool_id = SubElement(root, 'tool_id')
                tool_id.text = 'mdfreader'
                tool_vendor = SubElement(root, 'tool_vendor')
                tool_vendor.text = 'mdfreader is under GPL V3'
                tool_version = SubElement(root, 'tool_version')
                tool_version.text = '2.6'
            data = b''.join([tostring(root), b'\0'])
            # make sure block is multiple of 8
            data_lentgth = len(data)
            remain = data_lentgth % 8
            if not remain == 0:
                data = b''.join([data, b'\0' * (8 - (remain % 8))])
            _writeHeader(fid, b'##MD', 24 + len(data), 0)
            block_start = fid.write(data)
        return block_start


class DGBlock(dict):

    """ reads Data Group block and saves in class dict
    """

    def __init__(self, fid=None, pointer=None):
        if fid is not None:
            self.read(fid, pointer)

    def read(self, fid, pointer):
        fid.seek(pointer)
        self['pointer'] = pointer
        (self['id'],
         reserved,
         self['length'],
         self['link_count'],
         self['dg_dg_next'],
         self['dg_cg_first'],
         self['dg_data'],
         self['dg_md_comment'],
         self['dg_rec_id_size'],
         dg_reserved) = DGStruct.unpack(fid.read(64))
        if self['dg_md_comment']:  # comments exist
            self['Comment'] = {}
            self['Comment'].update(CommentBlock(fid, self['dg_md_comment']))

    def write(self, fid):
        pointers = {}
        # write block header
        currentPosition = _writeHeader(fid, b'##DG', 64, 4)
        pointers['block_start'] = currentPosition
        # link section
        currentPosition += 24
        pointers['DG'] = currentPosition
        pointers['CG'] = currentPosition + 8
        pointers['data'] = currentPosition + 8 * 2
        pointers['TX'] = currentPosition + 8 * 3
        # (Next Data group pointer, first channel group block pointer,
        # data block pointer, comment block pointer,
        # no recordID, reserved)
        dataBytes = (0, 0, 0, 0, 0, b'\0' * 7)
        fid.write(pack('<4QB7s', *dataBytes))
        return pointers


class CGBlock(dict):

    """ reads Channel Group block and saves in class dict
    """

    def __init__(self, fid=None, pointer=None):
        if fid is not None:
            self.read(fid, pointer)

    def read(self, fid, pointer):
        fid.seek(pointer)
        self['pointer'] = pointer
        (self['id'],
         self['reserved'],
         self['length'],
         self['link_count'],
         self['cg_cg_next'],
         self['cg_cn_first'],
         self['cg_tx_acq_name'],
         self['cg_si_acq_source'],
         self['cg_sr_first'],
         self['cg_md_comment'],
         self['cg_record_id'],
         self['cg_cycle_count'],
         self['cg_flags'],
         self['cg_path_separator'],
         self['cg_reserved'],
         self['cg_data_bytes'],
         self['cg_invalid_bytes']) = CGStruct.unpack(fid.read(104))
        if self['cg_md_comment']:  # comments exist
            self['Comment'] = {}
            self['Comment'].update(CommentBlock(fid, self['cg_md_comment']))
        if self['cg_tx_acq_name']:  # comments exist
            self['acq_name'] = {}
            self['acq_name'].update(CommentBlock(fid, self['cg_tx_acq_name']))

    def write(self, fid, cg_cycle_count, cg_data_bytes):
        pointers = {}
        # write block header
        currentPosition = _writeHeader(fid, b'##CG', 104, 6)
        # link section
        pointers['block_start'] = currentPosition
        currentPosition += 24
        pointers['CG'] = currentPosition
        pointers['CN'] = currentPosition + 8
        # pointers['ACN'] = currentPosition + 16
        # pointers['ACS'] = currentPosition + 24
        # pointers['SR'] = currentPosition + 32
        pointers['TX'] = currentPosition + 40
        pointers['cg_data_bytes'] = currentPosition + 72
        # (Next channel group pointer, first channel block pointer,
        # acquisition name pointer, acquisition source pointer,
        # sample reduction pointer, comment pointer,
        # no recordID, cycle count, no flags,
        # no character specified for path separator,
        # reserved, number of bytes taken by data in record,
        # no invalid bytes)
        dataBytes = (0, 0, 0, 0, 0, 0, 0,
                     cg_cycle_count, 0, 0,
                     b'\0' * 4,
                     cg_data_bytes, 0)
        fid.write(pack('<8Q2H4s2I', *dataBytes))
        return pointers


class CNBlock(dict):

    """ reads Channel block and saves in class dict
    """

    def read(self, fid, pointer):
        if pointer != 0 and pointer is not None:
            fid.seek(pointer)
            (self['id'],
             reserved,
             self['length'],
             self['link_count']) = HeaderStruct.unpack(fid.read(24))
            self['pointer'] = pointer
            # data section
            fid.seek(pointer + self['length'] - 72)
            (self['cn_type'],
             self['cn_sync_type'],
             self['cn_data_type'],
             self['cn_bit_offset'],
             self['cn_byte_offset'],
             self['cn_bit_count'],
             self['cn_flags'],
             self['cn_invalid_bit_pos'],
             self['cn_precision'],
             self['cn_reserved'],
             self['cn_attachment_count'],
             self['cn_val_range_min'],
             self['cn_val_range_max'],
             self['cn_limit_min'],
             self['cn_limit_max'],
             self['cn_limit_ext_min'],
             self['cn_limit_ext_max']) = CNStruct.unpack(fid.read(72))
            # Channel Group block : Links
            fid.seek(pointer + 24)
            (self['cn_cn_next'],
             self['cn_composition'],
             self['cn_tx_name'],
             self['cn_si_source'],
             self['cn_cc_conversion'],
             self['cn_data'],
             self['cn_md_unit'],
             self['cn_md_comment']) = unpack('<8Q', fid.read(64))
            if self['cn_attachment_count'] > 0:
                self['cn_at_reference'] = \
                    _mdfblockread(fid, LINK, self['cn_attachment_count'])
                self['attachment'] = {}
                if self['cn_attachment_count'] > 1:
                    for at in range(self['cn_attachment_count']):
                        self['attachment'][at] = \
                            ATBlock(fid, self['cn_at_reference'][at])
                else:
                    self['attachment'][0] = \
                        ATBlock(fid, self['cn_at_reference'])
            if self['link_count'] > (8 + self['cn_attachment_count']):
                self['cn_default_x'] = _mdfblockread(fid, LINK, 3)
            else:
                self['cn_default_x'] = None
            if self['cn_md_comment']:  # comments exist
                self['Comment'] = {}
                self['Comment'].update(CommentBlock(fid, self['cn_md_comment'], MDType='CN'))
            if self['cn_md_unit']:  # comments exist
                self['unit'] = {}
                self['unit'].update(CommentBlock(fid, self['cn_md_unit'], 'unit'))
            if self['cn_tx_name']:  # comments exist
                self['name'] = CommentBlock(fid, self['cn_tx_name'], MDType='CN')['name']
            self['name'] = self['name'].replace(':', '')

    def write(self, fid):
        pointers = {}
        # write block header
        # no attachement and default X
        currentPosition = _writeHeader(fid, b'##CN', 160, 8)
        # link section
        pointers['block_start'] = currentPosition
        currentPosition += 24
        pointers['CN'] = currentPosition
        pointers['CN_Comp'] = currentPosition + 8
        pointers['TX'] = currentPosition + 16
        pointers['SI'] = currentPosition + 24
        pointers['CC'] = currentPosition + 32
        pointers['data'] = currentPosition + 40
        pointers['Unit'] = currentPosition + 48
        pointers['Comment'] = currentPosition + 56
        # (Next channel block pointer, composition of channel pointer,
        # TXBlock pointer for channel name, source SIBlock pointer,
        # Conversion Channel CCBlock pointer, channel data pointer,
        # channel unit comment block pointer, channel comment block pointer,
        # no attachments and default_x
        # channel type, 0 normal, 2 master
        # sync type, 0 None, 1 time, 2 angle, 3 distance, 4 index
        # data type, bit offset, byte offset, bit count, no flags,
        # precision, reserved, attachments count,
        # val range min, val range max, val limit min, val limit max,
        # val limit ext min, val limit ext max)
        dataBytes = (0, 0, 0, 0, 0, 0, 0, 0,
                     self['cn_type'], self['cn_sync_type'],
                     self['cn_data_type'], self['cn_bit_offset'],
                     self['cn_byte_offset'], self['cn_bit_count'],
                     self['cn_flags'], 0, 0, 0, 0,
                     self['cn_val_range_min'], self['cn_val_range_max'],
                     0, 0, 0, 0)
        fid.write(pack('<8Q4B4I2BH6d', *dataBytes))
        return pointers


class CCBlock(dict):

    """ reads Channel Conversion block and saves in class dict
    """

    def __init__(self, fid=None, pointer=None):
        if fid is not None:
            self.read(fid, pointer)

    def read(self, fid, pointer):
        # block header
        if pointer != 0 and pointer is not None:
            fid.seek(pointer)
            self['pointer'] = pointer
            (self['id'],
             reserved,
             self['length'],
             self['link_count'],
             self['cc_tx_name'],
             self['cc_md_unit'],
             self['cc_md_comment'],
             self['cc_cc_inverse']) = CCStruct1.unpack(fid.read(56))
            if self['link_count'] - 4 > 0:  # can be no links for cc_ref
                self['cc_ref'] = _mdfblockread(fid, LINK,
                                               self['link_count'] - 4)
            # data section
            (self['cc_type'],
             self['cc_precision'],
             self['cc_flags'],
             self['cc_ref_count'],
             self['cc_val_count'],
             self['cc_phy_range_min'],
             self['cc_phy_range_max']) = CCStruct2.unpack(fid.read(24))
            if self['cc_val_count']:
                self['cc_val'] = _mdfblockread(fid, REAL, self['cc_val_count'])
            if self['cc_type'] == 3:  # reads Algebraic formula
                pointer = self['cc_ref']
                self['cc_ref'] = {}
                self['cc_ref'].update(CommentBlock(fid, pointer))
            elif self['cc_type']in (7, 8, 9, 10):  # text list
                self['cc_ref'] = list(self['cc_ref'])
                for i in range(self['cc_ref_count']):
                    fid.seek(self['cc_ref'][i])
                    # find if TX/MD or another CCBlock
                    ID = unpack('4s', fid.read(4))[0]
                    # for algebraic formulae
                    if ID in ('##TX', '##MD', b'##TX', b'##MD'):
                        temp = CommentBlock(fid, self['cc_ref'][i])
                        self['cc_ref'][i] = temp['Comment']
                    elif ID in ('##CC', b'##CC'):  # for table conversion
                        # much more complicated nesting conversions !!!
                        self['cc_ref'][i] = CCBlock(fid, self['cc_ref'][i])
            if self['cc_md_comment']:  # comments exist
                self['Comment'] = {}
                self['Comment'].update(CommentBlock(fid, self['cc_md_comment'],
                                               MDType='CC'))
            if self['cc_md_unit']:  # comments exist
                self['unit'] = {}
                self['unit'].update(CommentBlock(fid, self['cc_md_unit']))
            if self['cc_tx_name']:  # comments exist
                self['name'] = {}
                self['name'].update(CommentBlock(fid, self['cc_tx_name']))
        else:  # no conversion
            self['cc_type'] = 0


class CABlock(dict):

    """ reads Channel Array block and saves in class dict
    """

    def __init__(self, fid, pointer):
        # block header
        if pointer != 0 and pointer is not None:
            fid.seek(pointer)
            (self['id'],
             reserved,
             self['length'],
             self['link_count']) = HeaderStruct.unpack(fid.read(24))
            self['pointer'] = pointer
            # reads data section
            fid.seek(pointer + 24 + self['link_count'] * 8)
            (self['ca_type'],
             self['ca_storage'],
             self['ca_ndim'],
             self['ca_flags'],
             self['ca_byte_offset_base'],
             self['ca_invalid_bit_pos_base']) = unpack('2BHIiI', fid.read(16))
            self['ca_dim_size'] = _mdfblockread(fid, UINT64, self['ca_ndim'])
            try:  # more than one dimension, processing dict
                self['SNd'] = 0
                self['PNd'] = 1
                for x in self['ca_dim_size']:
                    self['SNd'] += x
                    self['PNd'] *= x
            except:  # only one dimension, processing int
                self['SNd'] = self['ca_dim_size']
                self['PNd'] = self['SNd']
            if 1 << 5 & self['ca_flags']:  # bit5
                self['ca_axis_value'] = \
                    _mdfblockread(fid, REAL, self['SNd'])
            if self['ca_storage'] >= 1:
                self['ca_cycle_count'] = \
                    _mdfblockread(fid, UINT64, self['PNd'])
            # Channel Conversion block : Links
            fid.seek(pointer + 24)
            # point to CN for array of structures or CA for array of array
            self['ca_composition'] = _mdfblockread(fid, LINK, 1)
            if self['ca_storage'] == 2:
                self['ca_data'] = _mdfblockread(fid, LINK, self['PNd'])
            if 1 << 0 & self['ca_flags']:  # bit 0
                self['ca_dynamic_size'] = \
                    _mdfblockread(fid, LINK, self['ca_ndim'] * 3)
            if 1 << 1 & self['ca_flags']:  # bit 1
                self['ca_input_quantity'] = \
                    _mdfblockread(fid, LINK, self['ca_ndim'] * 3)
            if 1 << 2 & self['ca_flags']:  # bit 2
                self['ca_output_quantity'] = \
                    _mdfblockread(fid, LINK, 3)
            if 1 << 3 & self['ca_flags']:  # bit 3
                self['ca_comparison_quantity'] = _mdfblockread(fid, LINK, 3)
            if 1 << 4 & self['ca_flags']:  # bit 4
                self['ca_cc_axis_conversion'] = \
                    _mdfblockread(fid, LINK, self['ca_ndim'])
            if 1 << 4 & self['ca_flags'] and not 1 << 5 & self['ca_flags']:
                # bit 4 and 5
                self['ca_axis'] = _mdfblockread(fid, LINK, self['ca_ndim'] * 3)
            # nested arrays
            if self['ca_composition']:
                self['CABlock'] = CABlock(fid, self['ca_composition'])


class ATBlock(dict):

    """ reads Attachment block and saves in class dict
    """

    def __init__(self, fid, pointer):
        # block header
        if pointer != 0 and pointer is not None:
            fid.seek(pointer)
            (self['id'],
             reserved,
             self['length'],
             self['link_count'],
             self['at_at_next'],
             self['at_tx_filename'],
             self['at_tx_mimetype'],
             self['at_md_comment'],
             self['at_flags'],
             self['at_creator_index'],
             at_reserved,
             self['at_md5_checksum'],
             self['at_original_size'],
             self['at_embedded_size']) = unpack('<4sI6Q2HI16s2Q', fid.read(96))
            if self['at_embedded_size'] > 0:
                self['at_embedded_data'] = fid.read(self['at_embedded_size'])
            if self['at_md_comment']:  # comments exist
                self['Comment'] = {}
                self['Comment'].update(CommentBlock(fid, self['at_md_comment']))


class EVBlock(dict):

    """ reads Event block and saves in class dict
    """

    def __init__(self, fid, pointer):
        # block header
        if pointer != 0 and pointer is not None:
            fid.seek(pointer)
            (self['id'],
             reserved,
             self['length'],
             self['link_count']) = HeaderStruct.unpack(fid.read(24))
            # data section
            fid.seek(pointer + self['length'] - 32)
            (self['ev_type'],
             self['ev_sync_type'],
             self['ev_range_type'],
             self['ev_cause'],
             self['ev_flags'],
             ev_reserved,
             self['ev_scope_count'],
             self['ev_attachment_count'],
             self['ev_creator_index'],
             self['ev_sync_base_value'],
             self['ev_sync_factor']) = unpack('<5B3sI2Hqd', fid.read(32))
            # link section
            fid.seek(pointer + 24)
            (self['ev_ev_next'],
             self['ev_ev_parent'],
             self['ev_ev_range'],
             self['ev_tx_name'],
             self['ev_md_comment']) = unpack('<5Q', fid.read(40))
            self['ev_scope'] = _mdfblockread(fid, LINK, self['ev_scope_count'])
            # post treatment
            if self['ev_cause'] == 0:
                self['ev_cause'] = 'OTHER'
            elif self['ev_cause'] == 1:
                self['ev_cause'] == 'ERROR'
            elif self['ev_cause'] == 2:
                self['ev_cause'] == 'TOOL'
            elif self['ev_cause'] == 3:
                self['ev_cause'] == 'SCRIPT'
            elif self['ev_cause'] == 4:
                self['ev_cause'] == 'USER'
            if self['ev_attachment_count'] > 0:
                self['ev_at_reference'] = \
                    _mdfblockread(fid, LINK, self['ev_attachment_count'])
                for at in range(self['ev_attachment_count']):
                    self['attachment'][at] = \
                        ATBlock(fid, self['ev_at_reference'][at])
            if self['ev_md_comment']:  # comments exist
                self['Comment'] = {}
                self['Comment'].update(CommentBlock(fid, self['ev_md_comment']))
            if self['ev_tx_name']:  # comments exist
                self['name'] = {}
                self['name'].update(CommentBlock(fid, self['ev_tx_name']))


class SRBlock(dict):

    """ reads Sample Reduction block and saves in class dict
    """

    def __init__(self, fid, pointer):
        # block header
        if pointer != 0 and pointer is not None:
            fid.seek(pointer)
            (self['id'],
             reserved,
             self['length'],
             self['link_count'],
             self['sr_sr_next'],
             self['sr_data'],
             self['sr_cycle_count'],
             self['sr_interval'],
             self['sr_sync_type'],
             self['sr_flags'],
             sr_reserved) = SRStruct.unpack(fid.read(64))


class SIBlock(dict):

    """ reads Source Information block and saves in class dict
    """

    def read(self, fid, pointer):
        # block header
        if pointer != 0 and pointer is not None:
            fid.seek(pointer)
            (self['id'],
             reserved,
             self['length'],
             self['link_count'],
             self['si_tx_name'],
             self['si_tx_path'],
             self['si_md_comment'],
             self['si_type'],
             self['si_bus_type'],
             self['si_flags'],
             si_reserved) = SIStruct.unpack(fid.read(56))
            if self['si_type'] == 0:
                self['si_type'] = 'OTHER'  # unknown
            elif self['si_type'] == 1:
                self['si_type'] = 'ECU'
            elif self['si_type'] == 2:
                self['si_type'] = 'BUS'
            elif self['si_type'] == 3:
                self['si_type'] = 'I/O'
            elif self['si_type'] == 4:
                self['si_type'] = 'TOOL'
            elif self['si_type'] == 5:
                self['si_type'] = 'USER'
            if self['si_bus_type'] == 0:
                self['si_bus_type'] = 'NONE'
            elif self['si_bus_type'] == 1:
                self['si_bus_type'] = 'OTHER'
            elif self['si_bus_type'] == 2:
                self['si_bus_type'] = 'CAN'
            elif self['si_bus_type'] == 3:
                self['si_bus_type'] = 'LIN'
            elif self['si_bus_type'] == 4:
                self['si_bus_type'] = 'MOST'
            elif self['si_bus_type'] == 5:
                self['si_bus_type'] = 'FLEXRAY'
            elif self['si_bus_type'] == 6:
                self['si_bus_type'] = 'K_LINE'
            elif self['si_bus_type'] == 7:
                self['si_bus_type'] = 'ETHERNET'
            elif self['si_bus_type'] == 8:
                self['si_bus_type'] = 'USB'
            # post treatment
            self['source_name'] = {}
            self['source_name'].update(CommentBlock(fid, self['si_tx_name']))
            self['source_path'] = {}
            self['source_path'].update(CommentBlock(fid, self['si_tx_path']))
            self['comment'] = {}
            self['comment'].update(CommentBlock(fid, self['si_md_comment'],
                                                MDType='SI'))


class DLBlock(dict):

    """ reads Data List block
    """

    def __init__(self, fid, link_count):
        # block header is already read
        self['dl_dl_next'] = unpack('<Q', fid.read(8))[0]
        self['dl_data'] = {}
        self['dl_data'][0] = unpack('<{}Q'.format(link_count - 1),
                                    fid.read(8 * (link_count - 1)))
        (self['dl_flags'],
         dl_reserved,
         self['dl_count']) = unpack('<B3sI', fid.read(8))
        if self['dl_flags']:  # equal length datalist
            self['dl_equal_length'] = unpack('<Q', fid.read(8))[0]
        else:  # datalist defined by byte offset
            self['dl_offset'] = unpack('<{}Q'.format(self['dl_count']),
                                       fid.read(8 * self['dl_count']))


class DZBlock(dict):

    """ reads Data List block
    """

    def __init__(self, fid):
        # block header is already read
        (self['dz_org_block_type'],
         self['dz_zip_type'],
         dz_reserved,
         self['dz_zip_parameter'],
         self['dz_org_data_length'],
         self['dz_data_length']) = DZStruct.unpack(fid.read(24))


class HLBlock(dict):

    """ reads Header List block
    """

    def __init__(self, fid):
        (self['hl_dl_first'],
         self['hl_flags'],
         self['hl_zip_type'],
         hl_reserved) = HLStruct.unpack(fid.read(16))


class info4(dict):
    __slots__ = ['fileName', 'fid', 'zipfile']
    """ information block parser fo MDF file version 4.x

    Attributes
    --------------
    fileName : str
        name of file
    fid
        file identifier
    zipfile
        flag to indicate the mdf4 is packaged in a zip

    Notes
    --------
    mdfinfo(FILENAME) contains a dict of structures, for
    each data group, containing key information about all channels in each
    group. FILENAME is a string that specifies the name of the MDF file.
    Either file name or fid should be given.
    General dictionary structure is the following

    - mdfinfo['HD'] header block
    - mdfinfo['DG'][dataGroup] Data Group block
    - mdfinfo['CG'][dataGroup][channelGroup] Channel Group block
    Channel block including text blocks for comment and identifier
    - mdfinfo['CN'][dataGroup][channelGroup][channel]
    Channel conversion information
    - mdfinfo['CC'][dataGroup][channelGroup][channel]
    """

    def __init__(self, fileName=None, fid=None, minimal=0):
        """ info4 class constructor

        Parameters
        ----------------
        fileName : str
            file name
        fid : float
            file identifier
        minimal : int
            0 will load every metadata
            1 will load DG, CG, CN and CC (for noDataLoading)
            2 will load only DG (for normal reading)

        Notes
        ---------
        Either fileName or fid can be used as argument
        """
        self['ID'] = {}  # Identifier Block
        self['HD'] = {}  # Header Block
        self['FH'] = {}
        self['CH'] = {}
        self['DG'] = {}  # Data Group Block
        self['CG'] = {}  # Channel Group Block
        self['CN'] = {}  # Channel Block
        self['CC'] = {}  # Conversion block
        self['AT'] = {}  # Attachment block
        self['allChannelList'] = set()  # all channels
        self.fileName = fileName
        self.fid = None
        if fid is None and fileName is not None:
            # Open file
            (self.fid, self.fileName, self.zipfile) = _open_MDF(self.fileName)
        if self.fileName is not None and fid is None:
            self.readinfo(self.fid, minimal)
            # Close the file
            self.fid.close()
            if self.zipfile:  # temporary uncompressed file, to be removed
                remove(fileName)
        elif self.fileName is None and fid is not None:
            # called by mdfreader.mdfinfo
            self.readinfo(fid, minimal)

    def readinfo(self, fid, minimal):
        """ read all file blocks except data

        Parameters
        ----------------
        fid : float
            file identifier
        minimal: falg
            to activate minimum content reading for raw data fetching
        """
        # reads IDBlock
        self['ID'].update(IDBlock(fid))

        # reads Header HDBlock
        self['HD'].update(HDBlock(fid))

        if not minimal:
            # print('reads File History blocks, always exists', file=stderr)
            fh = 0  # index of fh blocks
            self['FH'][fh] = {}
            self['FH'][fh] .update(FHBlock(fid, self['HD']['hd_fh_first']))
            while self['FH'][fh]['fh_fh_next']:
                self['FH'][fh + 1] = {}
                self['FH'][fh + 1].update(FHBlock(fid, self['FH'][fh]['fh_fh_next']))
                fh += 1

        # print('reads Channel Hierarchy blocks', file=stderr)
        if self['HD']['hd_ch_first']:
            ch = 0
            self['CH'][ch] = {}
            self['CH'][ch].update(CHBlock(fid, self['HD']['hd_ch_first']))
            while self['CH'][ch]['ch_ch_next']:
                self['CH'][ch].update(CHBlock(fid, self['CH'][ch]['ch_ch_next']))
                ch += 1

        # reads Attachment block
        if not minimal:
            self['AT'] = self.readATBlock(fid, self['HD']['hd_at_first'])

        # reads Event Block
        if self['HD']['hd_ev_first'] and not minimal:
            ev = 0
            self['EVBlock'] = {}
            self['EVBlock'][ev] = EVBlock(fid, self['HD']['hd_ev_first'])
            while self['EVBlock'][ev]['ev_ev_next']:
                ev += 1
                self['EVBlock'][ev] = EVBlock(fid, self['EVBlock'][ev - 1]['ev_ev_next'])

        # reads Data Group Blocks and recursively the other related blocks
        self.readDGBlock(fid, None, minimal)

    def readDGBlock(self, fid, channelNameList=False, minimal=0):
        """reads Data Group Blocks

        Parameters
        ----------------
        fid : float
            file identifier
        channelNameList : bool
            Flag to reads only channel blocks for listChannels4 method
        minimal: falg
            to activate minimum content reading for raw data fetching
        """
        self['ChannelNamesByDG'] = {}
        if self['HD']['hd_dg_first']:
            dg = 0
            self['DG'][dg] = {}
            self['DG'][dg].update(DGBlock(fid, self['HD']['hd_dg_first']))
            self['ChannelNamesByDG'][dg] = set()
            if minimal < 2:
                # reads Channel Group blocks
                self.readCGBlock(fid, dg, channelNameList, minimal)
            while self['DG'][dg]['dg_dg_next']:
                dg += 1
                self['DG'][dg] = {}
                self['DG'][dg].update(DGBlock(fid, self['DG'][dg - 1]['dg_dg_next']))
                self['ChannelNamesByDG'][dg] = set()
                if minimal < 2:
                    # reads Channel Group blocks
                    self.readCGBlock(fid, dg, channelNameList, minimal)

    def readCGBlock(self, fid, dg, channelNameList=False, minimal=0):
        """reads Channel Group blocks

        Parameters
        ----------------
        fid : float
            file identifier
        dg : int
            data group number
        channelNameList : bool
            Flag to reads only channel blocks for listChannels4 method
        minimal: falg
            to activate minimum content reading for raw data fetching
        """
        if self['DG'][dg]['dg_cg_first']:
            cg = 0
            self['CN'][dg] = {}
            self['CN'][dg][cg] = {}
            self['CC'][dg] = {}
            self['CC'][dg][cg] = {}
            self['CG'][dg] = {}
            self['CG'][dg][cg] = {}
            self['CG'][dg][cg].update(CGBlock(fid, self['DG'][dg]['dg_cg_first']))
            VLSDCGBlock = []
            vlsd = False  # flag for at least one VLSD channel in CG

            if not channelNameList and minimal:
                # reads Source Information Block
                temp = SIBlock()
                temp.read(fid, self['CG'][dg][cg]['cg_si_acq_source'])
                if temp is not None:
                    self['CG'][dg][cg]['SI'] = dict()
                    self['CG'][dg][cg]['SI'].update(temp)

                # reads Sample Reduction Block
                self['CG'][dg][cg]['SR'] = self.readSRBlock(fid, self['CG'][dg][cg]['cg_sr_first'])

            if not self['CG'][dg][cg]['cg_flags'] & 0b1:  # if not a VLSD channel group
                # reads Channel Block
                vlsd = self.readCNBlock(fid, dg, cg, channelNameList, minimal)
                if vlsd:
                    # VLSD needs to rename and append records but with python 2.x impossible,
                    # convert name to compatible python identifier
                    for cn in self['CN'][dg][cg]:
                        self['CN'][dg][cg][cn]['name'] = _convertName(
                            self['CN'][dg][cg][cn]['name'].encode('latin-1', ' ignore'))
            else:
                VLSDCGBlock.append(cg)

            while self['CG'][dg][cg]['cg_cg_next']:
                cg += 1
                vlsd = False  # flag for at least one VLSD channel in CG
                self['CG'][dg][cg] = {}
                self['CG'][dg][cg].update(CGBlock(fid, self['CG'][dg][cg - 1]['cg_cg_next']))
                self['CN'][dg][cg] = {}
                self['CC'][dg][cg] = {}
                if not channelNameList and not minimal:
                    # reads Source Information Block
                    temp = SIBlock()
                    temp.read(fid, self['CG'][dg][cg]['cg_si_acq_source'])
                    if temp is not None:
                        self['CG'][dg][cg]['SI'] = dict()
                        self['CG'][dg][cg]['SI'].update(temp)

                    # reads Sample Reduction Block
                    self['CG'][dg][cg]['SR'] = self.readSRBlock(fid, self['CG'][dg][cg]['cg_sr_first'])

                if not self['CG'][dg][cg]['cg_flags'] & 0b1:  # if not a VLSD channel group
                    # reads Channel Block
                    vlsd = self.readCNBlock(fid, dg, cg, channelNameList, minimal)
                    if vlsd:
                        # VLSD needs to rename and append records but with python 2.x impossible,
                        # convert name to compatible python identifier
                        for cn in self['CN'][dg][cg]:
                            self['CN'][dg][cg][cn]['name'] = _convertName(
                                self['CN'][dg][cg][cn]['name'].encode('latin-1', ' ignore'))
                else:
                    VLSDCGBlock.append(cg)

            if VLSDCGBlock:  # VLSD CG Block exiting
                self['VLSD_CG'] = {}

            # Matching VLSD CGBlock with corresponding channel
            for VLSDcg in VLSDCGBlock:
                VLSDCGBlockAdress = self['CG'][dg][VLSDcg]['pointer']
                for cg in self['CG'][dg]:
                    if cg not in VLSDCGBlock:
                        for cn in self['CN'][dg][cg]:
                            if VLSDCGBlockAdress == self['CN'][dg][cg][cn]['cn_data']:
                                # found matching channel with VLSD CGBlock
                                temp = {}
                                temp['cg_cn'] = (cg, cn)
                                self['VLSD_CG'][self['CG'][dg][VLSDcg]['cg_record_id']] = temp
                                break

            # reorder channel blocks and related blocks(CC, SI, AT, CA) based on byte offset
            # this reorder is meant to improve performance while parsing records using core.records.fromfile
            # as it will not use cn_byte_offset
            # first, calculate new mapping/order
            nChannel = len(self['CN'][dg][cg])
            Map = zeros(shape=nChannel, dtype=[('index', 'u4'), ('bit_offset', 'u4')])
            for cn in range(nChannel):
                Map[cn] = (cn, self['CN'][dg][cg][cn]['cn_byte_offset'] * 8 + self['CN'][dg][cg][cn]['cn_bit_offset'])
            orderedMap = sort(Map, order='bit_offset')

            toChangeIndex = Map == orderedMap
            for cn in range(nChannel):
                if not toChangeIndex[cn]:
                    # offset all indexes of indexes to be moved
                    self['CN'][dg][cg][cn + nChannel] = self['CN'][dg][cg].pop(cn)
                    self['CC'][dg][cg][cn + nChannel] = self['CC'][dg][cg].pop(cn)
            for cn in range(nChannel):
                if not toChangeIndex[cn]:
                    # change to ordered index
                    self['CN'][dg][cg][cn] = self['CN'][dg][cg].pop(orderedMap[cn][0] + nChannel)
                    self['CC'][dg][cg][cn] = self['CC'][dg][cg].pop(orderedMap[cn][0] + nChannel)

    def readCNBlock(self, fid, dg, cg, channelNameList=False, minimal=0):
        """reads Channel blocks

        Parameters
        ----------------
        fid : float
            file identifier
        dg : int
            data group number
        cg : int
            channel group number in data group
        channelNameList : bool
            Flag to reads only channel blocks for listChannels4 method
        minimal: falg
            to activate minimum content reading for raw data fetching
        """
        cn = 0
        vlsd = False
        self['CN'][dg][cg][cn] = {}
        self['CC'][dg][cg][cn] = {}
        self['CN'][dg][cg][cn] = {}
        temp = CNBlock()
        temp.read(fid, self['CG'][dg][cg]['cg_cn_first'])
        if temp is not None:
            self['CN'][dg][cg][cn].update(temp)
        MLSDChannels = []
        # check for MLSD
        if self['CN'][dg][cg][cn]['cn_type'] == 5:
            MLSDChannels.append(cn)
        # check if already existing channel name
        if self['CN'][dg][cg][cn]['name'] in self['ChannelNamesByDG'][dg]:  # for unsorted data
            self['CN'][dg][cg][cn]['name'] = u'{0}_{1}_{2}_{3}'.format(self['CN'][dg][cg][cn]['name'], dg, cg, cn)
        elif self['CN'][dg][cg][cn]['name'] in self['allChannelList']:
            # doublon name or master channel
            self['CN'][dg][cg][cn]['name'] = u'{0}_{1}'.format(self['CN'][dg][cg][cn]['name'], dg)
        if self['CN'][dg][cg][cn]['cn_type'] == 1 and PythonVersion < 3:
            # VLSD needs to rename and append records but with python 2.x impossible,
            # convert name to compatible python identifier
            vlsd = True
        self['ChannelNamesByDG'][dg].add(self['CN'][dg][cg][cn]['name'])
        self['allChannelList'].add(self['CN'][dg][cg][cn]['name'])

        if self['CG'][dg][cg]['cg_cn_first']:  # Can be NIL for VLSD
            # reads Channel Conversion Block
            self['CC'][dg][cg][cn] = {}
            temp = CCBlock()
            temp.read(fid, self['CN'][dg][cg][cn]['cn_cc_conversion'])
            self['CC'][dg][cg][cn].update(temp)
            if not channelNameList:
                if not minimal:
                    # reads Channel Source Information
                    temp = SIBlock()
                    temp.read(fid, self['CN'][dg][cg][cn]['cn_si_source'])
                    if temp is not None:
                        self['CN'][dg][cg][cn]['SI'] = dict()
                        self['CN'][dg][cg][cn]['SI'].update(temp)

                # reads Channel Array Block
                if self['CN'][dg][cg][cn]['cn_composition']:
                    # composition but can be either structure of channels or array
                    fid.seek(self['CN'][dg][cg][cn]['cn_composition'])
                    id = fid.read(4)
                    if id in ('##CA', b'##CA'):
                        self['CN'][dg][cg][cn]['CABlock'] = \
                            CABlock(fid, self['CN'][dg][cg][cn]['cn_composition'])
                    elif id in ('##CN', b'##CN'):
                        self['CN'][dg][cg][cn]['CN'] = {}
                        temp = CNBlock()
                        temp.read(fid, self['CN'][dg][cg][cn]['cn_composition'])
                        self['CN'][dg][cg][cn]['CN'].update(temp)
                    else:
                        raise('unknown channel composition')

                if not minimal:
                    # reads Attachment Block
                    if self['CN'][dg][cg][cn]['cn_attachment_count'] > 1:
                        for at in range(self['CN'][dg][cg][cn]['cn_attachment_count']):
                            self['CN'][dg][cg][cn]['attachment'][at].update(self.readATBlock(fid, self['CN'][dg][cg][cn]['cn_at_reference'][at]))
                    elif self['CN'][dg][cg][cn]['cn_attachment_count'] == 1:
                        self['CN'][dg][cg][cn]['attachment'][0].update(self.readATBlock(fid, self['CN'][dg][cg][cn]['cn_at_reference']))

            while self['CN'][dg][cg][cn]['cn_cn_next']:
                cn = cn + 1
                self['CN'][dg][cg][cn] = {}
                temp = CNBlock()
                temp.read(fid, self['CN'][dg][cg][cn - 1]['cn_cn_next'])
                self['CN'][dg][cg][cn].update(temp)
                # check for MLSD
                if self['CN'][dg][cg][cn]['cn_type'] == 5:
                    MLSDChannels.append(cn)
                # reads Channel Conversion Block
                self['CC'][dg][cg][cn] = {}
                temp = CCBlock(fid, self['CN'][dg][cg][cn]['cn_cc_conversion'])
                self['CC'][dg][cg][cn].update(temp)
                if not channelNameList:
                    if not minimal:
                        # reads Channel Source Information
                        temp = SIBlock()
                        temp.read(fid, self['CN'][dg][cg][cn]['cn_si_source'])
                        if temp is not None:
                            self['CN'][dg][cg][cn]['SI'] = dict()
                            self['CN'][dg][cg][cn]['SI'].update(temp)

                    # check if already existing channel name
                    if self['CN'][dg][cg][cn]['name'] in self['ChannelNamesByDG'][dg]:  # for unsorted data
                        self['CN'][dg][cg][cn]['name'] = \
                            u'{0}_{1}_{2}_{3}'.format(self['CN'][dg][cg][cn]['name'], dg, cg, cn)
                    elif self['CN'][dg][cg][cn]['name'] in self['allChannelList']:
                        # doublon name or master channel
                        self['CN'][dg][cg][cn]['name'] = u'{0}_{1}'.format(self['CN'][dg][cg][cn]['name'], dg)
                    if self['CN'][dg][cg][cn]['cn_type'] == 1 and PythonVersion < 3:
                        # VLSD needs to rename and append records but with python 2.x impossible,
                        # convert name to compatible python identifier
                        vlsd = True
                    self['ChannelNamesByDG'][dg].add(self['CN'][dg][cg][cn]['name'])
                    self['allChannelList'].add(self['CN'][dg][cg][cn]['name'])

                    # reads Channel Array Block
                    if self['CN'][dg][cg][cn]['cn_composition']:
                        # composition but can be either structure of channels or array
                        fid.seek(self['CN'][dg][cg][cn]['cn_composition'])
                        id = fid.read(4)
                        if id in ('##CA', b'##CA'):
                            self['CN'][dg][cg][cn]['CABlock'] = \
                                CABlock(fid, self['CN'][dg][cg][cn]['cn_composition'])
                        elif id in ('##CN', b'##CN'):
                            self['CN'][dg][cg][cn]['CN'] = {}
                            temp = CNBlock()
                            temp.read(fid, self['CN'][dg][cg][cn]['cn_composition'])
                            self['CN'][dg][cg][cn]['CN'].update(temp)
                        else:
                            raise('unknown channel composition')

                    if not minimal:
                        # reads Attachment Block
                        if self['CN'][dg][cg][cn]['cn_attachment_count'] > 1:
                            for at in range(self['CN'][dg][cg][cn]['cn_attachment_count']):
                                print(self['CN'][dg][cg][cn]['cn_at_reference'][at], file=stderr)
                                self['CN'][dg][cg][cn]['attachment'][at].update(self.readATBlock(fid, self['CN'][dg][cg][cn]['cn_at_reference'][at]))
                        elif self['CN'][dg][cg][cn]['cn_attachment_count'] == 1:
                            self['CN'][dg][cg][cn]['attachment'][0].update(self.readATBlock(fid, self['CN'][dg][cg][cn]['cn_at_reference']))

        MLSDChannels = self.readComposition(fid, dg, cg, MLSDChannels, channelNameList=False)

        if MLSDChannels:
            self['MLSD'] = {}
            self['MLSD'][dg] = {}
            self['MLSD'][dg][cg] = {}
        for MLSDcn in MLSDChannels:
            for cn in self['CN'][dg][cg]:
                if self['CN'][dg][cg][cn]['pointer'] == self['CN'][dg][cg][MLSDcn]['cn_data']:
                    self['MLSD'][dg][cg][MLSDcn] = cn
                    break
        return vlsd

    def cleanDGinfo(self, dg):
        """ delete CN,CC and CG blocks related to data group

        Parameters
        ----------------
        dg : int
            data group number
        """
        try:
            self['CN'][dg] = {}
        except KeyError:
            pass
        try:
            self['CC'][dg] = {}
        except KeyError:
            pass
        try:
            self['CA'][dg] = {}
        except KeyError:
            pass
        try:
            self['CG'][dg] = {}
        except KeyError:
            pass

    def readComposition(self, fid, dg, cg, MLSDChannels,
                        channelNameList=False):
        """check for composition of channels, arrays or structures

        Parameters
        ----------------
        fid : float
            file identifier
        dg : int
            data group number
        cg : int
            channel group number in data group
        MLSDChannels : list of int
            channel numbers
        channelNameList : bool
            Flag to reads only channel blocks for listChannels4 method

        Returns
        -----------
        MLSDChannels list of appended Maximum Length Sampling Data channels
        """
        chan = max(self['CN'][dg][cg].keys()) + 1
        for cn in list(self['CN'][dg][cg].keys()):
            if self['CN'][dg][cg][cn]['cn_composition']:
                fid.seek(self['CN'][dg][cg][cn]['cn_composition'])
                ID = unpack('4s', fid.read(4))[0]
                if ID in ('##CN', b'##CN'):  # Structures
                    self['CN'][dg][cg][chan] = CNBlock()
                    self['CN'][dg][cg][chan].read(fid,
                        self['CN'][dg][cg][cn]['cn_composition'])
                    self['CC'][dg][cg][chan] = \
                        CCBlock(fid, self['CN']
                                [dg][cg][chan]['cn_cc_conversion'])
                    if self['CN'][dg][cg][chan]['cn_type'] == 5:
                        MLSDChannels.append(chan)
                    while self['CN'][dg][cg][chan]['cn_cn_next']:
                        chan += 1
                        self['CN'][dg][cg][chan] = CNBlock()
                        self['CN'][dg][cg][chan]\
                            .read(fid, self['CN']
                                  [dg][cg][chan - 1]['cn_cn_next'])
                        self['CC'][dg][cg][chan] = \
                            CCBlock(fid, self['CN']
                                    [dg][cg][chan]['cn_cc_conversion'])
                        if self['CN'][dg][cg][chan]['cn_type'] == 5:
                            MLSDChannels.append(chan)
                    # makes the channel virtual
                    self['CN'][dg][cg][cn]['cn_type'] = 6
                elif ID in ('##CA', b'##CA'):  # arrays
                    pass
                else:
                    print('unknown channel composition', file=stderr)
        return MLSDChannels

    def readSRBlock(self, fid, pointer):
        """reads Sample Reduction Blocks

        Parameters
        ----------------
        fid : float
            file identifier
        pointer : int
            position of SRBlock in file

        Returns
        -----------
        Sample Reduction Blocks in a dict
        """
        if pointer > 0:
            sr = 0
            srBlocks = {}
            srBlocks[0] = SRBlock(fid, pointer)
            while srBlocks[sr]['sr_sr_next'] > 0:
                sr += 1
                srBlocks[sr] = SRBlock(fid, srBlocks[sr - 1]['sr_sr_next'])
            return srBlocks

    def readATBlock(selfself, fid, pointer):
        """reads Attachment blocks

        Parameters
        ----------------
        fid : float
            file identifier
        pointer : int
            position of ATBlock in file

        Returns
        -----------
        Attachments Blocks in a dict
        """
        if pointer > 0:
            at = 0
            atBlocks = {}
            if type(pointer) in (tuple, list):
                pointer = pointer[0]
            atBlocks[0] = ATBlock(fid, pointer)
            while atBlocks[at]['at_at_next'] > 0:
                at += 1
                atBlocks[at] = (ATBlock(fid, atBlocks[at - 1]['at_at_next']))
            return atBlocks

    def listChannels4(self, fileName=None, fid=None):
        """ Read MDF file and extract its complete structure

        Parameters
        ----------------
        fileName : str
            file name

        Returns
        -----------
        list of channel names contained in file
        """
        if fileName is not None:
            self.fileName = fileName
        # Open file
        if fid is None and fileName is not None:
            # Open file
            (fid, fileName, zipfile) = _open_MDF(self.fileName)
        channelNameList = []
        # reads Header HDBlock
        self['HD'].update(HDBlock(fid))

        # reads Data Group, channel groups and channel Blocks
        # recursively but not the other metadata block
        self.readDGBlock(fid, True)

        for dg in self['DG']:
            for cg in self['CG'][dg]:
                for cn in self['CN'][dg][cg]:
                    channelNameList.append(self['CN'][dg][cg][cn]['name'])

        # CLose the file
        fid.close()
        return channelNameList


def _generateDummyMDF4(info, channelList):
    """ computes MasterChannelList and dummy mdf dict from an info object

    Parameters
    ----------------
    info : info object
        information structure of file

    channelList : list of str
        channel list

    Returns
    -----------
    a dict which keys are master channels in files with values a list of related channels of the raster
    """
    MasterChannelList = {}
    allChannelList = set()
    mdfdict = {}
    for dg in info['DG']:
        master = ''
        mastertype = 0
        for cg in info['CG'][dg]:
            channelNameList = []
            for cn in info['CN'][dg][cg]:
                if info['CN'][dg][cg][cn]['cn_data_type'] == 13:
                    for name in ('ms', 'minute', 'hour', 'day', 'month', 'year'):
                        channelNameList.append(name)
                        allChannelList.add(name)
                        mdfdict[name] = {}
                        mdfdict[name][dataField] = None
                        mdfdict[name][descriptionField] = name
                        mdfdict[name][unitField] = name
                elif info['CN'][dg][cg][cn]['cn_data_type'] == 14:
                    for name in ('ms', 'days'):
                        channelNameList.append(name)
                        allChannelList.add(name)
                        mdfdict[name] = {}
                        mdfdict[name][dataField] = None
                        mdfdict[name][descriptionField] = name
                        mdfdict[name][unitField] = name
                else:
                    name = info['CN'][dg][cg][cn]['name']
                    if name in allChannelList:
                        name = ''.join([name, '_{}'.format(dg)])
                    if channelList is None or name in channelList:
                        channelNameList.append(name)
                        allChannelList.add(name)
                        # create mdf channel
                        mdfdict[name] = {}
                        mdfdict[name][dataField] = None
                        if 'description' in info['CN'][dg][cg][cn]:
                            mdfdict[name][descriptionField] = info['CN'][dg][cg][cn]['description']
                        else:
                            mdfdict[name][descriptionField] = ''
                        if 'unit' in info['CN'][dg][cg][cn]:
                            mdfdict[name][unitField] = info['CN'][dg][cg][cn]['unit']
                        else:
                            mdfdict[name][unitField] = ''
                        mdfdict[name][masterField] = 0  # default is time
                        mdfdict[name][masterTypeField] = None
                    if info['CN'][dg][cg][cn]['cn_sync_type']:
                        # master channel of cg
                        master = name
                        mastertype = info['CN'][dg][cg][cn]['cn_sync_type']
            for chan in channelNameList:
                mdfdict[chan][masterField] = master
                mdfdict[chan][masterTypeField] = mastertype
        MasterChannelList[master] = channelNameList
    return MasterChannelList, mdfdict

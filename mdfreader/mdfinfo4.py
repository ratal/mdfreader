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
from __future__ import absolute_import  # for consistency between python 2 and 3
from __future__ import print_function
from struct import calcsize, unpack, pack, Struct
from os import remove
from sys import version_info
from warnings import warn
from zlib import compress, decompress
from numpy import sort, zeros, array, append
import time
from xml.etree.ElementTree import Element, SubElement, \
    tostring, register_namespace
from lxml import objectify
from .mdf import _open_MDF, dataField, descriptionField, unitField, \
    masterField, masterTypeField, idField, _convertName

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
DZStruct = Struct('<2s2BI2Q')
HLStruct = Struct('<QHB5s')
SIStruct = Struct('<4sI5Q3B5s')
DGStruct = Struct('<4sI6QB7s')
CGStruct = Struct('<4sI10Q2H3I')
CNStruct = Struct('<4B4IBcH6d')
CCStruct1 = Struct('<4sI6Q')
CCStruct2 = Struct('<2B3H2d')
SRStruct = Struct('<4sI5Qd2B6s')

chunk_size_writing = 4194304  # write by chunk of 4Mb, can be tuned for best performance


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


def _calculate_block_start(current_position):
    """ converts pointer poistion into one being multiple of 8

        Parameters
        ----------------
        current_position : int
            position in file

        Returns
        -----------
        position multiple of 8
    """
    remain = current_position % 8
    if not remain == 0:
        return current_position - remain + 8
    else:
        return current_position


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
        if self['id_ver'] > 410 and b'UnFin' in self['id_file']:
            warn('  ! unfinalised file')
            if self['id_unfi_flags'] & 1:
                warn('Update of cycle counters for CG/CA blocks required')
            if self['id_unfi_flags'] & (1 << 1):
                warn('Update of cycle counters for SR blocks required')
            if self['id_unfi_flags'] & (1 << 2):
                warn('Update of length for last DT block required')
            if self['id_unfi_flags'] & (1 << 3):
                warn('Update of length for last RD block required')
            if self['id_unfi_flags'] & (1 << 4):
                warn('Update of last DL block in each chained list of '
                     'DL blocks required')
            if self['id_unfi_flags'] & (1 << 5):
                warn('Update of cg_data_bytes and cg_inval_bytes in VLSD '
                     'CG block required')
            if self['id_unfi_flags'] & (1 << 6):
                warn('Update of offset values for VLSD channel required '
                     'in case a VLSD CG block is used')

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
            comment = CommentBlock()
            comment.read(fid=fid, pointer=self['hd_md_comment'], MDType='HD')
            self['Comment'].update(comment)

    def write(self, fid):
        # write block header
        fid.seek(self['block_start'])
        # link section
        # (first Data group pointer, first file history block pointer,
        # pointer to hierarchy Block file, pointer to attachment Block,
        # pointer to event Block, pointer to comment Block,
        # time in ns, timezone offset in min, time daylight offset in min,
        # time flags, time class, hd flags, reserved, start angle in radians
        # start distance in meters)
        dataBytes = (b'##HD', 0, 104, 6,
                     self['DG'], self['FH'], 0, 0, 0, 0,
                     int(time.time()*1E9),
                     int(time.timezone/60),
                     int(time.daylight*60),
                     2, 0, 0, b'\0', 0, 0)
        fid.write(pack('<4sI2Q7Q2h3Bs2d', *dataBytes))


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
            comment = CommentBlock()
            comment.read(fid=fid, pointer=self['fh_md_comment'], MDType='FH')
            self['Comment'].update(comment)

    def write(self, fid):
        # write block header
        fid.seek(self['block_start'])
        # link section
        # (No next FH, comment block pointer
        # time in ns, timezone offset in min, time daylight offest in min,
        # time flags, reserved)
        dataBytes = (b'##FH', 0, 56, 2,
                     0, self['MD'],
                     int(time.time()*1E9),
                     int(time.timezone/60),
                     int(time.daylight*60),
                     2, b'\0' * 3)
        fid.write(pack('<4sI2Q3Q2hB3s', *dataBytes))


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
            comment= CommentBlock()
            comment.read(fid=fid, pointer=self['ch_md_comment'])
            self['Comment'].update(comment)
        if self['ch_tx_name']:  # text block containing name of hierarchy level
            self['ch_name_level'] = {}
            comment = CommentBlock()
            comment.read(fid=fid, pointer=self['ch_tx_name'])
            self['ch_name_level'].update(comment)


class CommentBlock(dict):
    """ reads or writes Comment block and saves in class dict
    """

    def read(self, **kargs):
        """ reads Comment block and saves in class dict
        Parameters
        ----------
        fid: file identifier
        pointer: int
            position in file
        MDType: str
            describes metadata type, ('CN', 'unit', 'FH', 'SI', 'HD', 'CC', 'EV')

        Notes
        --------
        Can read xml (MD metadata) or text (TX) comments from several
        kind of blocks
        """

        if kargs['pointer'] > 0:
            kargs['fid'].seek(kargs['pointer'])
            # block header
            (self['id'],
             reserved,
             self['length'],
             self['link_count']) = HeaderStruct.unpack(kargs['fid'].read(24))
            if self['id'] in ('##MD', b'##MD'):
                # Metadata block
                # removes normal 0 at end
                # self['Comment'] = None
                try:
                    xml_tree = objectify.fromstring(kargs['fid'].read(self['length'] - 24).rstrip(b'\x00'))
                except:
                    warn('xml metadata malformed')
                    xml_tree = None
                # specific action per comment block type,
                # #extracts specific tags from xml
                if 'MDType' in kargs:
                    if kargs['MDType'] == 'CN':  # channel comment
                        try:
                            self['description'] = xml_tree.TX.text
                        except AttributeError:
                            warn('Could not parse CN block TX tag')
                        try:
                            self['names'] = xml_tree.names.text
                        except AttributeError:
                            pass  # optional
                        if 'minimal' in kargs and kargs['minimal'] is False:
                            # not really used for the moment
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
                    elif kargs['MDType'] == 'unit':  # channel comment
                        try:
                            self['unit'] = xml_tree.TX
                        except AttributeError:
                            warn('Could not parse unit TX tag')
                    elif kargs['MDType'] == 'HD':  # header comment
                        try:
                            self['TX'] = xml_tree.TX.text
                        except AttributeError:
                            warn('Could not parse HD block TX tag')
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
                    elif kargs['MDType'] == 'FH':  # File History comment
                        try:
                            self['TX'] = xml_tree.TX.text
                        except AttributeError:
                            warn('Could not parse FH block TX tag')
                        try:
                            self['tool_id'] = xml_tree.tool_id.text
                        except AttributeError:
                            warn('Could not parse FH block tool_id tag')
                        try:
                            self['tool_vendor'] = xml_tree.tool_vendor.text
                        except AttributeError:
                            warn('Could not parse extract HD block tool_vendor tag')
                        try:
                            self['tool_version'] = xml_tree.tool_version.text
                        except AttributeError:
                            warn('Could not parse extract HD block tool_vendor tag')
                        try:
                            self['user_name'] = xml_tree.user_name.text
                        except AttributeError:
                            pass  # optional
                    elif kargs['MDType'] == 'SI':
                        try:
                            self['TX'] = xml_tree.TX.text
                        except AttributeError:
                            warn('Could not parse SI block TX tag')
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
                    elif kargs['MDType'] == 'CC':
                        try:
                            self['TX'] = xml_tree.TX.text
                        except AttributeError:
                            warn('Could not parse CC block TX tag')
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
                    elif kargs['MDType'] == 'EV':
                        try:
                            self['TX'] = xml_tree.TX.text
                        except AttributeError:
                            warn('Could not parse EV block TX tag')
                        try:
                            self['pre_trigger_interval'] = xml_tree.pre_trigger_interval.text
                        except AttributeError:
                            pass  # optional
                        try:
                            self['post_trigger_interval'] = xml_tree.post_trigger_interval.text
                        except AttributeError:
                            pass  # optional
                        try:
                            self['formula'] = xml_tree.formula.text
                        except AttributeError:
                            pass  # optional
                        try:
                            self['timeout'] = xml_tree.timeout.text
                        except AttributeError:
                            pass  # optional
                    else:
                        if kargs['MDType'] is not None:
                            warn('No recognized MDType {}'.format(kargs['MDType']))
            elif self['id'] in ('##TX', b'##TX'):
                if 'MDType' in kargs and kargs['MDType'] == 'CN':  # channel comment
                    self['name'] = kargs['fid'].read(self['length'] - 24).rstrip(b'\x00').decode('UTF-8', 'ignore')
                else:
                    self['Comment'] = kargs['fid'].read(self['length'] - 24).rstrip(b'\x00').decode('UTF-8', 'ignore')

    def load(self, data, MDType):
        if MDType == 'TX':
            data = b''.join([data.encode('utf-8', 'replace'), b'\0'])
            block_id = b'##TX'
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
            block_id = b'##MD'
        # make sure block is multiple of 8
        remain = len(data) % 8
        if not remain == 0:
            data = b''.join([data, b'\0' * (8 - (remain % 8))])
        self['block_length'] = 24 + len(data)
        self['Header'] = (block_id, 0, self['block_length'], 0)
        self['data'] = data

    def write(self, fid):
        fid.write(HeaderStruct.pack(*self['Header']))
        fid.write(self['data'])


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
            comment = CommentBlock()
            comment.read(fid=fid, pointer=self['dg_md_comment'])
            self['Comment'].update(comment)

    def write(self, fid):
        # write block header
        fid.seek(self['block_start'])
        # link section
        # (Next Data group pointer, first channel group block pointer,
        # data block pointer, comment block pointer,
        # no recordID, reserved)
        dataBytes = (b'##DG', 0, 64, 4,
                     self['DG'], self['CG'],
                     self['data'], 0, 0, b'\0' * 7)
        fid.write(pack('<4sI2Q4QB7s', *dataBytes))


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
            self['Comment'] = CommentBlock()
            self['Comment'].read(fid=fid, pointer=self['cg_md_comment'])
        if self['cg_tx_acq_name']:  # comments exist
            self['acq_name'] = {}
            comment = CommentBlock()
            comment.read(fid=fid, pointer=self['cg_tx_acq_name'])
            self['acq_name'].update(comment)

    def write(self, fid):
        # write block header
        #fid.seek(self['block_start'])
        # link section
        # (Next channel group pointer, first channel block pointer,
        # acquisition name pointer, acquisition source pointer,
        # sample reduction pointer, comment pointer,
        # no recordID, cycle count, no flags,
        # no character specified for path separator,
        # reserved, number of bytes taken by data in record,
        # no invalid bytes)
        dataBytes = (b'##CG', 0, 104, 6,
                     0, self['CN'], 0, 0, 0, 0, 0,
                     self['cg_cycle_count'], 0, 0,
                     b'\0' * 4,
                     self['cg_data_bytes'], 0)
        fid.write(pack('<4sI2Q8Q2H4s2I', *dataBytes))


class CNBlock(dict):

    """ reads Channel block and saves in class dict
    """

    def read(self, **kargs):
        if kargs['pointer'] != 0 and kargs['pointer'] is not None:
            kargs['fid'].seek(kargs['pointer'])
            (self['id'],
             reserved,
             self['length'],
             self['link_count']) = HeaderStruct.unpack(kargs['fid'].read(24))
            self['pointer'] = kargs['pointer']
            # data section
            kargs['fid'].seek(kargs['pointer'] + self['length'] - 72)
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
             cn_limit_min,
             cn_limit_max,
             cn_limit_ext_min,
             cn_limit_ext_max) = CNStruct.unpack(kargs['fid'].read(72))
            # Channel Group block : Links
            kargs['fid'].seek(kargs['pointer'] + 24)
            (self['cn_cn_next'],
             self['cn_composition'],
             self['cn_tx_name'],
             self['cn_si_source'],
             self['cn_cc_conversion'],
             self['cn_data'],
             self['cn_md_unit'],
             self['cn_md_comment']) = unpack('<8Q', kargs['fid'].read(64))
            if self['cn_attachment_count'] > 0:
                self['cn_at_reference'] = \
                    _mdfblockread(kargs['fid'], LINK, self['cn_attachment_count'])
                self['attachment'] = {}
                if self['cn_attachment_count'] > 1:
                    for at in range(self['cn_attachment_count']):
                        self['attachment'][at] = \
                            ATBlock(kargs['fid'], self['cn_at_reference'][at])
                else:
                    self['attachment'][0] = \
                        ATBlock(kargs['fid'], self['cn_at_reference'])
            if self['link_count'] > (8 + self['cn_attachment_count']):
                self['cn_default_x'] = _mdfblockread(kargs['fid'], LINK, 3)
            else:
                self['cn_default_x'] = None
            if self['cn_md_comment']:  # comments exist
                self['Comment'] = CommentBlock()
                self['Comment'].read(fid=kargs['fid'], pointer=self['cn_md_comment'], MDType='CN', minimal=True)
            if self['cn_md_unit']:  # comments exist
                self['unit'] = CommentBlock()
                self['unit'].read(fid=kargs['fid'], pointer=self['cn_md_unit'], MDType='unit')
            if self['cn_tx_name']:  # comments exist
                comment = CommentBlock()
                comment.read(fid=kargs['fid'], pointer=self['cn_tx_name'], MDType='CN')
                self['name'] = comment['name']

    def write(self, fid):
        # write block header
        # fid.seek(self['block_start'])
        # no attachement and default X
        # link section
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
        dataBytes = (b'##CN', 0, 160, 8,
                     self['CN'], self['Composition'], self['TX'], 0, 0, 0, self['Unit'], self['Comment'],
                     self['cn_type'], self['cn_sync_type'],
                     self['cn_data_type'], self['cn_bit_offset'],
                     self['cn_byte_offset'], self['cn_bit_count'],
                     self['cn_flags'], 0, 0, 0, 0,
                     self['cn_val_range_min'], self['cn_val_range_max'],
                     0, 0, 0, 0)
        fid.write(pack('<4sI2Q8Q4B4I2BH6d', *dataBytes))


class CCBlock(dict):

    """ reads Channel Conversion block and saves in class dict
    """

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
                cc = CommentBlock()
                cc.read(fid=fid, pointer=pointer)
                self['cc_ref'].update(cc)
            elif self['cc_type']in (7, 8, 9, 10):  # text list
                self['cc_ref'] = list(self['cc_ref'])
                for i in range(self['cc_ref_count']):
                    fid.seek(self['cc_ref'][i])
                    # find if TX/MD or another CCBlock
                    ID = unpack('4s', fid.read(4))[0]
                    # for algebraic formulae
                    if ID in ('##TX', '##MD', b'##TX', b'##MD'):
                        temp = CommentBlock()
                        temp.read(fid=fid, pointer=self['cc_ref'][i])
                        self['cc_ref'][i] = temp['Comment']
                    elif ID in ('##CC', b'##CC'):  # for table conversion
                        # much more complicated nesting conversions !!!
                        cc = CCBlock()
                        cc.read(fid, self['cc_ref'][i])
                        self['cc_ref'][i] = cc
            if self['cc_md_comment']:  # comments exist
                self['Comment'] = CommentBlock()
                self['Comment'].read(fid=fid, pointer=self['cc_md_comment'], MDType='CC')
            if self['cc_md_unit']:  # comments exist
                self['unit'] = CommentBlock()
                self['unit'].read(fid=fid, pointer=self['cc_md_unit'])
            if self['cc_tx_name']:  # comments exist
                self['name'] = CommentBlock()
                self['name'].read(fid=fid, pointer=self['cc_tx_name'])
        else:  # no conversion
            self['cc_type'] = 0


class CABlock(dict):

    """ reads Channel Array block and saves in class dict
    """

    def read(self, fid, pointer):
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
                self['CABlock'] = CABlock()
                self['CABlock'].read(fid, self['ca_composition'])

    def load(self):
        self['block_length'] = 48 + 8 * self['ndim']

    def write(self, fid):
        # default CN template
        dataBytes = (b'##CA', 0, 48 + 8 * self['ndim'], 1,
                     0,
                     0, 0, self['ndim'], 0, 0, 0)
        fid.write(pack('<4sI2Q1Q2BHIiI', *dataBytes))
        fid.write(pack('<{}Q'.format(self['ndim']), *self['ndim_size']))


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
                comment = CommentBlock()
                comment.read(fid=fid, pointer=self['at_md_comment'])
                self['Comment'].update(comment)


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
                self['Comment'] = CommentBlock()
                self['Comment'].read(fid=fid, pointer=self['ev_md_comment'], MDType='EV')
            if self['ev_tx_name']:  # comments exist
                self['name'] = CommentBlock()
                self['name'].read(fid=fid, pointer=self['ev_tx_name'])


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
            self['source_name'] = CommentBlock()
            self['source_name'].read(fid=fid, pointer=self['si_tx_name'])
            self['source_path'] = CommentBlock()
            self['source_path'].read(fid=fid, pointer=self['si_tx_path'])
            self['comment'] = CommentBlock()
            self['comment'].read(fid=fid, pointer=self['si_md_comment'], MDType='SI')


class DTBlock(dict):
    def load(self, record_byte_offset, nRecords, pointer):
        self['datablocks_length'] = 24 + record_byte_offset * nRecords
        self['pointer'] = pointer
        self['end_position'] = _calculate_block_start(self['pointer'] + self['datablocks_length'])

    def write(self, fid, data):
        fid.write(HeaderStruct.pack(b'##DT', 0, self['datablocks_length'], 0))
        # dumps data
        fid.write(data)
        return self['end_position']


class DLBlock(dict):

    """ reads Data List block
    """

    def read(self, fid, link_count):
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

    def write(self, fid, chunks, position):
        self['block_start'] = position
        number_DL = len(chunks)
        self['block_length'] = 40 + 16 * number_DL
        dl_offset = zeros(shape=number_DL, dtype='<u8')
        if number_DL > 1:
            for counter in range(1, number_DL):
                (nrecord_chunk, chunk_size) = chunks[counter]
                dl_offset[counter] = dl_offset[counter - 1] + chunk_size
        dataBytes = (b'##DL', 0, self['block_length'], number_DL + 1, 0)
        fid.write(pack('<4sI3Q'.format(number_DL, number_DL), *dataBytes))
        fid.write(pack('{0}Q'.format(number_DL), *zeros(shape=number_DL, dtype='<u8')))
        fid.write(pack('<2I', 0, number_DL))
        fid.write(pack('{0}Q'.format(number_DL), *dl_offset))


class DZBlock(dict):

    """ reads Data List block
    """

    def read(self, fid):
        # block header is already read
        (self['dz_org_block_type'],
         self['dz_zip_type'],
         dz_reserved,
         self['dz_zip_parameter'],
         self['dz_org_data_length'],
         self['dz_data_length']) = DZStruct.unpack(fid.read(24))

    @staticmethod
    def decompress_datablock(block, zip_type, zip_parameter, org_data_length):
        """ decompress datablock.

        Parameters
        --------------
        block : bytes
            raw data compressed
        zip_type : int
            0 for non transposed, 1 for transposed data
        zip_parameter : int
            first dimension of matrix to be transposed
        org_data_length : int
            uncompressed data length

        Returns
        ---------
        uncompressed raw data
        """

        if PythonVersion > 3:  # decompress data
            block = decompress(block)
        else:
            block = decompress(bytes(block))  # decompress data
        if zip_type == 1:  # data bytes transposed
            M = org_data_length // zip_parameter
            temp = array(memoryview(block[:M * zip_parameter]))
            tail = array(memoryview(block[M * zip_parameter:]))
            temp = temp.reshape(zip_parameter, M).T.ravel()
            if len(tail) > 0:
                temp = append(temp, tail)
            block = temp.tostring()
        return block

    def write(self, fid, data, record_length):
        fid.seek(self['block_start'])
        # compress data
        # uses data transposition for max compression
        org_data_length = len(data)
        M = org_data_length // record_length
        temp = array(memoryview(data[:M * record_length]))
        tail = array(memoryview(data[M * record_length:]))
        temp = temp.reshape(M, record_length).T.ravel()
        if len(tail) > 0:
            temp = append(temp, tail)
        # compress transposed data
        compressed_data = compress(temp.tostring())
        dz_data_length = len(compressed_data)
        self['DZBlock_length'] = 48 + dz_data_length
        # writes header
        fid.write(HeaderStruct.pack(b'##DZ', 0, self['DZBlock_length'], 0))
        dataBytes = (b'DT', 1, 0, record_length, org_data_length, dz_data_length)
        fid.write(DZStruct.pack(*dataBytes))
        # writes compressed data
        fid.write(compressed_data)
        return self['block_start'] + self['DZBlock_length']


class HLBlock(dict):

    """ reads Header List block
    """

    def read(self, fid):
        (self['hl_dl_first'],
         self['hl_flags'],
         self['hl_zip_type'],
         hl_reserved) = HLStruct.unpack(fid.read(16))

    def load(self, record_byte_offset, nRecords, position):
        self['record_length'] = record_byte_offset
        self['block_start'] = position
        self['block_length'] = 40
        # calculate data chunks
        nchunks = (record_byte_offset * nRecords) // chunk_size_writing + 1
        chunk_length = (record_byte_offset * nRecords) // nchunks
        nrecord_chunk = chunk_length // record_byte_offset
        self['chunks'] = [(nrecord_chunk, record_byte_offset * nrecord_chunk)] * nchunks
        nrecord_chunk = nRecords - nrecord_chunk * nchunks
        if nrecord_chunk > 0:
            self['chunks'].append((nrecord_chunk, record_byte_offset * nrecord_chunk))

    def write(self, fid, data):
        fid.write(HeaderStruct.pack(b'##HL', 0, self['block_length'], 1))
        # uses not equal length blocks, transposed compressed data
        fid.write(HLStruct.pack(self['block_start'] + self['block_length'], 0, 1, b'\x00'*5))
        pointer = self['block_start'] + self['block_length']
        DL = DLBlock()
        DL.write(fid, self['chunks'], pointer)
        pointer += DL['block_length']
        dl_data = zeros(shape=len(self['chunks']), dtype='<u8')
        data_pointer = 0
        for counter, (nrecord_chunk, chunk_size) in enumerate(self['chunks']):
            DZ = DZBlock()
            DZ['block_start'] = _calculate_block_start(pointer)
            dl_data[counter] = DZ['block_start']
            pointer = DZ.write(fid, data[data_pointer: data_pointer + chunk_size], self['record_length'])
            data_pointer += chunk_size
        # writes links to all DZBlocks
        # write dl_data
        fid.seek(DL['block_start'] + 32)
        fid.write(dl_data.tobytes())
        fid.seek(pointer)
        return _calculate_block_start(pointer)


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
    - mdfinfo['CC'][dataGroup][channelGroup][channel]"""

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
        Either fileName or fid can be used as argument"""

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
            # warn('reads File History blocks, always exists')
            fh = 0  # index of fh blocks
            self['FH'][fh] = {}
            self['FH'][fh] .update(FHBlock(fid, self['HD']['hd_fh_first']))
            while self['FH'][fh]['fh_fh_next']:
                self['FH'][fh + 1] = {}
                self['FH'][fh + 1].update(FHBlock(fid, self['FH'][fh]['fh_fh_next']))
                fh += 1

        # warn('reads Channel Hierarchy blocks')
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
            self['EV'] = {}
            self['EV'][ev] = EVBlock(fid, self['HD']['hd_ev_first'])
            while self['EV'][ev]['ev_ev_next']:
                ev += 1
                self['EV'][ev] = EVBlock(fid, self['EV'][ev - 1]['ev_ev_next'])

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
        temp.read(fid=fid, pointer=self['CG'][dg][cg]['cg_cn_first'])
        if temp is not None:
            self['CN'][dg][cg][cn].update(temp)
        MLSDChannels = []
        # check for MLSD
        if self['CN'][dg][cg][cn]['cn_type'] == 5:
            MLSDChannels.append(cn)
        # check if already existing channel name
        if self['CN'][dg][cg][cn]['name'] in self['ChannelNamesByDG'][dg]:  # for unsorted data
            # reads Channel Source Information
            if self['CN'][dg][cg][cn]['cn_si_source']:
                temp = SIBlock()
                temp.read(fid, self['CN'][dg][cg][cn]['cn_si_source'])
                if temp['si_tx_name'] > 0:
                    source_name = temp['source_name']['Comment']
                else:
                    source_name = cn
            else:
                source_name = cn
            self['CN'][dg][cg][cn]['name'] = \
                u'{0}_{1}_{2}_{3}'.format(self['CN'][dg][cg][cn]['name'], dg, cg, source_name)
        elif self['CN'][dg][cg][cn]['name'] in self['allChannelList']:
            # doublon name or master channel
            if self['CN'][dg][cg][cn]['cn_si_source']:
                temp = SIBlock()
                temp.read(fid, self['CN'][dg][cg][cn]['cn_si_source'])
                if temp['si_tx_name'] > 0:
                    source_name = temp['source_name']['Comment']
                else:
                    source_name = dg
            else:
                source_name = dg
            self['CN'][dg][cg][cn]['name'] = u'{0}_{1}_{2}'.format(self['CN'][dg][cg][cn]['name'], dg, source_name)
        if self['CN'][dg][cg][cn]['cn_type'] == 1 and PythonVersion < 3:
            # VLSD needs to rename and append records but with python 2.x impossible,
            # convert name to compatible python identifier
            vlsd = True
        self['ChannelNamesByDG'][dg].add(self['CN'][dg][cg][cn]['name'])
        self['allChannelList'].add(self['CN'][dg][cg][cn]['name'])

        if self['CG'][dg][cg]['cg_cn_first']:  # Can be NIL for VLSD
            # reads Channel Conversion Block
            self['CC'][dg][cg][cn] = CCBlock()
            self['CC'][dg][cg][cn].read(fid, self['CN'][dg][cg][cn]['cn_cc_conversion'])
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
                        self['CN'][dg][cg][cn]['CABlock'] = CABlock()
                        self['CN'][dg][cg][cn]['CABlock'].read(fid, self['CN'][dg][cg][cn]['cn_composition'])
                    elif id in ('##CN', b'##CN'):
                        self['CN'][dg][cg][cn]['CN'] = {}
                        temp = CNBlock()
                        temp.read(fid=fid, pointer=self['CN'][dg][cg][cn]['cn_composition'])
                        self['CN'][dg][cg][cn]['CN'].update(temp)
                    else:
                        raise('unknown channel composition')

                if not minimal:
                    # reads Attachment Block
                    if self['CN'][dg][cg][cn]['cn_attachment_count'] > 1:
                        for at in range(self['CN'][dg][cg][cn]['cn_attachment_count']):
                            self['CN'][dg][cg][cn]['attachment'][at].update(self.readATBlock(fid, self['CN'][dg][cg][cn]['cn_at_reference'][at]))
                    elif self['CN'][dg][cg][cn]['cn_attachment_count'] == 1:
                        self['CN'][dg][cg][cn]['attachment'][0].update(
                            self.readATBlock(fid, self['CN'][dg][cg][cn]['cn_at_reference']))

            while self['CN'][dg][cg][cn]['cn_cn_next']:
                cn = cn + 1
                self['CN'][dg][cg][cn] = {}
                temp = CNBlock()
                temp.read(fid=fid, pointer=self['CN'][dg][cg][cn - 1]['cn_cn_next'])
                self['CN'][dg][cg][cn].update(temp)
                # check for MLSD
                if self['CN'][dg][cg][cn]['cn_type'] == 5:
                    MLSDChannels.append(cn)
                # reads Channel Conversion Block
                self['CC'][dg][cg][cn] = CCBlock()
                self['CC'][dg][cg][cn].read(fid, self['CN'][dg][cg][cn]['cn_cc_conversion'])
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
                        # reads Channel Source Information
                        if self['CN'][dg][cg][cn]['cn_si_source']:
                            temp = SIBlock()
                            temp.read(fid, self['CN'][dg][cg][cn]['cn_si_source'])
                            if temp['si_tx_name'] > 0:
                                source_name = temp['source_name']['Comment']
                            else:
                                source_name = cn
                        else:
                            source_name = cn
                        self['CN'][dg][cg][cn]['name'] = \
                            u'{0}_{1}_{2}_{3}'.format(self['CN'][dg][cg][cn]['name'], dg, cg, source_name)
                    elif self['CN'][dg][cg][cn]['name'] in self['allChannelList']:
                        # doublon name or master channel
                        # reads Channel Source Information
                        if self['CN'][dg][cg][cn]['cn_si_source']:
                            temp = SIBlock()
                            temp.read(fid, self['CN'][dg][cg][cn]['cn_si_source'])
                            if temp['si_tx_name'] > 0:
                                source_name = temp['source_name']['Comment']
                            else:
                                source_name = dg
                        else:
                            source_name = dg
                        self['CN'][dg][cg][cn]['name'] = u'{0}_{1}_{2}'.format(self['CN'][dg][cg][cn]['name'], dg,
                                                                               source_name)
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
                            self['CN'][dg][cg][cn]['CABlock'] = CABlock()
                            self['CN'][dg][cg][cn]['CABlock'].read(fid, self['CN'][dg][cg][cn]['cn_composition'])
                        elif id in ('##CN', b'##CN'):
                            self['CN'][dg][cg][cn]['CN'] = {}
                            temp = CNBlock()
                            temp.read(fid=fid, pointer=self['CN'][dg][cg][cn]['cn_composition'])
                            self['CN'][dg][cg][cn]['CN'].update(temp)
                        else:
                            raise('unknown channel composition')

                    if not minimal:
                        # reads Attachment Block
                        if self['CN'][dg][cg][cn]['cn_attachment_count'] > 1:
                            for at in range(self['CN'][dg][cg][cn]['cn_attachment_count']):
                                self['CN'][dg][cg][cn]['attachment'][at].update(self.readATBlock(fid, self['CN'][dg][cg][cn]['cn_at_reference'][at]))
                        elif self['CN'][dg][cg][cn]['cn_attachment_count'] == 1:
                            self['CN'][dg][cg][cn]['attachment'][0].update(
                                self.readATBlock(fid, self['CN'][dg][cg][cn]['cn_at_reference']))

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
                    self['CN'][dg][cg][chan].read(fid=fid,
                                                  pointer=self['CN'][dg][cg][cn]['cn_composition'])
                    self['CC'][dg][cg][chan] = CCBlock()
                    self['CC'][dg][cg][chan].read(fid, self['CN'][dg][cg][chan]['cn_cc_conversion'])

                    if self['CN'][dg][cg][chan]['cn_type'] == 5:
                        MLSDChannels.append(chan)
                    while self['CN'][dg][cg][chan]['cn_cn_next']:
                        chan += 1
                        self['CN'][dg][cg][chan] = CNBlock()
                        self['CN'][dg][cg][chan].read(fid=fid, pointer=self['CN'][dg][cg][chan - 1]['cn_cn_next'])

                        self['CC'][dg][cg][chan] = CCBlock()
                        self['CC'][dg][cg][chan].read(fid, self['CN'][dg][cg][chan]['cn_cc_conversion'])
                        if self['CN'][dg][cg][chan]['cn_type'] == 5:
                            MLSDChannels.append(chan)
                    # makes the channel virtual
                    self['CN'][dg][cg][cn]['cn_type'] = 6
                elif ID in ('##CA', b'##CA'):  # arrays
                    pass
                else:
                    warn('unknown channel composition')
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
        ChannelNamesByDG = set()
        for cg in info['CG'][dg]:
            channelNameList = []
            for cn in info['CN'][dg][cg]:
                if info['CN'][dg][cg][cn]['cn_data_type'] == 13:
                    for name in ('ms', 'minute', 'hour', 'day', 'month', 'year'):
                        channelNameList.append(name)
                        allChannelList.add(name)
                        ChannelNamesByDG.add(name)
                        mdfdict[name] = {}
                        mdfdict[name][dataField] = None
                        mdfdict[name][descriptionField] = name
                        mdfdict[name][unitField] = name
                        mdfdict[name][idField] = (dg, cg, cn)
                elif info['CN'][dg][cg][cn]['cn_data_type'] == 14:
                    for name in ('ms', 'days'):
                        channelNameList.append(name)
                        allChannelList.add(name)
                        ChannelNamesByDG.add(name)
                        mdfdict[name] = {}
                        mdfdict[name][dataField] = None
                        mdfdict[name][descriptionField] = name
                        mdfdict[name][unitField] = name
                        mdfdict[name][idField] = (dg, cg, cn)
                else:
                    name = info['CN'][dg][cg][cn]['name']
                    if name in ChannelNamesByDG:
                        name = u'{0}_{1}_{2}_{3}'.format(name, dg, cg, cn)
                    elif name in allChannelList:
                        name = u'{0}_{1}'.format(name, dg)
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
                        mdfdict[name][idField] = (dg, cg, cn)
                    if info['CN'][dg][cg][cn]['cn_sync_type']:
                        # master channel of cg
                        master = name
                        mastertype = info['CN'][dg][cg][cn]['cn_sync_type']
            for chan in channelNameList:
                mdfdict[chan][masterField] = master
                mdfdict[chan][masterTypeField] = mastertype
        MasterChannelList[master] = channelNameList
    return MasterChannelList, mdfdict

# -*- coding: utf-8 -*-
""" Measured Data Format blocks parser for version 4.x

Created on Sun Dec 15 12:57:28 2013

:Author: `Aymeric Rateau <https://github.com/ratal/mdfreader>`__

mdfinfo4
--------------------------
"""
from struct import calcsize, unpack, pack, Struct
from os import remove
from warnings import warn
from zlib import compress, decompress
from numpy import zeros, array, append
from math import isnan
from time import time
from sys import getsizeof
from collections import OrderedDict
from xml.etree.ElementTree import Element, SubElement, \
    tostring, register_namespace
from lxml import objectify
from .mdf import _open_mdf, dataField, descriptionField, unitField, \
    masterField, masterTypeField, idField, _convert_name

# datatypes
_LINK = '<Q'
_REAL = '<d'
_BOOL = '<h'
_UINT8 = '<B'
_BYTE = '<c'
_INT16 = '<h'
_UINT16 = '<H'
_UINT32 = '<I'
_INT32 = '<i'
_UINT64 = '<Q'
_INT64 = '<q'
_CHAR = '<c'

_HeaderStruct = Struct('<4sI2Q')
_DZStruct = Struct('<2s2BI2Q')
_HLStruct = Struct('<QHB5s')
_SIStruct = Struct('<4sI5Q3B5s')
_DGStruct = Struct('<4sI6QB7s')
_CGStruct1 = Struct('<4sI8Q')
_CGStruct2 = Struct('<2Q2H3I')
_CNStruct1 = Struct('<4sI10Q')
_CNStruct2 = Struct('<4B4IBcH6d')
_CCStruct1 = Struct('<4sI6Q')
_CCStruct2 = Struct('<2B3H2d')
_SRStruct = Struct('<4sI5Qd2B6s')

# SI Types
si_type = {0: 'OTHER', 1: 'ECU', 2: 'BUS',
           3: 'I/O', 4: 'TOOL', 5: 'USER'}
si_bus_type = {0: 'NONE', 1: 'OTHER', 2: 'CAN', 3: 'LIN',
               4: 'MOST', 5: 'FLEXRAY', 6: 'K_LINE', 7: 'ETHERNET', 8: 'USB'}
# EV cause
ev_cause = {0: 'OTHER', 1: 'ERROR', 2: 'TOOL', 3: 'SCRIPT', 4: 'USER'}

# precompiled path for xmls to speed up parsing
CN_TX = objectify.ObjectPath('CNcomment.TX')
CN_names = objectify.ObjectPath('CNcomment.names')
CN_axis_monotony = objectify.ObjectPath('CNcomment.axis_monotony')
CN_raster = objectify.ObjectPath('CNcomment.raster')
CN_formula = objectify.ObjectPath('CNcomment.formula')
CN_linker_name = objectify.ObjectPath('CNcomment.linker_name')
CN_linker_address = objectify.ObjectPath('CNcomment.linker_address')
CN_address = objectify.ObjectPath('CNcomment.address')
CN_unit_TX = objectify.ObjectPath('CNunit.TX')
CC_TX = objectify.ObjectPath('CCcomment.TX')
CC_names = objectify.ObjectPath('CCcomment.names')
CC_COMPU_METHOD = objectify.ObjectPath('CCcomment.COMPU_METHOD')
CC_formula = objectify.ObjectPath('CCcomment.formula')
CC_unit_TX = objectify.ObjectPath('CCunit.TX')
SI_TX = objectify.ObjectPath('SIcomment.TX')
SI_names = objectify.ObjectPath('SIcomment.names')
SI_path = objectify.ObjectPath('SIcomment.path')
SI_bus = objectify.ObjectPath('SIcomment.bus')
SI_protocol = objectify.ObjectPath('SIcomment.protocol')
EV_TX = objectify.ObjectPath('EVcomment.TX')
EV_pre_trigger_interval = objectify.ObjectPath('EVcomment.pre_trigger_interval')
EV_post_trigger_interval = objectify.ObjectPath('EVcomment.post_trigger_interval')
EV_formula = objectify.ObjectPath('EVcomment.formula')
EV_timeout = objectify.ObjectPath('EVcomment.timeout')

chunk_size_writing = 4194304  # write by chunk of 4Mb, can be tuned for best performance


def _load_header(fid, pointer):
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
         temp['link_count']) = _HeaderStruct.unpack(fid.read(24))
        temp['pointer'] = pointer
        return temp
    else:
        return None


def _mdf_block_read(fid, data_type, count):
    """ converts a byte array of length count to a given data Type

    Parameters
    ----------------
    data_type : str
        C format data type
    count : int
        number of elements to sequentially read

    Returns
    -----------
    array of values of 'Type' parameter
    """
    value = fid.read(calcsize(data_type) * count)
    if value:
        if count == 1:
            return unpack(data_type, value)[0]
        else:
            if '<' in data_type or '>' in data_type:
                endian = data_type[0]
                data_type = data_type.strip('<>')
                return unpack(endian + count * data_type, value)
            else:
                return unpack(count * data_type, value)
    else:
        return None


def _calculate_block_start(current_position):
    """ converts pointer position into one being multiple of 8

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
        head = (b'MDF     ', self['id_vers'], b'MDFreadr', b'\0' * 4, self['id_ver'],
                b'\0' * 30, 0, 0)
        fid.write(pack('<8s8s8s4sH30s2H', *head))


class HDBlock(dict):

    """ reads Header block and save in class dict
    """

    def __init__(self, fid=None):
        self['time'] = time()
        if fid is not None:
            self.read(fid)

    def read(self, fid=None):
        fid.seek(64)
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
         self['hd_start_distance']) = unpack('<4sI9Q2h4B2d', fid.read(104))
        if self['hd_md_comment']:  # if comments exist
            self['Comment'] = {}
            comment = CommentBlock()
            comment.read_cm_hd(fid=fid, pointer=self['hd_md_comment'])
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
        data_bytes = (b'##HD', 0, 104, 6,
                      self['DG'], self['FH'], 0, 0, 0, 0,
                      int(self['time'] * 1E9),
                      0, 0, 1, 0, 0, b'\0', 0, 0)
        try:
            fid.write(pack('<4sI2Q7Q2h3Bs2d', *data_bytes))
        except:
            warn('time is Nan or malformed')
            data_bytes = (b'##HD', 0, 104, 6,
                          self['DG'], self['FH'], 0, 0, 0, 0,
                          int(time() * 1E9),
                          0, 0, 1, 0, 0, b'\0', 0, 0)
            fid.write(pack('<4sI2Q7Q2h3Bs2d', *data_bytes))


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
            comment.read_cm_fh(fid=fid, pointer=self['fh_md_comment'])
            self['Comment'].update(comment)

    def write(self, fid):
        # write block header
        fid.seek(self['block_start'])
        # link section
        # (No next FH, comment block pointer
        # time in ns, timezone offset in min, time daylight offest in min,
        # time flags, reserved)
        data_bytes = (b'##FH', 0, 56, 2,
                      0, self['MD'],
                      int(time()*1E9),
                      0, 0, 1, b'\0' * 3)
        fid.write(pack('<4sI2Q3Q2hB3s', *data_bytes))


class CHBlock(dict):

    """ reads Channel Hierarchy block and saves in class dict
    """

    def __init__(self, fid, pointer):
        # block header
        self.update(_load_header(fid, pointer))
        # Channel hierarchy block
        (self['ch_ch_next'],
         self['ch_ch_first'],
         self['ch_tx_name'],
         self['ch_md_comment']) = unpack('<4Q', fid.read(32))
        n_links = self['link_count'] - 4
        self['ch_element'] = unpack('<{}Q'.format(n_links),
                                    fid.read(n_links * 8))
        (self['ch_element_count'],
         self['ch_type'],
         ch_reserved) = unpack('<IB3s', fid.read(8))
        if self['ch_md_comment']:  # comments exist
            self['Comment'] = {}
            comment = CommentBlock()
            comment.read_cm_ch(fid=fid, pointer=self['ch_md_comment'])
            self['Comment'].update(comment)
        if self['ch_tx_name']:  # text block containing name of hierarchy level
            self['ch_name_level'] = {}
            comment = CommentBlock()
            comment.read_tx(fid=fid, pointer=self['ch_tx_name'])
            self['ch_name_level'].update(comment)


class CommentBlock(dict):
    """ reads or writes Comment block and saves in class dict
    """

    def __init__(self):
        self = {}

    def read_tx(self, fid, pointer):
        """ reads TX block

        Parameters
        ----------
        fid:
            file identifier
        pointer: int
            position in file
        """
        if pointer > 0:
            fid.seek(pointer)
            # block header
            (self['id'],
             reserved,
             self['length'],
             self['link_count']) = _HeaderStruct.unpack(fid.read(24))
            self['Comment'] = fid.read(self['length'] - 24).rstrip(b'\x00').decode('UTF-8', 'ignore')

    def read_cm_header(self, fid, pointer):
        """ reads Comment block header

        Parameters
        ----------
        fid:
            file identifier
        pointer: int
            position in file
        """
        if pointer > 0:
            fid.seek(pointer)
            # block header
            (self['id'],
             reserved,
             self['length'],
             self['link_count']) = _HeaderStruct.unpack(fid.read(24))

    def read_xml(self, fid):
        """ reads Comment block xml and objectifies it

        Parameters
        ----------
        fid:
            file identifier

        Metadata block
        removes normal 0 at end
        """
        try:
            xml_string = fid.read(self['length'] - 24).rstrip(b'\x00')
            xml_tree = objectify.fromstring(xml_string)
        except:
            warn('xml metadata malformed')
            xml_tree = None
        return xml_tree

    def read_cm_hd(self, fid, pointer):
        """ reads Comment block from header block

         Parameters
         ----------
         fid:
             file identifier
         pointer: int
             position in file
         """
        if pointer > 0:
            self.read_cm_header(fid, pointer)
            if self['id'] in ('##MD', b'##MD'):
                xml_tree = self.read_xml(fid)
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
                    for t in tmp.e:
                        self[t.attrib.values()[0]] = t.text
                except AttributeError:
                    pass  # optional
            elif self['id'] in ('##TX', b'##TX'):
                self['Comment'] = fid.read(self['length'] - 24).rstrip(b'\x00').decode('UTF-8', 'ignore')

    def read_cm_fh(self, fid, pointer):
        """ reads Comment block from file history block

         Parameters
         ----------
         fid:
             file identifier
         pointer: int
             position in file
         """
        if pointer > 0:
            self.read_cm_header(fid, pointer)
            if self['id'] in ('##MD', b'##MD'):
                xml_tree = self.read_xml(fid)
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
            elif self['id'] in ('##TX', b'##TX'):
                self['Comment'] = fid.read(self['length'] - 24).rstrip(b'\x00').decode('UTF-8', 'ignore')

    def read_cm_ch(self, fid, pointer):
        """ reads Comment block from file channel hierarchy block

         Parameters
         ----------
         fid:
             file identifier
         pointer: int
             position in file
         """
        if pointer > 0:
            self.read_cm_header(fid, pointer)
            if self['id'] in ('##MD', b'##MD'):
                xml_tree = self.read_xml(fid)
                if xml_tree is not None:
                    try:
                        self['TX'] = xml_tree.TX.text
                    except AttributeError:
                        warn('Could not parse CH block TX tag')
                    try:
                        self['names'] = xml_tree.names.text
                    except AttributeError:
                        warn('Could not parse CH block names tag')
            elif self['id'] in ('##TX', b'##TX'):
                self['Comment'] = fid.read(self['length'] - 24).rstrip(b'\x00').decode('UTF-8', 'ignore')

    def read_cm_at(self, fid, pointer):
        """ reads Comment block from attachment block

         Parameters
         ----------
         fid:
             file identifier
         pointer: int
             position in file
         """
        if pointer > 0:
            self.read_cm_header(fid, pointer)
            if self['id'] in (b'##MD', '##MD'):
                xml_tree = self.read_xml(fid)
                try:
                    self['TX'] = xml_tree.TX.text
                except AttributeError:
                    warn('Could not parse AT block TX tag')
            elif self['id'] in (b'##TX', '##TX'):
                self['Comment'] = fid.read(self['length'] - 24).rstrip(b'\x00').decode('UTF-8', 'ignore')

    def read_cm_ev(self, fid, pointer):
        """ reads Comment block from event block

         Parameters
         ----------
         fid:
             file identifier
         pointer: int
             position in file
         """
        if pointer > 0:
            self.read_cm_header(fid, pointer)
            if self['id'] in ('##MD', b'##MD'):
                xml_tree = self.read_xml(fid)
                if xml_tree is not None:
                    try:
                        self['TX'] = EV_TX(xml_tree).text
                    except AttributeError:
                        warn('Could not parse EV block TX tag')
                    try:
                        self['pre_trigger_interval'] = EV_pre_trigger_interval(xml_tree).text
                    except AttributeError:
                        pass  # optional
                    try:
                        self['post_trigger_interval'] = EV_post_trigger_interval(xml_tree).text
                    except AttributeError:
                        pass  # optional
                    try:
                        self['formula'] = EV_formula(xml_tree).text
                    except AttributeError:
                        pass  # optional
                    try:
                        self['timeout'] = EV_timeout(xml_tree).text
                    except AttributeError:
                        pass  # optional
            elif self['id'] in ('##TX', b'##TX'):
                self['Comment'] = fid.read(self['length'] - 24).rstrip(b'\x00').decode('UTF-8', 'ignore')

    def read_cm_dg(self, fid, pointer):
        """ reads Comment block from data group block

         Parameters
         ----------
         fid:
             file identifier
         pointer: int
             position in file
         """
        if pointer > 0:
            self.read_cm_header(fid, pointer)
            if self['id'] in ('##MD', b'##MD'):
                xml_tree = self.read_xml(fid)
                try:
                    self['TX'] = xml_tree.TX.text
                except AttributeError:
                    warn('Could not parse AT block TX tag')
            elif self['id'] in ('##TX', b'##TX'):
                self['Comment'] = fid.read(self['length'] - 24).rstrip(b'\x00').decode('UTF-8', 'ignore')

    def read_cm_cg(self, fid, pointer):
        """ reads Comment block from channel group block

         Parameters
         ----------
         fid:
             file identifier
         pointer: int
             position in file
         """
        if pointer > 0:
            self.read_cm_header(fid, pointer)
            if self['id'] in ('##MD', b'##MD'):
                xml_tree = self.read_xml(fid)
                try:
                    self['TX'] = xml_tree.TX.text
                except AttributeError:
                    warn('Could not parse CG block TX tag')
                try:
                    self['names'] = xml_tree.names.text
                except AttributeError:
                    warn('Could not parse CG block names tag')
            elif self['id'] in ('##TX', b'##TX'):
                self['Comment'] = fid.read(self['length'] - 24).rstrip(b'\x00').decode('UTF-8', 'ignore')

    def read_cm_si(self, fid, pointer):
        """ reads Comment block from source information block

         Parameters
         ----------
         fid:
             file identifier
         pointer: int
             position in file
         """
        if pointer > 0:
            self.read_cm_header(fid, pointer)
            if self['id'] in ('##MD', b'##MD'):
                xml_tree = self.read_xml(fid)
                if xml_tree is not None:
                    try:
                        self['TX'] = SI_TX(xml_tree).text
                    except AttributeError:
                        warn('Could not parse SI block TX tag')
                    try:
                        self['names'] = SI_names(xml_tree).text
                    except AttributeError:
                        pass  # optional
                    try:
                        self['path'] = SI_path(xml_tree).text
                    except AttributeError:
                        pass  # optional
                    try:
                        self['bus'] = SI_bus(xml_tree).text
                    except AttributeError:
                        pass  # optional
                    try:
                        self['protocol'] = SI_protocol(xml_tree).text
                    except AttributeError:
                        pass  # optional
            elif self['id'] in ('##TX', b'##TX'):
                self['Comment'] = fid.read(self['length'] - 24).rstrip(b'\x00').decode('UTF-8', 'ignore')

    def read_cm_cn(self, fid, pointer, minimal=True):
        """ reads Comment block from channel block

        Parameters
        ----------
        fid:
            file identifier
        pointer: int
            position in file
        minimal: boolean
           flag to reduce metadata parsing
        """
        if pointer > 0:
            self.read_cm_header(fid, pointer)
            if self['id'] in ('##MD', b'##MD'):
                xml_tree = self.read_xml(fid)
                if xml_tree is not None:
                    try:
                        self['description'] = CN_TX(xml_tree).text
                    except AttributeError:
                        warn('Could not parse CN block TX tag')
                    try:
                        self['names'] = CN_names(xml_tree).text
                    except AttributeError:
                        pass  # optional
                    if minimal is False:
                        # not really used for the moment
                        try:
                            self['axis_monotony'] = CN_axis_monotony(xml_tree).text
                        except AttributeError:
                            pass  # optional
                        try:
                            self['raster'] = CN_raster(xml_tree).text
                        except AttributeError:
                            pass  # optional
                        try:
                            self['formula'] = CN_formula(xml_tree).text
                        except AttributeError:
                            pass  # optional
                        try:
                            self['linker_name'] = CN_linker_name(xml_tree).text
                        except AttributeError:
                            pass  # optional
                        try:
                            self['linker_address'] = CN_linker_address(xml_tree).text
                        except AttributeError:
                            pass  # optional
                        try:
                            self['address'] = CN_address(xml_tree).text
                        except AttributeError:
                            pass  # optional
            elif self['id'] in ('##TX', b'##TX'):
                self['name'] = fid.read(self['length'] - 24).rstrip(b'\x00').decode('UTF-8', 'ignore')

    def read_cm_cn_unit(self, fid, pointer):
        """ reads Comment block for channel unit

         Parameters
         ----------
         fid:
             file identifier
         pointer: int
             position in file
         """
        if pointer > 0:
            self.read_cm_header(fid, pointer)
            if self['id'] in ('##MD', b'##MD'):
                xml_tree = self.read_xml(fid)
                try:
                    self['unit'] = CN_unit_TX(xml_tree).text
                except AttributeError:
                    warn('Could not parse unit TX tag')
            elif self['id'] in ('##TX', b'##TX'):
                self['Comment'] = fid.read(self['length'] - 24).rstrip(b'\x00').decode('UTF-8', 'ignore')

    def read_cm_cc(self, fid, pointer):
        """ reads Comment block from channel conversion block

         Parameters
         ----------
         fid:
             file identifier
         pointer: int
             position in file
         """
        if pointer > 0:
            self.read_cm_header(fid, pointer)
            if self['id'] in ('##MD', b'##MD'):
                xml_tree = self.read_xml(fid)
                if xml_tree is not None:
                    try:
                        self['TX'] = CC_TX(xml_tree).text
                    except AttributeError:
                        warn('Could not parse CC block TX tag')
                    try:
                        self['names'] = CC_names(xml_tree).text
                    except AttributeError:
                        pass  # optional
                    try:
                        self['COMPU_METHOD'] = CC_COMPU_METHOD(xml_tree).text
                    except AttributeError:
                        pass  # optional
                    try:
                        self['formula'] = CC_formula(xml_tree).text
                    except AttributeError:
                        pass  # optional
            elif self['id'] in ('##TX', b'##TX'):
                self['Comment'] = fid.read(self['length'] - 24).rstrip(b'\x00').decode('UTF-8', 'ignore')

    def read_cm_cc_unit(self, fid, pointer):
        """ reads Comment block for channel conversion unit

         Parameters
         ----------
         fid:
             file identifier
         pointer: int
             position in file
         """
        if pointer > 0:
            self.read_cm_header(fid, pointer)
            if self['id'] in ('##MD', b'##MD'):
                xml_tree = self.read_xml(fid)
                try:
                    self['unit'] = CC_unit_TX(xml_tree).text
                except AttributeError:
                    warn('Could not parse unit TX tag')
            elif self['id'] in ('##TX', b'##TX'):
                self['Comment'] = fid.read(self['length'] - 24).rstrip(b'\x00').decode('UTF-8', 'ignore')

    def load(self, data, md_type):
        if md_type == 'TX':
            data = b''.join([data.encode('utf-8', 'replace'), b'\0'])
            block_id = b'##TX'
        else:
            register_namespace('', 'http://www.asam.net/mdf/v4')
            if md_type == 'HD':
                root = Element('HDcomment')
                root.set('xmlns', 'http://www.asam.net/mdf/v4')
                tx = SubElement(root, 'TX')
                tx.text = data['comment']
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
            elif md_type == 'CN':
                pass
            elif md_type == 'FH':
                root = Element('FHcomment')
                root.set('xmlns', 'http://www.asam.net/mdf/v4')
                tx = SubElement(root, 'TX')
                tx.text = data['comment']
                tool_id = SubElement(root, 'tool_id')
                tool_id.text = 'mdfreader'
                tool_vendor = SubElement(root, 'tool_vendor')
                tool_vendor.text = 'mdfreader is under GPL V3'
                tool_version = SubElement(root, 'tool_version')
                tool_version.text = '2.8'
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
        fid.write(_HeaderStruct.pack(*self['Header']))
        fid.write(self['data'])


class DGBlock(dict):

    """ reads Data Group block and saves in class dict
    """

    def __init__(self, fid=None, pointer=None):
        if fid is not None:
            self.read_dg(fid, pointer)

    def read_dg(self, fid, pointer):
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
         dg_reserved) = _DGStruct.unpack(fid.read(64))
        self['data_block_header'] = _load_header(fid, self['dg_data'])
        if self['dg_md_comment']:  # comments exist
            self['Comment'] = {}
            comment = CommentBlock()
            comment.read_cm_dg(fid=fid, pointer=self['dg_md_comment'])
            self['Comment'].update(comment)

    def write(self, fid):
        # write block header
        fid.seek(self['block_start'])
        # link section
        # (Next Data group pointer, first channel group block pointer,
        # data block pointer, comment block pointer,
        # no recordID, reserved)
        data_bytes = (b'##DG', 0, 64, 4,
                      self['DG'], self['CG'],
                      self['data'], 0, 0, b'\0' * 7)
        fid.write(pack('<4sI2Q4QB7s', *data_bytes))


class CGBlock(dict):

    """ reads Channel Group block and saves in class dict
    """

    def __init__(self, fid=None, pointer=None):
        if fid is not None:
            self.read_cg(fid, pointer)

    def read_cg(self, fid, pointer):
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
         self['cg_md_comment']) = _CGStruct1.unpack(fid.read(72))
        if self['link_count'] > 6:  # linked to another cg containing master
            self['cg_cg_master'] = unpack(_LINK, fid.read(8))[0]
        # data section
        (self['cg_record_id'],
         self['cg_cycle_count'],
         self['cg_flags'],
         self['cg_path_separator'],
         self['cg_reserved'],
         self['cg_data_bytes'],
         self['cg_invalid_bytes']) = _CGStruct2.unpack(fid.read(32))
        if self['cg_md_comment']:  # comments exist
            self['Comment'] = CommentBlock()
            self['Comment'].read_cm_cg(fid=fid, pointer=self['cg_md_comment'])
        if self['cg_tx_acq_name']:  # comments exist
            self['acq_name'] = {}
            comment = CommentBlock()
            comment.read_tx(fid=fid, pointer=self['cg_tx_acq_name'])
            self['acq_name'].update(comment)
        if self['cg_si_acq_source']:  # comments exist
            self['acq_source'] = {}
            si = SIBlock()
            si.read_si(fid=fid, pointer=self['cg_si_acq_source'])
            self['acq_source'].update(si)

    def write(self, fid):
        # write block header
        # fid.seek(self['block_start'])
        # link section
        # (Next channel group pointer, first channel block pointer,
        # acquisition name pointer, acquisition source pointer,
        # sample reduction pointer, comment pointer,
        # no recordID, cycle count, no flags,
        # no character specified for path separator,
        # reserved, number of bytes taken by data in record,
        # no invalid bytes)
        if 'cg_cg_master' in self and self['cg_cg_master'] is not None:
            data_bytes = (b'##CG', 0, self['length'], 7,
                          0, self['CN'], 0, 0, 0, 0, self['cg_cg_master'], 0,
                          self['cg_cycle_count'], 8, 0,
                          b'\0' * 4,
                          self['cg_data_bytes'], self['cg_inval_bytes'])
            fid.write(pack('<4sI2Q9Q2H4s2I', *data_bytes))
        else:
            data_bytes = (b'##CG', 0, self['length'], 6,
                          0, self['CN'], 0, 0, 0, 0, 0,
                          self['cg_cycle_count'], 0, 0,
                          b'\0' * 4,
                          self['cg_data_bytes'], self['cg_inval_bytes'])
            fid.write(pack('<4sI2Q8Q2H4s2I', *data_bytes))


class CNBlock(dict):

    """ reads Channel block and saves in class dict
    """

    def read_cn(self, **kargs):
        if kargs['pointer'] != 0 and kargs['pointer'] is not None:
            self['pointer'] = kargs['pointer']
            kargs['fid'].seek(kargs['pointer'])
            links = []
            (self['id'],
             reserved,
             self['length'],
             self['link_count'],
             self['cn_cn_next'],
             self['cn_composition'],
             self['cn_tx_name'],
             self['cn_si_source'],
             self['cn_cc_conversion'],
             self['cn_data'],
             self['cn_md_unit'],
             self['cn_md_comment']) = _CNStruct1.unpack(kargs['fid'].read(88))
            if self['link_count'] > 8:
                links = _mdf_block_read(kargs['fid'], _LINK, self['link_count'] - 8)
            # data section
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
             cn_limit_ext_max) = _CNStruct2.unpack(kargs['fid'].read(72))
            if self['cn_attachment_count'] > 0:
                self['attachment'] = {}
                if self['cn_attachment_count'] > 1:
                    self['cn_at_reference'] = links[:self['cn_attachment_count']]
                    for at in range(self['cn_attachment_count']):
                        self['attachment'][at] = \
                            ATBlock(kargs['fid'], self['cn_at_reference'][at])
                else:
                    self['cn_at_reference'] = links
                    self['attachment'][0] = \
                        ATBlock(kargs['fid'], self['cn_at_reference'])
            if self['link_count'] > (8 + self['cn_attachment_count']):
                self['cn_default_x'] = links[self['cn_attachment_count']:]
            else:
                self['cn_default_x'] = None
            if self['cn_md_comment']:  # comments exist
                self['Comment'] = CommentBlock()
                self['Comment'].read_cm_cn(fid=kargs['fid'], pointer=self['cn_md_comment'], minimal=True)
            if self['cn_md_unit']:  # comments exist
                self['unit'] = CommentBlock()
                self['unit'].read_cm_cn_unit(fid=kargs['fid'], pointer=self['cn_md_unit'])
                if self['cn_sync_type'] and not self['unit']:
                    # no units but already known by spec
                    if self['cn_sync_type'] == 1:
                        self['unit'] = 's'
                    elif self['cn_sync_type'] == 2:
                        self['unit'] = 'rad'
                    elif self['cn_sync_type'] == 3:
                        self['unit'] = 'm'
            if self['cn_tx_name']:  # comments exist
                comment = CommentBlock()
                comment.read_tx(fid=kargs['fid'], pointer=self['cn_tx_name'])
                self['name'] = comment['Comment']

    def write(self, fid):
        # write block header
        # fid.seek(self['block_start'])
        # no attachment and default X
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
        data_bytes = (b'##CN', 0, 160, 8,
                      self['CN'], self['Composition'], self['TX'], 0, 0, 0, self['Unit'], self['Comment'],
                      self['cn_type'], self['cn_sync_type'],
                      self['cn_data_type'], self['cn_bit_offset'],
                      self['cn_byte_offset'], self['cn_bit_count'],
                      self['cn_flags'], 0, 0, 0, 0,
                      self['cn_val_range_min'], self['cn_val_range_max'],
                      0, 0, 0, 0)
        fid.write(pack('<4sI2Q8Q4B4I2BH6d', *data_bytes))


class CCBlock(dict):

    """ reads Channel Conversion block and saves in class dict
    """

    def read_cc(self, fid, pointer):
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
             self['cc_cc_inverse']) = _CCStruct1.unpack(fid.read(56))
            if self['link_count'] - 4 > 0:  # can be no links for cc_ref
                self['cc_ref'] = _mdf_block_read(fid, _LINK, self['link_count'] - 4)
            # data section
            (self['cc_type'],
             self['cc_precision'],
             self['cc_flags'],
             self['cc_ref_count'],
             self['cc_val_count'],
             self['cc_phy_range_min'],
             self['cc_phy_range_max']) = _CCStruct2.unpack(fid.read(24))
            if self['cc_val_count']:
                self['cc_val'] = _mdf_block_read(fid, _REAL, self['cc_val_count'])
            if self['cc_type'] == 3:  # reads Algebraic formula
                pointer = self['cc_ref']
                self['cc_ref'] = {}
                cc = CommentBlock()
                cc.read_tx(fid=fid, pointer=pointer)
                self['cc_ref'].update(cc)
            elif self['cc_type']in (7, 8, 9, 10, 11):  # text list
                self['cc_ref'] = list(self['cc_ref'])
                for i in range(self['cc_ref_count']):
                    fid.seek(self['cc_ref'][i])
                    # find if TX/MD or another CCBlock
                    identifier = unpack('4s', fid.read(4))[0]
                    # for algebraic formulae
                    if identifier in ('##TX', '##MD', b'##TX', b'##MD'):
                        temp = CommentBlock()
                        temp.read_tx(fid=fid, pointer=self['cc_ref'][i])
                        self['cc_ref'][i] = temp['Comment']
                    elif identifier in ('##CC', b'##CC'):  # for table conversion
                        # much more complicated nesting conversions !!!
                        cc = CCBlock()
                        cc.read_cc(fid, self['cc_ref'][i])
                        self['cc_ref'][i] = cc
            if self['cc_md_comment']:  # comments exist
                self['Comment'] = CommentBlock()
                self['Comment'].read_cm_cc(fid=fid, pointer=self['cc_md_comment'])
            if self['cc_md_unit']:  # comments exist
                self['unit'] = CommentBlock()
                self['unit'].read_cm_cc_unit(fid=fid, pointer=self['cc_md_unit'])
            if self['cc_tx_name']:  # comments exist
                self['name'] = CommentBlock()
                self['name'].read_tx(fid=fid, pointer=self['cc_tx_name'])
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
             self['link_count']) = _HeaderStruct.unpack(fid.read(24))
            self['pointer'] = pointer
            # reads data section
            fid.seek(pointer + 24 + self['link_count'] * 8)
            (self['ca_type'],
             self['ca_storage'],
             self['ca_ndim'],
             self['ca_flags'],
             self['ca_byte_offset_base'],
             self['ca_invalid_bit_pos_base']) = unpack('2BHIiI', fid.read(16))
            self['ca_dim_size'] = _mdf_block_read(fid, _UINT64, self['ca_ndim'])
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
                    _mdf_block_read(fid, _REAL, self['SNd'])
            if self['ca_storage'] >= 1:
                self['ca_cycle_count'] = \
                    _mdf_block_read(fid, _UINT64, self['PNd'])
            # Channel Conversion block : Links
            fid.seek(pointer + 24)
            # point to CN for array of structures or CA for array of array
            self['ca_composition'] = _mdf_block_read(fid, _LINK, 1)
            if self['ca_storage'] == 2:
                self['ca_data'] = _mdf_block_read(fid, _LINK, self['PNd'])
            if 1 << 0 & self['ca_flags']:  # bit 0
                self['ca_dynamic_size'] = \
                    _mdf_block_read(fid, _LINK, self['ca_ndim'] * 3)
            if 1 << 1 & self['ca_flags']:  # bit 1
                self['ca_input_quantity'] = \
                    _mdf_block_read(fid, _LINK, self['ca_ndim'] * 3)
            if 1 << 2 & self['ca_flags']:  # bit 2
                self['ca_output_quantity'] = \
                    _mdf_block_read(fid, _LINK, 3)
            if 1 << 3 & self['ca_flags']:  # bit 3
                self['ca_comparison_quantity'] = _mdf_block_read(fid, _LINK, 3)
            if 1 << 4 & self['ca_flags']:  # bit 4
                self['ca_cc_axis_conversion'] = \
                    _mdf_block_read(fid, _LINK, self['ca_ndim'])
            if 1 << 4 & self['ca_flags'] and not 1 << 5 & self['ca_flags']:
                # bit 4 and 5
                self['ca_axis'] = _mdf_block_read(fid, _LINK, self['ca_ndim'] * 3)
            # nested arrays
            if self['ca_composition']:
                self['CABlock'] = CABlock()
                self['CABlock'].read(fid, self['ca_composition'])

    def load(self, byte_offset_base):
        self['block_length'] = 48 + 8 * self['ndim']
        self['byte_offset_base'] = byte_offset_base

    def write(self, fid):
        # default CN template
        data_bytes = (b'##CA', 0, 48 + 8 * self['ndim'], 1,
                      0,
                      0, 0, self['ndim'], 0, self['byte_offset_base'], 0)
        fid.write(pack('<4sI2Q1Q2BHIiI', *data_bytes))
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
                comment.read_cm_at(fid=fid, pointer=self['at_md_comment'])
                self['Comment'].update(comment)
            if self['at_tx_filename']:  # file name
                temp = CommentBlock()
                temp.read_tx(fid, self['at_tx_filename'])
                self['at_tx_filename'] = temp['Comment']
            if self['at_tx_mimetype']:  # file name mime type
                temp = CommentBlock()
                temp.read_tx(fid, self['at_tx_mimetype'])
                self['at_tx_mimetype'] = temp['Comment']


class EVBlock(dict):

    """ reads Event block and saves in class dict
    """

    def __init__(self, fid, pointer):
        # block header
        if pointer != 0 and pointer is not None:
            fid.seek(pointer)
            links = []
            (self['id'],
             reserved,
             self['length'],
             self['link_count']) = _HeaderStruct.unpack(fid.read(24))
            # link section
            (self['ev_ev_next'],
             self['ev_ev_parent'],
             self['ev_ev_range'],
             self['ev_tx_name'],
             self['ev_md_comment']) = unpack('<5Q', fid.read(40))
            if self['link_count'] > 5:
                links = _mdf_block_read(fid, _LINK, self['link_count'] - 5)
            # data section
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
            # variable links
            if self['ev_scope_count'] > 0:
                if self['ev_scope_count'] > 1:
                    self['ev_scope'] = links[:self['ev_scope_count']]
                else:
                    self['ev_scope'] = links
            if self['ev_attachment_count'] > 0:
                if self['ev_attachment_count'] > 1:
                    self['ev_at_reference'] = links[self['ev_scope_count']:
                                                    self['ev_scope_count'] + self['ev_attachment_count']]
                else:
                    self['ev_at_reference'] = links
            # post treatment
            try:
                self['ev_cause'] = ev_cause[self['ev_cause']]
            except KeyError:
                warn('unexpected ev cause')
            if self['ev_md_comment']:  # comments exist
                self['Comment'] = CommentBlock()
                self['Comment'].read_cm_ev(fid=fid, pointer=self['ev_md_comment'])
            if self['ev_tx_name']:  # comments exist
                temp = CommentBlock()
                temp.read_tx(fid=fid, pointer=self['ev_tx_name'])
                self['ev_tx_name'] = temp['Comment']
            if self['ev_flags'] & (1 << 1):  # Bit 1 set
                self['ev_tx_group_name'] = links[-1:]


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
             sr_reserved) = _SRStruct.unpack(fid.read(64))


class SIBlock(dict):

    """ reads Source Information block and saves in class dict
    """

    def __init__(self):
        self = {}

    def read_si(self, fid, pointer):
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
             si_reserved) = _SIStruct.unpack(fid.read(56))
            try:
                self['si_type'] = si_type[self['si_type']]
            except KeyError:
                warn('unexpected SI type')
            try:
                self['si_bus_type'] = si_bus_type[self['si_bus_type']]
            except KeyError:
                warn('unexpected SI bus type')
            # post treatment
            self['source_name'] = CommentBlock()
            self['source_name'].read_tx(fid=fid, pointer=self['si_tx_name'])
            self['source_path'] = CommentBlock()
            self['source_path'].read_tx(fid=fid, pointer=self['si_tx_path'])
            self['comment'] = CommentBlock()
            self['comment'].read_cm_si(fid=fid, pointer=self['si_md_comment'])


class DTBlock(dict):
    def load(self, record_byte_offset, nRecords, pointer):
        self['datablocks_length'] = 24 + record_byte_offset * nRecords
        self['pointer'] = pointer
        self['end_position'] = _calculate_block_start(self['pointer'] + self['datablocks_length'])

    def write(self, fid, data):
        fid.seek(self['pointer'])
        fid.write(_HeaderStruct.pack(b'##DT', 0, self['datablocks_length'], 0))
        # dumps data
        fid.write(data)
        return self['end_position']


class DVBlock(dict):
    def load(self, record_byte_offset, nRecords, pointer):
        self['datablocks_length'] = 24 + record_byte_offset * nRecords
        self['pointer'] = pointer
        self['end_position'] = _calculate_block_start(self['pointer'] + self['datablocks_length'])

    def write(self, fid, data):
        fid.seek(self['pointer'])
        fid.write(_HeaderStruct.pack(b'##DV', 0, self['datablocks_length'], 0))
        # dumps data
        fid.write(data)
        return self['end_position']


class DIBlock(dict):
    def load(self, invalid_bytes, nRecords, pointer):
        self['datablocks_length'] = 24 + invalid_bytes * nRecords
        self['pointer'] = pointer
        self['end_position'] = _calculate_block_start(self['pointer'] + self['datablocks_length'])

    def write(self, fid, data):
        fid.seek(self['pointer'])
        fid.write(_HeaderStruct.pack(b'##DI', 0, self['datablocks_length'], 0))
        # dumps data
        fid.write(data)
        return self['end_position']


class LDBlock(dict):

    """ reads List Data block
    """

    def read_ld(self, fid, link_count):
        # block header is already read
        self['next'] = unpack('<Q', fid.read(8))[0]
        self['list_data'] = {}
        self['list_data'][0] = list(unpack('<{}Q'.format(link_count - 1),
                                        fid.read(8 * (link_count - 1))))
        (self['flags'],
         self['count']) = unpack('<2I', fid.read(8))
        if self['flags']:  # flags existing
            if self['flags'] & (1 << 31):  # invalid data present
                self['inval_data'] = {}
                self['inval_data'][0] = [self['list_data'][0].pop(i) for i in range((link_count - 1)//2, link_count - 1)]
            if self['flags'] & 0b1:  # equal length data list
                self['equal_sample_count'] = unpack('<Q', fid.read(8))[0]
            else:  # data list defined by byte offset
                self['sample_offset'] = unpack('<{}Q'.format(self['count']),
                                               fid.read(8 * self['count']))
            if self['flags'] & 0b10:  # time values
                self['time_values'] = unpack('<{}Q'.format(self['count']),
                                             fid.read(8 * self['count']))
            if self['flags'] & 0b100:  # angle values
                self['angle_values'] = unpack('<{}Q'.format(self['count']),
                                              fid.read(8 * self['count']))
            if self['flags'] & 0b1000:  # distance values
                self['distance_values'] = unpack('<{}Q'.format(self['count']),
                                                 fid.read(8 * self['count']))

    def load(self, record_byte_offset, n_records, position, invalid_bytes=0, column_oriented_flag=False):
        self['record_length'] = record_byte_offset
        self['invalid_bytes'] = invalid_bytes
        self['n_records'] = n_records
        self['block_start'] = position
        self['column_oriented_flag'] = column_oriented_flag
        self['block_length'] = 40

        # calculate data chunks
        n_chunks = (record_byte_offset * n_records) // chunk_size_writing + 1
        chunk_length = (record_byte_offset * n_records) // n_chunks
        n_record_chunk = chunk_length // record_byte_offset
        self['chunks'] = [(n_record_chunk, record_byte_offset * n_record_chunk)] * n_chunks
        n_record_chunk = n_records - n_record_chunk * n_chunks
        if n_record_chunk > 0:
            self['chunks'].append((n_record_chunk, record_byte_offset * n_record_chunk))

    def write(self, fid, data, invalid_data=None, compression_flag=False):
        number_ld = len(self['chunks'])
        ld_flags = 0
        if invalid_data is not None:
            number_invalid_ld = number_ld
            dl_invalid_data = zeros(shape=len(self['chunks']), dtype='<u8')
            ld_flags = 1 << 31
        else:
            number_invalid_ld = 0
        self['block_length'] = 40 + 16 * number_ld + 8 * number_invalid_ld
        ld_offset = zeros(shape=number_ld, dtype='<u8')
        if number_ld > 1:
            for counter in range(1, number_ld):
                (n_record_chunk, chunk_size) = self['chunks'][counter]
                ld_offset[counter] = ld_offset[counter - 1] + chunk_size
        data_bytes = (b'##LD', 0, self['block_length'], number_ld + number_invalid_ld + 1, 0)
        fid.write(pack('<4sI3Q', *data_bytes))
        fid.write(pack('{0}Q'.format(number_ld), *zeros(shape=number_ld, dtype='<u8')))
        if invalid_data is not None:
            fid.write(pack('{0}Q'.format(number_ld), *zeros(shape=number_ld, dtype='<u8')))
        fid.write(pack('<2I', ld_flags, number_ld))
        fid.write(pack('{0}Q'.format(number_ld), *ld_offset))

        # writing data blocks
        pointer = self['block_start'] + self['block_length']
        dl_data = zeros(shape=len(self['chunks']), dtype='<u8')
        data_pointer = 0
        data_invalid_pointer = 0
        if compression_flag:
            if self['column_oriented_flag']:
                dz_zip_type = 0  # for a column oriented channel, no meaning to have transposition
            else:
                # uses not equal length blocks, transposed compressed data
                dz_zip_type = 1

            for counter, (n_record_chunk, chunk_size) in enumerate(self['chunks']):
                position, dl_data = self.write_DZ(fid, pointer, data, dl_data, counter, data_pointer,
                                                  self['record_length'], chunk_size, dz_zip_type, b'DV')
                if position is not None:
                    pointer = position
                else:  # no compression, not enough data
                    pointer, dl_data = self.write_DIV(fid, pointer, DVBlock(), data, dl_data,
                                                      counter, data_pointer, self['record_length'],
                                                      chunk_size, n_record_chunk)
                if invalid_data is not None:
                    position, dl_invalid_data = self.write_DZ(fid, pointer, invalid_data, dl_invalid_data,
                                                              counter, data_invalid_pointer, self['invalid_bytes'],
                                                              n_record_chunk, dz_zip_type, b'DI')
                    if position is not None:
                        pointer = position
                    else:  # no compression
                        pointer, dl_invalid_data = self.write_DIV(fid, pointer, DIBlock(), invalid_data,
                                                                  dl_invalid_data, counter,
                                                                  data_invalid_pointer, self['invalid_bytes'],
                                                                  n_record_chunk * self['invalid_bytes'],
                                                                  n_record_chunk)
                data_pointer += chunk_size
                data_invalid_pointer += n_record_chunk
        else:
            for counter, (n_record_chunk, chunk_size) in enumerate(self['chunks']):
                pointer, dl_data = self.write_DIV(fid, pointer, DVBlock(), data, dl_data,
                                                  counter, data_pointer, self['record_length'],
                                                  chunk_size, n_record_chunk)
                if invalid_data is not None:
                    pointer, dl_invalid_data = self.write_DIV(fid, pointer, DIBlock(), invalid_data,
                                                              dl_invalid_data, counter, data_invalid_pointer,
                                                              self['invalid_bytes'],
                                                              n_record_chunk * self['invalid_bytes'],
                                                              n_record_chunk)
                data_pointer += chunk_size
                data_invalid_pointer += n_record_chunk

        # writes links to all Blocks
        # write dl_data
        fid.seek(self['block_start'] + 32)
        fid.write(dl_data.tobytes())
        if invalid_data is not None:
            fid.write(dl_invalid_data.tobytes())
        fid.seek(pointer)
        return _calculate_block_start(pointer)

    def write_DIV(self, fid, pointer, block, data, dl_data, counter, data_pointer,
                  record_length, chunk_size, n_records):
        block.load(record_length, n_records, pointer)
        dl_data[counter] = block['pointer']
        position = block.write(fid, data[data_pointer: data_pointer + chunk_size])
        return position, dl_data

    def write_DZ(self, fid, pointer, data, dl_data, counter, data_pointer, record_length, chunk_size,
                 dz_zip_type, dz_org_block_type):
        DZ = DZBlock()
        DZ['block_start'] = _calculate_block_start(pointer)
        DZ['dz_org_block_type'] = dz_org_block_type
        DZ['dz_zip_type'] = dz_zip_type
        dl_data[counter] = DZ['block_start']
        position = DZ.write(fid, data[data_pointer: data_pointer + chunk_size], record_length)
        return position, dl_data


class DLBlock(dict):

    """ reads List Data block
    """

    def read_dl(self, fid, link_count):
        # block header is already read
        self['next'] = unpack('<Q', fid.read(8))[0]
        self['list_data'] = {}
        self['list_data'][0] = unpack('<{}Q'.format(link_count - 1),
                                    fid.read(8 * (link_count - 1)))
        (self['flags'],
         dl_reserved,
         self['count']) = unpack('<B3sI', fid.read(8))
        if self['flags'] & 0b1:  # equal length data list
            self['equal_length'] = unpack('<Q', fid.read(8))[0]
        else:  # data list defined by byte offset
            self['offset'] = unpack('<{}Q'.format(self['count']),
                                    fid.read(8 * self['count']))
        if self['flags'] & 0b10:  # time values
            self['time_values'] = unpack('<{}Q'.format(self['count']),
                                         fid.read(8 * self['count']))
        if self['flags'] & 0b100:  # angle values
            self['angle_values'] = unpack('<{}Q'.format(self['count']),
                                          fid.read(8 * self['count']))
        if self['flags'] & 0b1000:  # distance values
            self['distance_values'] = unpack('<{}Q'.format(self['count']),
                                             fid.read(8 * self['count']))

    def write(self, fid, chunks):
        number_dl = len(chunks)
        self['block_length'] = 40 + 16 * number_dl
        dl_offset = zeros(shape=number_dl, dtype='<u8')
        if number_dl > 1:
            for counter in range(1, number_dl):
                (n_record_chunk, chunk_size) = chunks[counter]
                dl_offset[counter] = dl_offset[counter - 1] + chunk_size
        data_bytes = (b'##DL', 0, self['block_length'], number_dl + 1, 0)
        fid.write(pack('<4sI3Q', *data_bytes))
        fid.write(pack('{0}Q'.format(number_dl), *zeros(shape=number_dl, dtype='<u8')))
        fid.write(pack('<2I', 0, number_dl))
        fid.write(pack('{0}Q'.format(number_dl), *dl_offset))


class DZBlock(dict):

    """ reads Data List block
    """

    def read_dz(self, fid):
        # block header is already read
        (self['dz_org_block_type'],
         self['dz_zip_type'],
         dz_reserved,
         self['dz_zip_parameter'],
         self['dz_org_data_length'],
         self['dz_data_length']) = _DZStruct.unpack(fid.read(24))

    @staticmethod
    def decompress_data_block(block, zip_type, zip_parameter, org_data_length):
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

        # decompress data
        block = decompress(block)
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
        if self['dz_zip_type']:
            M = org_data_length // record_length
            temp = array(memoryview(data[:M * record_length]))
            tail = array(memoryview(data[M * record_length:]))
            temp = temp.reshape(M, record_length).T.ravel()
            if len(tail) > 0:
                temp = append(temp, tail)
            # compress transposed data
            compressed_data = compress(temp.tostring())
        else:
            compressed_data = compress(data)
            record_length = 0
        dz_data_length = len(compressed_data)
        if org_data_length > dz_data_length + 24:
            self['DZBlock_length'] = 48 + dz_data_length
            # writes header
            fid.write(_HeaderStruct.pack(b'##DZ', 0, self['DZBlock_length'], 0))
            data_bytes = (self['dz_org_block_type'], self['dz_zip_type'], 0,
                          record_length, org_data_length, dz_data_length)
            fid.write(_DZStruct.pack(*data_bytes))
            # writes compressed data
            fid.write(compressed_data)
            return self['block_start'] + self['DZBlock_length']
        else:
            # not enough data to compress
            return None


class HLBlock(dict):

    """ reads Header List block
    """

    def read_hl(self, fid):
        (self['hl_dl_first'],
         self['hl_flags'],
         self['hl_zip_type'],
         hl_reserved) = _HLStruct.unpack(fid.read(16))

    def load(self, record_byte_offset, n_records, position):
        self['record_length'] = record_byte_offset
        self['block_start'] = position
        self['block_length'] = 40
        # calculate data chunks
        n_chunks = (record_byte_offset * n_records) // chunk_size_writing + 1
        chunk_length = (record_byte_offset * n_records) // n_chunks
        n_record_chunk = chunk_length // record_byte_offset
        self['chunks'] = [(n_record_chunk, record_byte_offset * n_record_chunk)] * n_chunks
        n_record_chunk = n_records - n_record_chunk * n_chunks
        if n_record_chunk > 0:
            self['chunks'].append((n_record_chunk, record_byte_offset * n_record_chunk))

    def write(self, fid, data):
        fid.write(_HeaderStruct.pack(b'##HL', 0, self['block_length'], 1))
        # uses not equal length blocks, transposed compressed data
        DL = DLBlock()
        dz_org_block_type = b'DT'
        dz_zip_type = 1
        fid.write(_HLStruct.pack(self['block_start'] + self['block_length'], 0, dz_zip_type, b'\x00' * 5))
        pointer = self['block_start'] + self['block_length']
        DL['block_start'] = pointer
        DL.write(fid, self['chunks'])
        pointer += DL['block_length']
        dl_data = zeros(shape=len(self['chunks']), dtype='<u8')
        data_pointer = 0
        for counter, (n_record_chunk, chunk_size) in enumerate(self['chunks']):
            DZ = DZBlock()
            DZ['block_start'] = _calculate_block_start(pointer)
            DZ['dz_org_block_type'] = dz_org_block_type
            DZ['dz_zip_type'] = dz_zip_type
            dl_data[counter] = DZ['block_start']
            position = DZ.write(fid, data[data_pointer: data_pointer + chunk_size], self['record_length'])
            if position is not None:
                pointer = position
            else:  # not enough data to be compressed, back to normal DTBlock
                DT = DTBlock()
                DT.load(self['record_length'], n_record_chunk, pointer)
                dl_data[counter] = DT['pointer']
                pointer = DT.write(fid, data[data_pointer: data_pointer + chunk_size])
            data_pointer += chunk_size
        # writes links to all DZBlocks
        # write dl_data
        fid.seek(DL['block_start'] + 32)
        fid.write(dl_data.tobytes())
        fid.seek(pointer)
        return _calculate_block_start(pointer)


class Info4(dict):
    __slots__ = ['fileName', 'fid', 'filterChannelNames', 'zipfile', 'multi_sources']
    """ information block parser fo MDF file version 4.x

    Attributes
    --------------
    file_name : str
        name of file
    fid
        file identifier
    zipfile
        flag to indicate the mdf4 is packaged in a zip
        
    multisources : bool
    0 will filter without source message information
    1 will filter with source message information

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

    def __init__(self, file_name=None, fid=None, filter_channel_names=False, minimal=0, multi_sources=False):
        """ info4 class constructor

        Parameters
        ----------------
        file_name : str
            file name
        fid : float
            file identifier
        filter_channel_names : bool, optional
            flag to filter long channel names including module names separated by a '.'
        minimal : int
            0 will load every metadata
            1 will load DG, CG, CN and CC (for noDataLoading)
            2 will load only DG (for normal reading)

        multisources : bool, default false
            0 will filter without source message information
            1 will filter with source message information

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
        self['masters'] = dict()  # channels grouped by master
        self.filterChannelNames = filter_channel_names
        self.fileName = file_name
        self.fid = None
        if fid is None and file_name is not None:
            # Open file
            (self.fid, self.fileName, self.zipfile) = _open_mdf(self.fileName)
        if self.fileName is not None and fid is None:
            self.read_info(self.fid, minimal, multi_sources)
            # Close the file
            self.fid.close()
            if self.zipfile:  # temporary uncompressed file, to be removed
                remove(file_name)
        elif self.fileName is None and fid is not None:
            # called by mdfreader.mdfinfo
            self.read_info(fid, minimal)

    def read_info(self, fid, minimal, multi_sources=False):
        """ read all file blocks except data

        Parameters
        ----------------
        fid : identifier
            file identifier
        minimal: flag
            to activate minimum content reading for raw data fetching
        multi_sources: flag
            to read same signal with same name from different source message
        """
        # reads IDBlock
        self['ID'].update(IDBlock(fid))

        # reads Header HDBlock
        self['HD'].update(HDBlock(fid))

        # reads File History blocks, always exists
        self['FH'] = list()
        self['FH'].append(FHBlock(fid, self['HD']['hd_fh_first']))
        while self['FH'][-1]['fh_fh_next']:
            self['FH'].append(FHBlock(fid, self['FH'][-1]['fh_fh_next']))

        # reads Channel Hierarchy blocks
        if self['HD']['hd_ch_first'] and not minimal:
            self['CH'] = self.read_ch_block(fid, self['HD']['hd_ch_first'])

        # reads Attachment block
        if self['HD']['hd_at_first'] and not minimal == 1:
            self['AT'] = OrderedDict()
            pointer = self['HD']['hd_at_first']
            while pointer:
                self['AT'][pointer] = ATBlock(fid, pointer)
                pointer = self['AT'][pointer]['at_at_next']

        # reads Event Block
        if self['HD']['hd_ev_first'] and not minimal == 1:
            self['EV'] = self.read_ev_block(fid, self['HD']['hd_ev_first'])

        # reads Data Group Blocks and recursively the other related blocks
        self.read_dg_block(fid, False, minimal, multi_sources)

    def read_dg_block(self, fid, channel_name_list=False, minimal=0, multi_sources=False):
        """reads Data Group Blocks

        Parameters
        ----------------
        fid : float
            file identifier
        channel_name_list : bool
            Flag to reads only channel blocks for listChannels4 method
        minimal: flag
            to activate minimum content reading for raw data fetching
        multi_sources: flag
            to read same signal with same name from different source message
        """
        self['ChannelNamesByDG'] = {}
        if self['HD']['hd_dg_first']:
            dg = 0
            self['DG'][dg] = {}
            self['DG'][dg].update(DGBlock(fid, self['HD']['hd_dg_first']))
            self['ChannelNamesByDG'][dg] = set()
            if minimal < 2:
                # reads Channel Group blocks
                self.read_cg_blocks(fid, dg, channel_name_list, minimal, multi_sources)
            while self['DG'][dg]['dg_dg_next']:
                dg += 1
                self['DG'][dg] = {}
                self['DG'][dg].update(DGBlock(fid, self['DG'][dg - 1]['dg_dg_next']))
                self['ChannelNamesByDG'][dg] = set()
                if minimal < 2:
                    # reads Channel Group blocks
                    self.read_cg_blocks(fid, dg, channel_name_list, minimal, multi_sources)

    def read_cg_blocks(self, fid, dg, channel_name_list=False, minimal=0, multi_sources=False):
        """reads Channel Group blocks linked to same Data Block dg

        Parameters
        ----------------
        fid : float
            file identifier
        dg : int
            data group number
        channel_name_list : bool
            Flag to reads only channel blocks for listChannels4 method
        minimal: falg
            to activate minimum content reading for raw data fetching
        multi_sources: flag
            to read same signal with same name from different source message
        """
        if self['DG'][dg]['dg_cg_first']:
            cg = 0
            self['CN'][dg] = {}
            self['CN'][dg][cg] = dict()
            self['CC'][dg] = {}
            self['CC'][dg][cg] = dict()
            self['CG'][dg] = {}

            vlsd_cg_block = []

            vlsd_cg_block = self.read_cg_block(fid, dg, cg, self['DG'][dg]['dg_cg_first'],
                                               vlsd_cg_block, channel_name_list=False, minimal=0, multi_sources= multi_sources)

            if self['CN'][dg][cg] and self['CG'][dg][cg]['unique_channel_in_CG'] and \
                    not self['CG'][dg][cg]['cg_cg_next']:
                # unique channel in data group
                self['DG'][dg]['unique_channel_in_DG'] = True
            else:
                self['DG'][dg]['unique_channel_in_DG'] = False

            while self['CG'][dg][cg]['cg_cg_next']:
                cg += 1
                self['CG'][dg][cg] = {}
                self['CN'][dg][cg] = dict()
                self['CC'][dg][cg] = dict()
                vlsd_cg_block = self.read_cg_block(fid, dg, cg, self['CG'][dg][cg - 1]['cg_cg_next'],
                                                   vlsd_cg_block, channel_name_list=False, minimal=0, multi_sources= multi_sources)

            if vlsd_cg_block and 'VLSD_CG' not in self:  # VLSD CG Block exiting
                self['VLSD_CG'] = {}

            # Matching VLSD CGBlock with corresponding channel
            for VLSDcg in vlsd_cg_block:
                vlsd_cg_block_address = self['CG'][dg][VLSDcg]['pointer']
                for cg in self['CG'][dg]:
                    if cg not in vlsd_cg_block:
                        for cn in self['CN'][dg][cg]:
                            if vlsd_cg_block_address == self['CN'][dg][cg][cn]['cn_data']:
                                # found matching channel with VLSD CGBlock
                                self['VLSD_CG'][self['CG'][dg][VLSDcg]['cg_record_id']] = {'cg_cn': (cg, cn)}
                                break

    def read_cg_block(self, fid, dg, cg, pointer, vlsd_cg_block, channel_name_list=False, minimal=0, multi_sources=False):
        """reads one Channel Group block

        Parameters
        ----------------
        fid : float
            file identifier
        dg : int
            data group number
        cg : int
            channel group number
        channel_name_list : bool
            Flag to reads only channel blocks for listChannels4 method
        minimal: flag
            to activate minimum content reading for raw data fetching
        multi_sources: flag
            to read same signal with same name from different source message

        Returns
        -----------
        vlsd_cg_block : boolean
        """
        self['CG'][dg][cg] = {}
        self['CG'][dg][cg].update(CGBlock(fid, pointer))
        if not channel_name_list and minimal:
            # reads Source Information Block
            temp = SIBlock()
            temp.read_si(fid, self['CG'][dg][cg]['cg_si_acq_source'])
            if temp:
                self['CG'][dg][cg]['SI'] = dict()
                self['CG'][dg][cg]['SI'].update(temp)

            # reads Sample Reduction Block
            self['CG'][dg][cg]['SR'] = self.read_sr_block(fid, self['CG'][dg][cg]['cg_sr_first'])

        if not self['CG'][dg][cg]['cg_flags'] & 0b1:  # if not a VLSD channel group
            # reads Channel Block
            vlsd = self.read_cn_blocks(fid, dg, cg, channel_name_list, minimal, multi_sources)
            if vlsd:
                # VLSD needs to rename and append records but with python 2.x impossible,
                # convert name to compatible python identifier
                for cn in self['CN'][dg][cg]:
                    self['CN'][dg][cg][cn]['name'] = _convert_name(
                        self['CN'][dg][cg][cn]['name'].encode('latin-1', ' ignore'))
        else:
            vlsd_cg_block.append(cg)

        return vlsd_cg_block

    def read_cn_blocks(self, fid, dg, cg, channel_name_list=False, minimal=0, multi_sources=False):
        """reads Channel blocks link to CG Block

        Parameters
        ----------------
        fid : float
            file identifier
        dg : int
            data group number
        cg : int
            channel group number in data group
        channel_name_list : bool
            Flag to reads only channel blocks for listChannels4 method
        minimal: flag
            to activate minimum content reading for raw data fetching
        multi_sources: flag
            to read same signal with same name from different source message

        Returns
        -----------
        vlsd : boolean
        """

        vlsd = False
        mlsd_channels = []
        cn, mlsd_channels, vlsd = self.read_cn_block(fid, self['CG'][dg][cg]['cg_cn_first'],
                                                     dg, cg, mlsd_channels, vlsd, minimal, channel_name_list, multi_sources)
        if not self['CN'][dg][cg][cn]['cn_cn_next']:  # only one channel in CGBlock
            self['CG'][dg][cg]['unique_channel_in_CG'] = True
        else:
            self['CG'][dg][cg]['unique_channel_in_CG'] = False

        while self['CN'][dg][cg][cn]['cn_cn_next']:
            cn, mlsd_channels, vlsd = self.read_cn_block(fid, self['CN'][dg][cg][cn]['cn_cn_next'],
                                                         dg, cg, mlsd_channels, vlsd, minimal, channel_name_list, multi_sources)

        if mlsd_channels:
            if 'MLSD' not in self:
                self['MLSD'] = {}
            if dg not in self['MLSD']:
                self['MLSD'][dg] = {}
            self['MLSD'][dg][cg] = {}
        for MLSDcn in mlsd_channels:
            for cn in self['CN'][dg][cg]:
                if self['CN'][dg][cg][cn]['pointer'] == self['CN'][dg][cg][MLSDcn]['cn_data']:
                    self['MLSD'][dg][cg][MLSDcn] = cn
                    break
        return vlsd

    def read_cn_block(self, fid, pointer, dg, cg, mlsd_channels, vlsd, minimal, channel_name_list, multi_sources=False):
        """reads single Channel block

        Parameters
        ----------------
        fid : float
            file identifier
        pointer : int
            position in file
        dg : int
            data group number
        cg : int
            channel group number in data group
        mlsd_channels : list of int
            list of maximum length data channel numbers
        minimal: flag
            to activate minimum content reading for raw data fetching
        channel_name_list : bool
            Flag to reads only channel blocks for listChannels4 method
        multi_sources: flag
            to read same signal with same name from different source message

        Returns
        -----------
        cn : integer
            channel number
        MLSDChannels list of appended Maximum Length Sampling Data channels
        vlsd : boolean
        """
        temp = CNBlock()
        temp.read_cn(fid=fid, pointer=pointer)
        cn = temp['cn_byte_offset'] * 8 + temp['cn_bit_offset']
        self['CN'][dg][cg][cn] = {}
        self['CN'][dg][cg][cn].update(temp)
        # check for MLSD
        if self['CN'][dg][cg][cn]['cn_type'] == 5:
            mlsd_channels.append(cn)
        # reads Channel Conversion Block
        self['CC'][dg][cg][cn] = CCBlock()
        self['CC'][dg][cg][cn].read_cc(fid, self['CN'][dg][cg][cn]['cn_cc_conversion'])
        if not channel_name_list:
            if not minimal:
                # reads Channel Source Information
                temp = SIBlock()
                temp.read_si(fid, self['CN'][dg][cg][cn]['cn_si_source'])
                if temp:
                    self['CN'][dg][cg][cn]['SI'] = dict()
                    self['CN'][dg][cg][cn]['SI'].update(temp)

            # keep original non unique channel name
            self['CN'][dg][cg][cn]['orig_name'] = self['CN'][dg][cg][cn]['name']
            # check if already existing channel name
            self['CN'][dg][cg][cn]['name'] = \
                self._unique_channel_name(fid, self['CN'][dg][cg][cn]['name'], dg, cg, cn, multi_sources= multi_sources)
            if self.filterChannelNames:
                # filters channels modules
                self['CN'][dg][cg][cn]['name'] = self['CN'][dg][cg][cn]['name'].split('.')[-1]

            # reads Channel Array Block
            if self['CN'][dg][cg][cn]['cn_composition']:
                # composition but can be either structure of channels or array
                fid.seek(self['CN'][dg][cg][cn]['cn_composition'])
                id = fid.read(4)
                if id in ('##CA', b'##CA'):
                    self['CN'][dg][cg][cn]['CABlock'] = CABlock()
                    self['CN'][dg][cg][cn]['CABlock'].read(fid, self['CN'][dg][cg][cn]['cn_composition'])
                elif id in ('##CN', b'##CN'):
                    # channel composition
                    # stores minimal info from upper node channel
                    cn_upper_node = cn
                    cn_next_upper_node = self['CN'][dg][cg][cn]['cn_cn_next']
                    cn, mlsd_channels, vlsd = self.read_cn_block(fid, self['CN'][dg][cg][cn]['cn_composition'],
                                                                 dg, cg, mlsd_channels, vlsd,
                                                                 minimal, channel_name_list)
                    while self['CN'][dg][cg][cn]['cn_cn_next']:
                        cn, mlsd_channels, vlsd = self.read_cn_block(fid, self['CN'][dg][cg][cn]['cn_cn_next'],
                                                                     dg, cg, mlsd_channels, vlsd, minimal,
                                                                     channel_name_list)
                    # restores minimal info from upper node channel
                    cn = cn_upper_node
                    self['CN'][dg][cg][cn]['cn_cn_next'] = cn_next_upper_node
                else:
                    warn('unknown channel composition')

            if self['CN'][dg][cg][cn]['cn_data']:  # channel type specific signal data
                fid.seek(self['CN'][dg][cg][cn]['cn_data'])
                id = fid.read(4)
                if self['CN'][dg][cg][cn]['cn_type'] == 1:
                    if id in ('##SD', b'##SD', '##DZ', b'##DZ', '##DL', b'##DL', '##HL', b'##HL'):
                        # simple VLSD
                        try:
                            self['VLSD'][dg][cg].append(cn)
                        except KeyError:
                            if 'VLSD' not in self:
                                self['VLSD'] = {}
                            if dg not in self['VLSD']:
                                self['VLSD'][dg] = {}
                            if cg not in self['VLSD'][dg]:
                                self['VLSD'][dg][cg] = []
                            self['VLSD'][dg][cg].append(cn)
                    elif id in ('##CG', b'##CG'):
                        # VLSD CG Channel
                        # need to wait all CG blocks are parsed with this Data Group to make the link
                        pass
                elif self['CN'][dg][cg][cn]['cn_type'] == 5:
                    # MLSD Channel
                    mlsd_channels.append(cn)
                if not minimal and id in ('##EV', b'##EV'):
                    # Event signal structure
                    self['CN'][dg][cg][cn]['cn_data'] = self.read_ev_block(fid, self['CN'][dg][cg][cn]['cn_data'])
        return cn, mlsd_channels, vlsd

    def clean_dg_info(self, dg):
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

    @staticmethod
    def read_sr_block(fid, pointer):
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
            sr_blocks = dict()
            sr_blocks[0] = SRBlock(fid, pointer)
            while sr_blocks[sr]['sr_sr_next'] > 0:
                sr += 1
                sr_blocks[sr] = SRBlock(fid, sr_blocks[sr - 1]['sr_sr_next'])
            return sr_blocks

    @staticmethod
    def read_ev_block(fid, pointer):

        """ reads Events Blocks

        Parameters
        ----------------
        fid : identifier
            file identifier
        pointer : int
            position of EVBlock in file

        Returns
        -----------
        Event Blocks in a dict
        """

        ev_blocks = OrderedDict()
        while pointer:
            ev_blocks[pointer] = EVBlock(fid, pointer)
            pointer = ev_blocks[pointer]['ev_ev_next']
        return ev_blocks

    def read_ch_block(self, fid, pointer):
        """reads channel hierarchy Blocks

        Parameters
        ----------------
        fid : identifier
            file identifier
        pointer : int
            position of EVBlock in file

        Returns
        -----------
        channel hierarchy Blocks in a dict
        """
        ch_blocks = list()
        while pointer:
            ch_blocks.append(CHBlock(fid, pointer))
            if ch_blocks[-1]['ch_ch_first']:  # child CHBlock
                ch_blocks[-1]['child'] = self.read_ch_block(fid, ch_blocks[-1]['ch_ch_first'])
            pointer = ch_blocks[-1]['ch_ch_next']
        return ch_blocks

    def list_channels4(self, file_name=None, fid=None):
        """ Read MDF file and extract its complete structure

        Parameters
        ----------------
        file_name : str
            file name
        fid

        Returns
        -----------
        list of channel names contained in file
        """
        if file_name is not None:
            self.fileName = file_name
        # Open file
        if fid is None and file_name is not None:
            # Open file
            (fid, file_name, zipfile) = _open_mdf(self.fileName)
        channel_name_list = []
        # reads Header HDBlock
        self['HD'].update(HDBlock(fid))

        # reads Data Group, channel groups and channel Blocks
        # recursively but not the other metadata block
        self.read_dg_block(fid, True)

        for dg in self['DG']:
            for cg in self['CG'][dg]:
                for cn in self['CN'][dg][cg]:
                    channel_name_list.append(self['CN'][dg][cg][cn]['name'])

        # CLose the file
        fid.close()
        return channel_name_list

    def list_channels_sources4(self, file_name=None, fid=None):
        """ Read MDF file and extract its complete structure

        Parameters
        ----------------
        file_name : str
            file name
        fid

        Returns
        -----------
        list of channel names and sources messages contained in file
        """
        if file_name is not None:
            self.fileName = file_name
        # Open file
        if fid is None and file_name is not None:
            # Open file
            (fid, file_name, zipfile) = _open_mdf(self.fileName)
        channel_name_sources_list = []
        # reads Header HDBlock
        self['HD'].update(HDBlock(fid))

        # reads Data Group, channel groups and channel Blocks
        # recursively but not the other metadata block
        self.read_dg_block(fid, True)

        for dg in self['DG']:
            for cg in self['CG'][dg]:
                for cn in self['CN'][dg][cg]:
                    original_name = self['CN'][dg][cg][cn]['orig_name']
                    if len(self['CG'][dg][cg]['acq_source']['source_path']) == 0:
                        original_source = self['CG'][dg][cg]['acq_name']['Comment']
                    else:
                        original_source = self['CG'][dg][cg]['acq_source']['source_path']['Comment'] \
                        + '.' + self['CG'][dg][cg]['acq_name']['Comment']
                    channel_name_sources_list.append([original_name, original_source])

        # CLose the file
        fid.close()
        return channel_name_sources_list

    def _unique_channel_name(self, fid, name, dg, cg, cn, multi_sources=False):
        """ generate unique channel name

            Parameters
            ----------------
            fid
                file identifier

            name: str
                channel name to be checked

            dg : int
                data group number

            cg: int
                channel group number

            cn : int
                channel number number

            multi_sources: flag
            to read same signal with same name from different source message

            Returns
            -----------
            channel name made unique
        """
        # check if already existing channel name
        if name in self['ChannelNamesByDG'][dg]:  # for unsorted data
            if self['CN'][dg][cg][cn]['cn_si_source']:
                temp = SIBlock()
                temp.read_si(fid, self['CN'][dg][cg][cn]['cn_si_source'])
                if temp['si_tx_name'] > 0:
                    source_name = temp['source_name']['Comment']
                else:
                    source_name = cn
            else:
                source_name = cn
            if multi_sources:
                name = u'{0}'.format(name, dg, source_name)
            else:
                name = u'{0}_{1}_{2}'.format(name, dg, source_name)
        elif name in self['allChannelList']:  # for sorted data
            if self['CN'][dg][cg][cn]['cn_si_source']:
                temp = SIBlock()
                temp.read_si(fid, self['CN'][dg][cg][cn]['cn_si_source'])
                if temp['si_tx_name'] > 0:
                    source_name = temp['source_name']['Comment']
                else:
                    source_name = dg
            else:
                source_name = dg
            if multi_sources and name != 't':
                name = u'{0}'.format(name, dg, source_name)
            else:
                name = u'{0}_{1}_{2}'.format(name, dg, source_name)
        self['ChannelNamesByDG'][dg].add(name)
        self['allChannelList'].add(name)

        if self['CG'][dg][cg]['link_count'] > 6:  # cg master link
            try:
                self['masters'][self['CG'][dg][cg]['cg_cg_master']]['channels'].add(name)
            except KeyError:
                self['masters'][self['CG'][dg][cg]['cg_cg_master']] = dict()
                self['masters'][self['CG'][dg][cg]['cg_cg_master']]['channels'] = set()
                self['masters'][self['CG'][dg][cg]['cg_cg_master']]['name'] = 'master_{}'.format(dg)
            self['CN'][dg][cg][cn]['masterCG'] = self['CG'][dg][cg]['cg_cg_master']
        else:
            try:
                 self['masters'][self['CG'][dg][cg]['pointer']]['channels'].add(name)
            except KeyError:
                self['masters'][self['CG'][dg][cg]['pointer']] = dict()
                self['masters'][self['CG'][dg][cg]['pointer']]['channels'] = set()
                self['masters'][self['CG'][dg][cg]['pointer']]['name'] = 'master_{}'.format(dg)
            if self['CN'][dg][cg][cn]['cn_type'] in (2, 3):  # master channel
                self['masters'][self['CG'][dg][cg]['pointer']]['name'] = name
                self['masters'][self['CG'][dg][cg]['pointer']]['id'] = (dg, cg)
            self['CN'][dg][cg][cn]['masterCG'] = self['CG'][dg][cg]['pointer']
        return name

    def unique_id(self, ndg, ncg, ncn):
        """ generate unique id tuples

            Parameters
            ----------------
            ndg : int
                data group number

            ncg: int
                channel group number

            ncn : int
                channel number

            Returns
            -----------
            tuples: (data group number, channel group number, channel number),
                    (channel name, channel source, channel path),
                    (group name, group source, group path)
        """
        cn = self['CN'][ndg][ncg][ncn]['orig_name']
        try:  # SI block not always existing
            cs = self['CN'][ndg][ncg][ncn]['SI']['source_name']['Comment']
        except KeyError:
            cs = None
        try:
            cp = self['CN'][ndg][ncg][ncn]['SI']['source_path']['Comment']
        except KeyError:
            cp = None
        try:
            gn = self['CG'][ndg][ncg]['acq_name']['Comment']
        except KeyError:
            gn = None
        try:
            gs = self['CG'][ndg][ncg]['acq_source']['source_name']['Comment']
        except KeyError:
            gs = None
        try:
            gp = self['CG'][ndg][ncg]['acq_source']['source_path']['Comment']
        except KeyError:
            gp = None
        return (ndg, ncg, ncn), (cn, cs, cp), (gn, gs, gp)


def _generate_dummy_mdf4(info, channel_list):
    """ computes MasterChannelList and dummy mdf dict from an info object

    Parameters
    ----------------
    info : info object
        information structure of file

    channel_list : list of str
        channel list

    Returns
    -----------
    a dict which keys are master channels in files with values a list of related channels of the raster
    """
    MasterChannelList = {}
    all_channel_list = set()
    mdfdict = {}
    for dg in info['DG']:
        master = ''
        mastertype = 0
        channel_names_by_dg = set()
        for cg in info['CG'][dg]:
            channel_name_list = []
            for cn in info['CN'][dg][cg]:
                if info['CN'][dg][cg][cn]['cn_data_type'] == 13:
                    for name in ('ms', 'minute', 'hour', 'day', 'month', 'year'):
                        channel_name_list.append(name)
                        all_channel_list.add(name)
                        channel_names_by_dg.add(name)
                        mdfdict[name] = {}
                        mdfdict[name][dataField] = None
                        mdfdict[name][descriptionField] = name
                        mdfdict[name][unitField] = name
                        mdfdict[name][idField] = (dg, cg, cn)
                elif info['CN'][dg][cg][cn]['cn_data_type'] == 14:
                    for name in ('ms', 'days'):
                        channel_name_list.append(name)
                        all_channel_list.add(name)
                        channel_names_by_dg.add(name)
                        mdfdict[name] = {}
                        mdfdict[name][dataField] = None
                        mdfdict[name][descriptionField] = name
                        mdfdict[name][unitField] = name
                        mdfdict[name][idField] = (dg, cg, cn)
                else:
                    name = info['CN'][dg][cg][cn]['name']
                    if name in channel_names_by_dg:
                        name = u'{0}_{1}_{2}_{3}'.format(name, dg, cg, cn)
                    elif name in all_channel_list:
                        name = u'{0}_{1}'.format(name, dg)
                    if channel_list is None or name in channel_list:
                        channel_name_list.append(name)
                        all_channel_list.add(name)
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
                        mdfdict[name][idField] = info.unique_id(dg, cg, cn)
                    if info['CN'][dg][cg][cn]['cn_sync_type']:
                        # master channel of cg
                        master = name
                        mastertype = info['CN'][dg][cg][cn]['cn_sync_type']
            for chan in channel_name_list:
                mdfdict[chan][masterField] = master
                mdfdict[chan][masterTypeField] = mastertype
        MasterChannelList[master] = channel_name_list
    return MasterChannelList, mdfdict

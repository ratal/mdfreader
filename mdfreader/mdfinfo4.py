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
    Python version currently running, needed for compatibility of both python 2.6+ and 3.2+

mdfinfo4 module
--------------------------
"""
from __future__ import print_function

from struct import calcsize, unpack, pack
from os.path import dirname, abspath
from sys import version_info, stderr, path
_root = dirname(abspath(__file__))
path.append(_root)
from mdf import _open_MDF, dataField, descriptionField, unitField, masterField, masterTypeField

from numpy import sort, zeros
import time
from xml.etree.ElementTree import Element, SubElement, Comment, tostring, XMLParser, XML, register_namespace
from lxml import etree
namespace = '{http://www.asam.net/mdf/v4}'
_parsernsclean = etree.XMLParser(ns_clean=True)
_find_TX = etree.XPath('/TX') # efficient way to find TX in xml
_find_names = etree.XPath('/names')
_find_linker_name = etree.XPath('/linker_name')
_find_linker_address = etree.XPath('/linker_address')
_find_address = etree.XPath('/address')
_find_axis_monotony = etree.XPath('/axis_monotony')
_find_raster = etree.XPath('/raster')
_find_formula = etree.XPath('/formula')
_find_COMPU_METHOD = etree.XPath('/COMPU_METHOD')
_find_path = etree.XPath('/path')
_find_bus = etree.XPath('/bus')
_find_protocol = etree.XPath('/protocol')
_find_tool_id = etree.XPath('/tool_id')
_find_tool_vendor = etree.XPath('/tool_vendor')
_find_tool_version = etree.XPath('/tool_version')
_find_user_name = etree.XPath('/user_name')
_find_common_properties = etree.XPath('common_properties')

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


class MDFBlock(dict):

    """MDFBlock base class for the MDF related block classes

    Methods
    ------------
    loadHeader(fid, pointer)
        reads block's header and put in class dict
    mdfblockread( fid, type,  count )
        converts a byte array of length count to a given data type
    mdfblockreadCHAR( fid,  count )
        reads a character chain of length count encoded in latin.
    mdfblockreadBYTE( fid, count )
        reads an array of UTF-8 encoded bytes
    """

    def __init__(self):
        """ MDFBlock constructor.
        Reserves header keys in class dict (id, reserved,link_count,pointer)
        """
        self['id'] = {}
        self['reserved'] = {}
        self['length'] = {}
        self['link_count'] = {}
        self['pointer'] = 0

    def loadHeader(self, fid, pointer):
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
            self['id'] = self.mdfblockreadCHAR(fid, 4)
            self['reserved'] = self.mdfblockreadBYTE(fid, 4)
            self['length'] = self.mdfblockread(fid, UINT64, 1)
            self['link_count'] = self.mdfblockread(fid, UINT64, 1)
            self['pointer'] = pointer

    def writeHeader(self, fid, Id, block_length, link_count):
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
        self.writeChar(fid, Id)
        self.writeChar(fid, '\0' * 4)  # reserved
        block_length_pointer = fid.tell()
        fid.write(pack(UINT64, block_length))  # block size
        link_count_pointer = fid.tell()
        fid.write(pack(UINT64, link_count))  # number of links
        return current_position

    @staticmethod
    def mdfblockread(fid, type, count):
        """ converts a byte array of length count to a given data type

        Parameters
        ----------------
        type : str
            C format data type
        count : int
            number of elements to sequentially read

        Returns
        -----------
        array of values of 'type' parameter
        """
        value = fid.read(calcsize(type) * count)
        if value:
            if count == 1:
                return unpack(type, value)[0]
            else:
                if '<' in type or '>' in type:
                    endian = type[0]
                    type = type.strip('<>')
                    return list(unpack(endian+count*type, value))
                else:
                    return list(unpack(count*type, value))
        else:
            return None

    @staticmethod
    def mdfblockreadCHAR(fid, count):
        """ reads a character chain of length count encoded in latin. Removes trailing 0

        Parameters
        ----------------
        count : int
            number of characters to read

        Returns
        -----------
        a string of length count
        """
        if PythonVersion < 3:
            return fid.read(count).rstrip(b'\x00')
        else:
            return fid.read(count).decode('iso8859-1', 'replace').rstrip('\x00')

    @staticmethod
    def writeChar(fid, value, size=None):
        """ Writes a char in a block

        Parameters
        ----------------
        fid : float
            file identifier
        value : str
            char value to write
        size : int
            size that should take this char
        """
        if size is None:
            temp = value
        else:
            if len(value) > size:
                temp = value[:size]
            else:
                temp = value + '\0' * (size - len(value))
            temp += '\0'
        if PythonVersion >= 3:
            temp = temp.encode('latin1', 'replace')
        fid.write(pack('<' + str(len(temp)) + 's', temp))

    @staticmethod
    def mdfblockreadBYTE(fid, count):
        """ reads an array of UTF-8 encoded bytes. Removes trailing 0

        Parameters
        ----------------
        count : int
            number of bytes to read

        Returns
        -----------
        bytes array of length count
        """
        # UTF-8 encoded bytes
        if PythonVersion < 3:
            return fid.read(count).decode('UTF-8', 'ignore').rstrip(b'\x00')
        else:
            return fid.read(count).decode('UTF-8', 'ignore').rstrip('\x00')

    @staticmethod
    def writePointer(fid, pointer, value):
        """ Writes a value at pointer position and comes back to orgianl position

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
        fid.seek(currentPosition)


class IDBlock(MDFBlock):

    """ reads or writes ID Block
    """
    def __init__(self, fid=None):
        if fid is not None:
            self.read(fid)

    def read(self, fid):
        """ reads IDBlock
        """
        fid.seek(0)
        self['id_file'] = self.mdfblockreadCHAR(fid, 8)
        self['id_vers'] = self.mdfblockreadCHAR(fid, 8)
        self['id_prog'] = self.mdfblockreadCHAR(fid, 8)
        self['id_reserved1'] = self.mdfblockreadBYTE(fid, 4)
        self['id_ver'] = self.mdfblockread(fid, UINT16, 1)
        self['id_reserved2'] = self.mdfblockreadBYTE(fid, 30)
        self['id_unfi_flags'] = self.mdfblockread(fid, UINT16, 1)
        self['id_custom_unfi_flags'] = self.mdfblockread(fid, UINT16, 1)
        # treatment of unfinalised file
        if self['id_ver'] > 410 and 'UnFin' in self['id_file']:
            print('  /!\ unfinalised file', file=stderr)
            if self['id_unfi_flags'] & 1:
                print('Update of cycle counters for CG/CA blocks required', file=stderr)
            if self['id_unfi_flags'] & (1 << 1):
                print('Update of cycle counters for SR blocks required', file=stderr)
            if self['id_unfi_flags'] & (1 << 2):
                print('Update of length for last DT block required', file=stderr)
            if self['id_unfi_flags'] & (1 << 3):
                print('Update of length for last RD block required', file=stderr)
            if self['id_unfi_flags'] & (1 << 4):
                print('Update of last DL block in each chained list of DL blocks required', file=stderr)
            if self['id_unfi_flags'] & (1 << 5):
                print('Update of cg_data_bytes and cg_inval_bytes in VLSD CG block required', file=stderr)
            if self['id_unfi_flags'] & (1 << 6):
                print('Update of offset values for VLSD channel required in case a VLSD CG block is used', file=stderr)

    def write(self, fid):
        """ Writes IDBlock
        """
        self.writeChar(fid, 'MDF     ')
        self.writeChar(fid, '4.11    ')
        self.writeChar(fid, 'MDFreadr')
        fid.write(pack(UINT32, 0))  # reserved
        fid.write(pack(UINT16, 411))  # version 4.11
        self.writeChar(fid, '\0' * 30)  # reserved
        fid.write(pack(UINT16, 0))  # finalisation flag
        fid.write(pack(UINT16, 0))  # custom finalisation flag


class HDBlock(MDFBlock):

    """ reads Header block and save in class dict
    """

    def __init__(self, fid=None, pointer=64):
        if fid is not None:
            self.read(fid)

    def read(self, fid=None, pointer=64):
        # block header
        self.loadHeader(fid, pointer)
        # header block
        self['hd_dg_first'] = self.mdfblockread(fid, LINK, 1)
        self['hd_fh_first'] = self.mdfblockread(fid, LINK, 1)
        self['hd_ch_first'] = self.mdfblockread(fid, LINK, 1)
        self['hd_at_first'] = self.mdfblockread(fid, LINK, 1)
        self['hd_ev_first'] = self.mdfblockread(fid, LINK, 1)
        self['hd_md_comment'] = self.mdfblockread(fid, LINK, 1)
        self['hd_start_time_ns'] = self.mdfblockread(fid, UINT64, 1)
        self['hd_tz_offset_min'] = self.mdfblockread(fid, INT16, 1)
        self['hd_dst_offset_min'] = self.mdfblockread(fid, INT16, 1)
        self['hd_time_flags'] = self.mdfblockread(fid, UINT8, 1)
        self['hd_time_class'] = self.mdfblockread(fid, UINT8, 1)
        self['hd_flags'] = self.mdfblockread(fid, UINT8, 1)
        self['hd_reserved'] = self.mdfblockreadBYTE(fid, 1)
        self['hd_start_angle_rad'] = self.mdfblockread(fid, REAL, 1)
        self['hd_start_distance'] = self.mdfblockread(fid, REAL, 1)
        if self['hd_md_comment']:  # if comments exist
            self['Comment'] = CommentBlock(fid, self['hd_md_comment'], 'HD')

    def write(self, fid):
        # write block header
        self.writeHeader(fid, '##HD', 104, 6)
        # link section
        pointers = {}
        pointers['HD'] = {}
        pointers['HD']['DG'] = fid.tell()
        fid.write(pack(LINK, 0))  # first Data group pointer
        pointers['HD']['FH'] = fid.tell()
        fid.write(pack(LINK, 0))  # first file history block pointer
        pointers['HD']['CH'] = fid.tell()
        fid.write(pack(LINK, 0))  # pointer to hierarchy Block file
        pointers['HD']['AT'] = fid.tell()
        fid.write(pack(LINK, 0))  # pointer to attachment Block
        pointers['HD']['EV'] = fid.tell()
        fid.write(pack(LINK, 0))  # pointer to event Block
        pointers['HD']['MD'] = fid.tell()
        fid.write(pack(LINK, 0))  # pointer to comment Block
        # data section
        fid.write(pack(UINT64, int(time.time()*1E9)))  # time in ns
        fid.write(pack(INT16, int(time.timezone/60)))  # timezone offset in min
        fid.write(pack(INT16, int(time.daylight*60)))  # time daylight offest in min
        fid.write(pack(UINT8, 2))  # time flags
        fid.write(pack(UINT8, 0))  # time class
        fid.write(pack(UINT8, 0))  # hd flags
        self.writeChar(fid, '\0')  # reserved
        fid.write(pack(REAL, 0))  # start angle in radians
        fid.write(pack(REAL, 0))  # start distance in meters
        # writes Header comments
        #temp = CommentBlock()
        #pointers.update(temp.write(fid))
        return pointers


class FHBlock(MDFBlock):

    """ reads File History block and save in class dict
    """

    def __init__(self, fid=None, pointer=None):
        if fid is not None:
            self.read(fid, pointer)

    def read(self, fid, pointer):
        # block header
        self.loadHeader(fid, pointer)
        # history block
        self['fh_fh_next'] = self.mdfblockread(fid, LINK, 1)
        self['fh_md_comment'] = self.mdfblockread(fid, LINK, 1)
        self['fh_time_ns'] = self.mdfblockread(fid, UINT64, 1)
        self['fh_tz_offset_min'] = self.mdfblockread(fid, INT16, 1)
        self['fh_dst_offset_min'] = self.mdfblockread(fid, INT16, 1)
        self['fh_time_flags'] = self.mdfblockread(fid, UINT8, 1)
        self['fh_reserved'] = self.mdfblockreadBYTE(fid, 3)
        if self['fh_md_comment']:  # comments exist
            self['Comment'] = CommentBlock(fid, self['fh_md_comment'], 'FH')

    def write(self, fid):
        # write block header
        self.writeHeader(fid, '##FH', 56, 2)
        # link section
        pointers = {}
        pointers['FH'] = {}
        pointers['FH']['FH'] = fid.tell()
        fid.write(pack(LINK, 0))  # No next FH
        pointers['FH']['MD'] = fid.tell()
        fid.write(pack(LINK, 0))  # comment block
        # data section
        fid.write(pack(UINT64, int(time.time()*1E9)))  # time in ns
        fid.write(pack(INT16, int(time.timezone/60)))  # timezone offset in min
        fid.write(pack(INT16, int(time.daylight*60)))  # time daylight offest in min
        fid.write(pack(UINT8, 2))  # time flags
        self.writeChar(fid, '\0' * 3)  # reserved
        return pointers


class CHBlock(MDFBlock):

    """ reads Channel Hierarchy block and saves in class dict
    """

    def __init__(self, fid, pointer):
        # block header
        self.loadHeader(fid, pointer)
        # Channel hierarchy block
        self['ch_ch_next'] = self.mdfblockread(fid, LINK, 1)
        self['ch_ch_first'] = self.mdfblockread(fid, LINK, 1)
        self['ch_tx_name'] = self.mdfblockread(fid, LINK, 1)
        self['ch_md_comment'] = self.mdfblockread(fid, LINK, 1)
        self['ch_element'] = self.mdfblockread(fid, LINK, (self['link_count'] - 4) / 3)
        # data section
        self['ch_element_count'] = self.mdfblockread(fid, LINK, 1)
        self['ch_type'] = self.mdfblockread(fid, UINT8, 1)
        self['ch_reserved'] = self.mdfblockreadBYTE(fid, 3)
        if self['ch_md_comment']:  # comments exist
            self['Comment'] = CommentBlock(fid, self['ch_md_comment'])
        if self['ch_tx_name']:  # text block containing name of hierarchy level
            self['ch_name_level'] = CommentBlock(fid, self['ch_tx_name'])


class CommentBlock(MDFBlock):

    """ reads or writes Comment block and saves in class dict
    """

    def __init__(self, fid=None, pointer=None, MDType=None):
        if fid is not None:
            self.read(fid, pointer, MDType)

    def read(self, fid, pointer, MDType=None):
        """ reads Comment block and saves in class dict

        Notes
        --------
        Can read xml (MD metadata) or text (TX) comments from several kind of blocks
        """
        
        if pointer > 0:
            # block header
            self.loadHeader(fid, pointer)
            if self['id'] in ('##MD', b'##MD'):
                # Metadata block
                self['Comment'] = self.mdfblockreadBYTE(fid, self['length'] - 24)  # [:-1] # removes normal 0 at end
                try:
                    self['xml_tree'] = etree.fromstring(self['Comment'].encode('utf-8'), _parsernsclean)
                except:
                    print('xml metadata malformed', file=stderr)
                    self['xml_tree'] = None
                #self['xml'] = elementTreeToDict(self['xml_tree'])
                # specific action per comment block type, extracts specific tags from xml
                if MDType == 'CN':  # channel comment
                    self['description'] = self.extractXmlField(self['xml_tree'], _find_TX)
                    self['names'] = self.extractXmlField(self['xml_tree'], _find_names)
                    self['linker_name'] = self.extractXmlField(self['xml_tree'], _find_linker_name)
                    self['linker_address'] = self.extractXmlField(self['xml_tree'], _find_linker_address)
                    self['address'] = self.extractXmlField(self['xml_tree'], _find_address)
                    self['axis_monotony'] = self.extractXmlField(self['xml_tree'], _find_axis_monotony)
                    self['raster'] = self.extractXmlField(self['xml_tree'], _find_raster)
                    self['formula'] = self.extractXmlField(self['xml_tree'], _find_formula)
                elif MDType == 'unit':  # channel comment
                    self['unit'] = self.extractXmlField(self['xml_tree'], _find_TX)
                elif MDType == 'HD':  # header comment
                    self['TX'] = self.extractXmlField(self['xml_tree'], _find_TX)
                    tmp = self.extractXmlField(self['xml_tree'], _find_common_properties)
                    tmp = self['xml_tree'].find(namespace + 'common_properties')
                    if tmp is None:
                        tmp = self['xml_tree'].find('common_properties')
                    if tmp is not None:
                        for t in tmp:
                            self[t.attrib['name']] = t.text
                elif MDType == 'FH':  # File History comment
                    self['TX'] = self.extractXmlField(self['xml_tree'], _find_TX)
                    self['tool_id'] = self.extractXmlField(self['xml_tree'], _find_tool_id)
                    self['tool_vendor'] = self.extractXmlField(self['xml_tree'], _find_tool_vendor)
                    self['tool_version'] = self.extractXmlField(self['xml_tree'], _find_tool_version)
                    self['user_name'] = self.extractXmlField(self['xml_tree'], _find_user_name)
                elif MDType == 'SI':
                    self['TX'] = self.extractXmlField(self['xml_tree'], _find_TX)
                    self['names'] = self.extractXmlField(self['xml_tree'], _find_names)
                    self['path'] = self.extractXmlField(self['xml_tree'], _find_path)
                    self['bus'] = self.extractXmlField(self['xml_tree'], _find_bus)
                    self['protocol'] = self.extractXmlField(self['xml_tree'], _find_protocol)
                elif MDType == 'CC':
                    self['TX'] = self.extractXmlField(self['xml_tree'], _find_TX)
                    self['names'] = self.extractXmlField(self['xml_tree'], _find_names)
                    self['COMPU_METHOD'] = self.extractXmlField(self['xml_tree'], _find_COMPU_METHOD)
                    self['formula'] = self.extractXmlField(self['xml_tree'], _find_formula)
                else:
                    if MDType is not None:
                        print('No recognized MDType', file=stderr)
                        print(MDType, file=stderr)
            elif self['id'] in ('##TX', b'##TX'):
                if MDType == 'CN':  # channel comment
                    self['name'] = self.mdfblockreadBYTE(fid, self['length'] - 24)
                else:
                    self['Comment'] = self.mdfblockreadBYTE(fid, self['length'] - 24)

    def extractXmlField(self, xml_tree, find):
        """ Extract Xml field from a xml tree

        Parameters
        ----------------
        xml_tree : xml tree from xml.etree.ElementTree
        field : str

        Returns
        -----------
        field value in xml tree
        """
        try:
            ret = find(xml_tree)
            if ret:
                ret = ret[0].text
            else:
                ret = None
            return ret
        except:
            print('problem parsing metadata', file=stderr)
            return None

    def write(self, fid, data, MDType):
        
        if MDType == 'TX':
            data = data.encode('utf-8', 'replace')
            data += b'\0'
            # make sure block is multiple of 8
            data_lentgth = len(data)
            remain = data_lentgth % 8
            if not remain == 0:
                data += b'\0' * (8 - (remain % 8))
            block_start = self.writeHeader(fid, '##TX', 24 + len(data), 0)
            fid.write(data)
        else:
            register_namespace('', 'http://www.asam.net/mdf/v4')
            if MDType == 'HD':
                root = Element('HDcomment')
                root.set('xmlns', 'http://www.asam.net/mdf/v4')
                TX = SubElement(root, 'TX')
                TX.text = data['comment']
                common_properties = SubElement(root, 'common_properties')
                e = SubElement(common_properties, 'e', attrib={'name':'subject'})
                e.text =  data['subject']
                e = SubElement(common_properties, 'e', attrib={'name':'project'})
                e.text =  data['project']
                e = SubElement(common_properties, 'e', attrib={'name':'department'})
                e.text =  data['organisation']
                e = SubElement(common_properties, 'e', attrib={'name':'author'})
                e.text =  data['author']
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
            data = tostring(root)
            data += b'\0'
            # make sure block is multiple of 8
            data_lentgth = len(data)
            remain = data_lentgth % 8
            if not remain == 0:
                data += b'\0' * (8 - (remain % 8))
            self.writeHeader(fid, '##MD', 24 + len(data), 0)
            block_start = fid.write(data)
        return block_start


def elementTreeToDict(element):
    """ converts xml tree into dictionnary

    Parameters
    ----------------
    element : xml tree from xml.etree.ElementTree

    Returns
    -----------
    dict of xml tree flattened
    """
    return element.tag, \
            dict(map(elementTreeToDict, element)) or element.text


class DGBlock(MDFBlock):

    """ reads Data Group block and saves in class dict
    """

    def __init__(self, fid=None, pointer=None):
        if fid is not None:
            self.read(fid, pointer)

    def read(self, fid, pointer):
        # block header
        self.loadHeader(fid, pointer)
        # Data Group block : Links
        self['dg_dg_next'] = self.mdfblockread(fid, LINK, 1)
        self['dg_cg_first'] = self.mdfblockread(fid, LINK, 1)
        self['dg_data'] = self.mdfblockread(fid, LINK, 1)
        self['dg_md_comment'] = self.mdfblockread(fid, LINK, 1)
        # data section
        self['dg_rec_id_size'] = self.mdfblockread(fid, UINT8, 1)
        self['dg_reserved'] = self.mdfblockreadBYTE(fid, 7)
        if self['dg_md_comment']:  # comments exist
            self['Comment'] = CommentBlock(fid, self['dg_md_comment'])

    def write(self, fid):
        pointers = {}
        # write block header
        pointers['block_start'] = self.writeHeader(fid, '##DG', 64, 4)
        # link section
        pointers['DG'] = fid.tell()
        fid.write(pack(LINK, 0))  # Next Data group pointer
        pointers['CG'] = fid.tell()
        fid.write(pack(LINK, 0))  # first channel group block pointer
        pointers['data'] = fid.tell()
        fid.write(pack(LINK, 0))  # data block pointer
        pointers['TX'] = fid.tell()
        fid.write(pack(LINK, 0))  # comment block pointer
        # data section
        fid.write(pack(UINT8, 0))  # no recordID
        self.writeChar(fid, '\0' * 7)  # reserved
        return pointers


class CGBlock(MDFBlock):

    """ reads Channel Group block and saves in class dict
    """

    def __init__(self, fid=None, pointer=None):
        if fid is not None:
            self.read(fid, pointer)

    def read(self, fid, pointer):
        # block header
        self.loadHeader(fid, pointer)
        # Channel Group block : Links
        self['cg_cg_next'] = self.mdfblockread(fid, LINK, 1)
        self['cg_cn_first'] = self.mdfblockread(fid, LINK, 1)
        self['cg_tx_acq_name'] = self.mdfblockread(fid, LINK, 1)
        self['cg_si_acq_source'] = self.mdfblockread(fid, LINK, 1)
        self['cg_sr_first'] = self.mdfblockread(fid, LINK, 1)
        self['cg_md_comment'] = self.mdfblockread(fid, LINK, 1)
        # data section
        self['cg_record_id'] = self.mdfblockread(fid, UINT64, 1)
        self['cg_cycle_count'] = self.mdfblockread(fid, UINT64, 1)
        self['cg_flags'] = self.mdfblockread(fid, UINT16, 1)
        self['cg_path_separator'] = self.mdfblockread(fid, UINT16, 1)
        self['cg_reserved'] = self.mdfblockreadBYTE(fid, 4)
        self['cg_data_bytes'] = self.mdfblockread(fid, UINT32, 1)
        self['cg_invalid_bytes'] = self.mdfblockread(fid, UINT32, 1)
        if self['cg_md_comment']:  # comments exist
            self['Comment'] = CommentBlock(fid, self['cg_md_comment'])
        if self['cg_tx_acq_name']:  # comments exist
            self['acq_name'] = CommentBlock(fid, self['cg_tx_acq_name'])

    def write(self, fid, cg_cycle_count, cg_data_bytes):
        pointers = {}
        # write block header
        pointers['block_start'] = self.writeHeader(fid, '##CG', 104, 6)
        # link section
        pointers['CG'] = fid.tell()
        fid.write(pack(LINK, 0))  # Next channel group pointer
        pointers['CN'] = fid.tell()
        fid.write(pack(LINK, 0))  # first channel block pointer
        pointers['ACN'] = fid.tell()
        fid.write(pack(LINK, 0))  # acquisition name pointer
        pointers['ACS'] = fid.tell()
        fid.write(pack(LINK, 0))  # acquisition source pointer
        pointers['SR'] = fid.tell()
        fid.write(pack(LINK, 0))  # sample reduction pointer
        pointers['TX'] = fid.tell()
        fid.write(pack(LINK, 0))  # comment pointer
        # data section
        fid.write(pack(UINT64, 0))  # no recordID
        fid.write(pack(UINT64, cg_cycle_count))  # cycle count
        fid.write(pack(UINT16, 0))  # no flags
        fid.write(pack(UINT16, 0))  # no character specified for path separator
        self.writeChar(fid, '\0' * 4)  # reserved
        pointers['cg_data_bytes'] = fid.tell()
        fid.write(pack(UINT32, cg_data_bytes))  # number of bytes taken by data in record
        fid.write(pack(UINT32, 0))  # no invalid bytes
        return pointers


class CNBlock(MDFBlock):

    """ reads Channel block and saves in class dict
    """

    def __init__(self, fid=None, pointer=None):
        if fid is not None:
            self.read(fid, pointer)

    def read(self, fid, pointer):
        if pointer != 0 and pointer is not None:
            fid.seek(pointer)
            self['pointer'] = pointer
            self['id'] = self.mdfblockreadCHAR(fid, 4)
            self['reserved'] = self.mdfblockreadBYTE(fid, 4)
            self['length'] = self.mdfblockread(fid, UINT64, 1)
            self['link_count'] = self.mdfblockread(fid, UINT64, 1)
            # data section
            fid.seek(pointer + self['length'] - 72)
            self['cn_type'] = self.mdfblockread(fid, UINT8, 1)
            self['cn_sync_type'] = self.mdfblockread(fid, UINT8, 1)
            self['cn_data_type'] = self.mdfblockread(fid, UINT8, 1)
            self['cn_bit_offset'] = self.mdfblockread(fid, UINT8, 1)
            self['cn_byte_offset'] = self.mdfblockread(fid, UINT32, 1)
            self['cn_bit_count'] = self.mdfblockread(fid, UINT32, 1)
            self['cn_flags'] = self.mdfblockread(fid, UINT32, 1)
            self['cn_invalid_bit_pos'] = self.mdfblockread(fid, UINT32, 1)
            self['cn_precision'] = self.mdfblockread(fid, UINT8, 1)
            self['cn_reserved'] = self.mdfblockreadBYTE(fid, 1)
            self['cn_attachment_count'] = self.mdfblockread(fid, UINT16, 1)
            self['cn_val_range_min'] = self.mdfblockread(fid, REAL, 1)
            self['cn_val_range_max'] = self.mdfblockread(fid, REAL, 1)
            self['cn_limit_min'] = self.mdfblockread(fid, REAL, 1)
            self['cn_limit_max'] = self.mdfblockread(fid, REAL, 1)
            self['cn_limit_ext_min'] = self.mdfblockread(fid, REAL, 1)
            self['cn_limit_ext_max'] = self.mdfblockread(fid, REAL, 1)
            # Channel Group block : Links
            fid.seek(pointer + 24)
            self['cn_cn_next'] = self.mdfblockread(fid, LINK, 1)
            self['cn_composition'] = self.mdfblockread(fid, LINK, 1)
            self['cn_tx_name'] = self.mdfblockread(fid, LINK, 1)
            self['cn_si_source'] = self.mdfblockread(fid, LINK, 1)
            self['cn_cc_conversion'] = self.mdfblockread(fid, LINK, 1)
            self['cn_data'] = self.mdfblockread(fid, LINK, 1)
            self['cn_md_unit'] = self.mdfblockread(fid, LINK, 1)
            self['cn_md_comment'] = self.mdfblockread(fid, LINK, 1)
            if self['cn_attachment_count'] > 0:
                self['cn_at_reference'] = self.mdfblockread(fid, LINK, self['cn_attachment_count'])
                self['attachment'] = {}
                if self['cn_attachment_count'] > 1:
                    for at in range(self['cn_attachment_count']):
                        self['attachment'][at] = ATBlock(fid, self['cn_at_reference'][at])
                else:
                    self['attachment'][0] = ATBlock(fid, self['cn_at_reference'])
            if self['link_count'] > (8 + self['cn_attachment_count']):
                self['cn_default_x'] = self.mdfblockread(fid, LINK, 3)
            else:
                self['cn_default_x'] = None
            if self['cn_md_comment']:  # comments exist
                self['Comment'] = CommentBlock(fid, self['cn_md_comment'], MDType='CN')
            if self['cn_md_unit']:  # comments exist
                self['unit'] = CommentBlock(fid, self['cn_md_unit'], 'unit')
            if self['cn_tx_name']:  # comments exist
                self['name'] = CommentBlock(fid, self['cn_tx_name'], MDType='CN')['name']
            self['name'] = self['name'].replace(':','')

    def write(self, fid):
        pointers = {}
        # write block header
        pointers['block_start'] = self.writeHeader(fid, '##CN', 160, 8) # no attachement and default X
        # link section
        pointers['CN'] = fid.tell()
        fid.write(pack(LINK, 0))  # Next channel block pointer
        pointers['CN_Comp'] = fid.tell()
        fid.write(pack(LINK, 0))  # composition of channel pointer
        pointers['TX'] = fid.tell()
        fid.write(pack(LINK, 0))  # TXBlock pointer for channel name
        pointers['SI'] = fid.tell()
        fid.write(pack(LINK, 0))  # source SIBlock pointer
        pointers['CC'] = fid.tell()
        fid.write(pack(LINK, 0))  # Conversion Channel CCBlock pointer
        pointers['data'] = fid.tell()
        fid.write(pack(LINK, 0))  # channel data pointer
        pointers['Unit'] = fid.tell()
        fid.write(pack(LINK, 0))  # channel unit comment block pointer
        pointers['Comment'] = fid.tell()
        fid.write(pack(LINK, 0))  # channel comment block pointer
        # no attachment and default X pointers
        # data section
        fid.write(pack(UINT8, self['cn_type']))  # channel type, 0 normal, 2 master
        fid.write(pack(UINT8, self['cn_sync_type']))  # sync type, 0 None, 1 time, 2 angle, 3 distance, 4 index
        fid.write(pack(UINT8, self['cn_data_type']))  # data type
        fid.write(pack(UINT8, self['cn_bit_offset']))  # bit offset
        fid.write(pack(UINT32, self['cn_byte_offset']))  # byte offset
        fid.write(pack(UINT32, self['cn_bit_count']))  # bit count
        fid.write(pack(UINT32, 0))  # no flags
        fid.write(pack(UINT32, 0))  # no flags
        fid.write(pack(UINT8, 0))  # precision
        self.writeChar(fid, '\0')  # reserved
        fid.write(pack(UINT16, 0))  # attachments count
        fid.write(pack(REAL, 0))  # val range min
        fid.write(pack(REAL, 0))  # val range max
        fid.write(pack(REAL, 0))  # val limit min
        fid.write(pack(REAL, 0))  # val limit max
        fid.write(pack(REAL, 0))  # val limit ext min
        fid.write(pack(REAL, 0))  # val limit ext max
        return pointers


class CCBlock(MDFBlock):

    """ reads Channel Conversion block and saves in class dict
    """

    def __init__(self, fid=None, pointer=None):
        if fid is not None:
            self.read(fid, pointer)

    def read(self, fid, pointer):
        # block header
        if pointer != 0 and pointer is not None:
            fid.seek(pointer)
            self['id'] = self.mdfblockreadCHAR(fid, 4)
            self['reserved'] = self.mdfblockreadBYTE(fid, 4)
            self['length'] = self.mdfblockread(fid, UINT64, 1)
            self['link_count'] = self.mdfblockread(fid, UINT64, 1)
            # Channel Conversion block : Links
            self['cc_tx_name'] = self.mdfblockread(fid, LINK, 1)
            self['cc_md_unit'] = self.mdfblockread(fid, LINK, 1)
            self['cc_md_comment'] = self.mdfblockread(fid, LINK, 1)
            self['cc_cc_inverse'] = self.mdfblockread(fid, LINK, 1)
            self['cc_ref'] = self.mdfblockread(fid, LINK, self['link_count'] - 4)
            # data section
            self['cc_type'] = self.mdfblockread(fid, UINT8, 1)
            self['cc_precision'] = self.mdfblockread(fid, UINT8, 1)
            self['cc_flags'] = self.mdfblockread(fid, UINT16, 1)
            self['cc_ref_count'] = self.mdfblockread(fid, UINT16, 1)
            self['cc_val_count'] = self.mdfblockread(fid, UINT16, 1)
            self['cc_phy_range_min'] = self.mdfblockread(fid, REAL, 1)
            self['cc_phy_range_max'] = self.mdfblockread(fid, REAL, 1)
            if self['cc_val_count']:
                self['cc_val'] = self.mdfblockread(fid, REAL, self['cc_val_count'])
            if self['cc_type'] == 3:  # reads Algebraic formula
                self['cc_ref'] = CommentBlock(fid, self['cc_ref'])
            elif self['cc_type']in (7, 8, 9, 10):  # text list
                for i in range(self['cc_ref_count']):
                    fid.seek(self['cc_ref'][i])
                    ID = self.mdfblockreadCHAR(fid, 4)  # find if TX/MD or another CCBlock
                    if ID in ('##TX', '##MD', b'##TX', b'##MD'):  # for algebraic formulae
                        temp = CommentBlock(fid, self['cc_ref'][i])
                        self['cc_ref'][i] = temp['Comment']
                    elif ID in ('##CC', b'##CC'):  # for table conversion
                        # much more complicated nesting conversions !!!
                        self['cc_ref'][i] = CCBlock(fid, self['cc_ref'][i])
            if self['cc_md_comment']:  # comments exist
                self['Comment'] = CommentBlock(fid, self['cc_md_comment'], MDType='CC')
            if self['cc_md_unit']:  # comments exist
                self['unit'] = CommentBlock(fid, self['cc_md_unit'])
            if self['cc_tx_name']:  # comments exist
                self['name'] = CommentBlock(fid, self['cc_tx_name'])
        else:  # no conversion
            self['cc_type'] = 0


class CABlock(MDFBlock):

    """ reads Channel Array block and saves in class dict
    """

    def __init__(self, fid, pointer):
        # block header
        if pointer != 0 and pointer is not None:
            fid.seek(pointer)
            self['id'] = self.mdfblockreadCHAR(fid, 4)
            self['reserved'] = self.mdfblockreadBYTE(fid, 4)
            self['length'] = self.mdfblockread(fid, UINT64, 1)
            self['link_count'] = self.mdfblockread(fid, UINT64, 1)
            # reads data section
            fid.seek(pointer + 24 + self['link_count'] * 8)
            self['ca_type'] = self.mdfblockread(fid, UINT8, 1)
            self['ca_storage'] = self.mdfblockread(fid, UINT8, 1)
            self['ca_ndim'] = self.mdfblockread(fid, UINT16, 1)
            self['ca_flags'] = int(self.mdfblockread(fid, UINT32, 1))
            self['ca_byte_offset_base'] = self.mdfblockread(fid, INT32, 1)
            self['ca_invalid_bit_pos_base'] = self.mdfblockread(fid, UINT32, 1)
            self['ca_dim_size'] = self.mdfblockread(fid, UINT64, self['ca_ndim'])
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
                self['ca_axis_value'] = self.mdfblockread(fid, REAL, self['SNd'])
            if self['ca_storage'] >= 1:
                self['ca_cycle_count'] = self.mdfblockread(fid, UINT64, self['PNd'])
            # Channel Conversion block : Links
            fid.seek(pointer + 24)
            self['ca_composition'] = self.mdfblockread(fid, LINK, 1)  # point to CN for array of structures or CA for array of array
            if self['ca_storage'] == 2:
                self['ca_data'] = self.mdfblockread(fid, LINK, self['PNd'])
            if 1 << 0 & self['ca_flags']:  # bit 0
                self['ca_dynamic_size'] = self.mdfblockread(fid, LINK, self['ca_ndim'] * 3)
            if 1 << 1 & self['ca_flags']:  # bit 1
                self['ca_input_quantity'] = self.mdfblockread(fid, LINK, self['ca_ndim'] * 3)
            if 1 << 2 & self['ca_flags']:  # bit 2
                self['ca_output_quantity'] = self.mdfblockread(fid, LINK, 3)
            if 1 << 3 & self['ca_flags']:  # bit 3
                self['ca_comparison_quantity'] = self.mdfblockread(fid, LINK, 3)
            if 1 << 4 & self['ca_flags']:  # bit 4
                self['ca_cc_axis_conversion'] = self.mdfblockread(fid, LINK, self['ca_ndim'])
            if 1 << 4 & self['ca_flags'] and not 1 << 5 & self['ca_flags']:  # bit 4 and 5
                self['ca_axis'] = self.mdfblockread(fid, LINK, self['ca_ndim'] * 3)
            # nested arrays
            if self['ca_composition']:
                self['CABlock'] = CABlock(fid, self['ca_composition'])


class ATBlock(MDFBlock):

    """ reads Attachment block and saves in class dict
    """

    def __init__(self, fid, pointer):
        # block header
        if pointer != 0 and pointer is not None:
            fid.seek(pointer)
            self['id'] = self.mdfblockreadCHAR(fid, 4)
            self['reserved'] = self.mdfblockreadBYTE(fid, 4)
            self['length'] = self.mdfblockread(fid, UINT64, 1)
            self['link_count'] = self.mdfblockread(fid, UINT64, 1)
            # links
            self['at_at_next'] = self.mdfblockread(fid, LINK, 1)
            self['at_tx_filename'] = self.mdfblockread(fid, LINK, 1)
            self['at_tx_mimetype'] = self.mdfblockread(fid, LINK, 1)
            self['at_md_comment'] = self.mdfblockread(fid, LINK, 1)
            # data
            self['at_flags'] = self.mdfblockread(fid, UINT16, 1)
            self['at_creator_index'] = self.mdfblockread(fid, UINT16, 1)
            self['at_reserved'] = self.mdfblockreadBYTE(fid, 4)
            self['at_md5_checksum'] = self.mdfblockreadBYTE(fid, 16)
            self['at_original_size'] = self.mdfblockread(fid, UINT64, 1)
            self['at_embedded_size'] = self.mdfblockread(fid, UINT64, 1)
            if self['at_embedded_size'] > 0:
                self['at_embedded_data'] = self.mdfblockreadBYTE(fid, self['at_embedded_size'])
            if self['at_md_comment']:  # comments exist
                self['Comment'] = CommentBlock(fid, self['at_md_comment'])


class EVBlock(MDFBlock):

    """ reads Event block and saves in class dict
    """

    def __init__(self, fid, pointer):
        # block header
        if pointer != 0 and pointer is not None:
            fid.seek(pointer)
            self['id'] = self.mdfblockreadCHAR(fid, 4)
            self['reserved'] = self.mdfblockreadBYTE(fid, 4)
            self['length'] = self.mdfblockread(fid, UINT64, 1)
            self['link_count'] = self.mdfblockread(fid, UINT64, 1)
            # data section
            fid.seek(pointer + self['length'] - 32)
            self['ev_type'] = self.mdfblockread(fid, UINT8, 1)
            self['ev_sync_type'] = self.mdfblockread(fid, UINT8, 1)
            self['ev_range_type'] = self.mdfblockread(fid, UINT8, 1)
            self['ev_cause'] = self.mdfblockread(fid, UINT8, 1)
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
            self['ev_flags'] = self.mdfblockread(fid, UINT8, 1)
            self['at_reserved'] = self.mdfblockreadBYTE(fid, 3)
            self['ev_scope_count'] = self.mdfblockread(fid, UINT32, 1)
            self['ev_attachment_count'] = self.mdfblockread(fid, UINT16, 1)
            self['ev_creator_index'] = self.mdfblockread(fid, UINT16, 1)
            self['ev_sync_base_value'] = self.mdfblockread(fid, UINT64, 1)
            self['ev_sync_factor'] = self.mdfblockread(fid, REAL, 1)
            # link section
            fid.seek(pointer + 24)
            self['ev_ev_next'] = self.mdfblockread(fid, LINK, 1)
            self['ev_ev_parent'] = self.mdfblockread(fid, LINK, 1)
            self['ev_ev_range'] = self.mdfblockread(fid, LINK, 1)
            self['ev_tx_name'] = self.mdfblockread(fid, LINK, 1)
            self['ev_md_comment'] = self.mdfblockread(fid, LINK, 1)
            self['ev_scope'] = self.mdfblockread(fid, LINK, self['ev_scope_count'])
            # post treatment
            if self['ev_attachment_count'] > 0:
                self['ev_at_reference'] = self.mdfblockread(fid, LINK, self['ev_attachment_count'])
                for at in range(self['ev_attachment_count']):
                    self['attachment'][at] = ATBlock(fid, self['ev_at_reference'][at])
            if self['ev_md_comment']:  # comments exist
                self['Comment'] = CommentBlock(fid, self['ev_md_comment'])
            if self['ev_tx_name']:  # comments exist
                self['name'] = CommentBlock(fid, self['ev_tx_name'])


class SRBlock(MDFBlock):

    """ reads Sample Reduction block and saves in class dict
    """

    def __init__(self, fid, pointer):
        # block header
        if pointer != 0 and pointer is not None:
            fid.seek(pointer)
            self['id'] = self.mdfblockreadCHAR(fid, 4)
            self['reserved'] = self.mdfblockreadBYTE(fid, 4)
            self['length'] = self.mdfblockread(fid, UINT64, 1)
            self['link_count'] = self.mdfblockread(fid, UINT64, 1)
            # link section
            self['sr_sr_next'] = self.mdfblockread(fid, LINK, 1)
            self['sr_data'] = self.mdfblockread(fid, LINK, 1)
            # data section
            self['sr_cycle_count'] = self.mdfblockread(fid, UINT64, 1)
            self['sr_interval'] = self.mdfblockread(fid, REAL, 1)
            self['sr_sync_type'] = self.mdfblockread(fid, UINT8, 1)
            self['sr_flags'] = self.mdfblockread(fid, UINT8, 1)
            self['sr_reserved'] = self.mdfblockreadBYTE(fid, 6)


class SIBlock(MDFBlock):

    """ reads Source Information block and saves in class dict
    """

    def __init__(self, fid, pointer):
        # block header
        if pointer != 0 and pointer is not None:
            fid.seek(pointer)
            self['id'] = self.mdfblockreadCHAR(fid, 4)
            self['reserved'] = self.mdfblockreadBYTE(fid, 4)
            self['length'] = self.mdfblockread(fid, UINT64, 1)
            self['link_count'] = self.mdfblockread(fid, UINT64, 1)
            # link section
            self['si_tx_name'] = self.mdfblockread(fid, LINK, 1)
            self['si_tx_path'] = self.mdfblockread(fid, LINK, 1)
            self['si_md_comment'] = self.mdfblockread(fid, LINK, 1)
            # data section
            self['si_type'] = self.mdfblockread(fid, UINT8, 1)
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
            self['si_bus_type'] = self.mdfblockread(fid, UINT8, 1)
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
            self['si_flags'] = self.mdfblockread(fid, UINT8, 1)
            self['si_reserved'] = self.mdfblockreadBYTE(fid, 5)
            # post treatment
            self['source_name'] = CommentBlock(fid, self['si_tx_name'])
            self['source_path'] = CommentBlock(fid, self['si_tx_path'])
            self['comment'] = CommentBlock(fid, self['si_md_comment'], MDType='SI')


class info4(dict):

    """ information block parser fo MDF file version 4.x

    Attributes
    --------------
    fileName : str
        name of file

    Notes
    --------
    mdfinfo(FILENAME) contains a dict of structures, for
    each data group, containing key information about all channels in each
    group. FILENAME is a string that specifies the name of the MDF file.
    Either file name or fid should be given.
    General dictionary structure is the following

    - mdfinfo['HDBlock'] header block
    - mdfinfo['DGBlock'][dataGroup] Data Group block
    - mdfinfo['CGBlock'][dataGroup][channelGroup] Channel Group block
    - mdfinfo['CNBlock'][dataGroup][channelGroup][channel] Channel block including text blocks for comment and identifier
    - mdfinfo['CCBlock'][dataGroup][channelGroup][channel] Channel conversion information
    """

    def __init__(self, fileName=None, fid=None):
        """ info4 class constructor

        Parameters
        ----------------
        fileName : str
            file name
        fid : float
            file identifier

        Notes
        ---------
        Either fileName or fid can be used as argument
        """
        self['IDBlock'] = {}  # Identifier Block
        self['HDBlock'] = {}  # Header Block
        self['FHBlock'] = {}
        self['CHBlock'] = {}
        self['DGBlock'] = {}  # Data Group Block
        self['CGBlock'] = {}  # Channel Group Block
        self['CNBlock'] = {}  # Channel Block
        self['CCBlock'] = {}  # Conversion block
        self['ATBlock'] = {}  # Attachment block
        self.fileName = fileName
        if fid is None and fileName is not None:
            # Open file
            (self.fid, self.fileName, self.zipfile) = _open_MDF(self.fileName)
        if self.fileName is not None and fid is None:
            self.readinfo(self.fid)
            # Close the file
            self.fid.close()
            if self.zipfile: # temporary uncompressed file, to be removed
                remove(fileName)
        elif self.fileName is None and fid is not None:
            # called by mdfreader.mdfinfo
            self.readinfo(fid)

    def readinfo(self, fid):
        """ read all file blocks except data

        Parameters
        ----------------
        fid : float
            file identifier
        """
        # reads IDBlock
        self['IDBlock'].update(IDBlock(fid))

        # reads Header HDBlock
        self['HDBlock'].update(HDBlock(fid))

        #print('reads File History blocks, always exists', file=stderr)
        fh = 0  # index of fh blocks
        self['FHBlock'][fh] = {}
        self['FHBlock'][fh] .update(FHBlock(fid, self['HDBlock']['hd_fh_first']))
        while self['FHBlock'][fh]['fh_fh_next']:
            self['FHBlock'][fh + 1] = {}
            self['FHBlock'][fh + 1] .update(FHBlock(fid, self['FHBlock'][fh]['fh_fh_next']))
            fh += 1

        # print('reads Channel Hierarchy blocks', file=stderr)
        if self['HDBlock']['hd_ch_first']:
            ch = 0
            self['CHBlock'][ch] = {}
            self['CHBlock'][ch] .update(CHBlock(fid, self['HDBlock']['hd_ch_first']))
            while self['CHBlock'][ch]['ch_ch_next']:
                self['CHBlock'][ch] .update(CHBlock(fid, self['CHBlock'][ch]['ch_ch_next']))
                ch += 1

        # reads Attachment block
        self['ATBlock'] = self.readATBlock(fid, self['HDBlock']['hd_at_first'])

        # reads Event Block
        if self['HDBlock']['hd_ev_first']:
            ev = 0
            self['EVBlock'] = {}
            self['EVBlock'][ev] = EVBlock(fid, self['HDBlock']['hd_ev_first'])
            while self['EVBlock'][ev]['ev_ev_next']:
                ev += 1
                self['EVBlock'][ev] = EVBlock(fid, self['EVBlock'][ev - 1]['ev_ev_next'])

        # reads Data Group Blocks and recursively the other related blocks
        self.readDGBlock(fid)

    def readDGBlock(self, fid, channelNameList=False):
        """reads Data Group Blocks

        Parameters
        ----------------
        fid : float
            file identifier
        channelNameList : bool
            Flag to reads only channel blocks for listChannels4 method
        """
        self['ChannelNamesByDG'] = {}
        if self['HDBlock']['hd_dg_first']:
            dg = 0
            self['ChannelNamesByDG'][dg] = set()
            self['DGBlock'][dg] = {}
            self['DGBlock'][dg].update(DGBlock(fid, self['HDBlock']['hd_dg_first']))
            # reads Channel Group blocks
            self.readCGBlock(fid, dg, channelNameList)
            while self['DGBlock'][dg]['dg_dg_next']:
                dg += 1
                self['ChannelNamesByDG'][dg] = set()
                self['DGBlock'][dg] = {}
                self['DGBlock'][dg].update(DGBlock(fid, self['DGBlock'][dg - 1]['dg_dg_next']))
                # reads Channel Group blocks
                self.readCGBlock(fid, dg, channelNameList)

    def readCGBlock(self, fid, dg, channelNameList=False):
        """reads Channel Group blocks

        Parameters
        ----------------
        fid : float
            file identifier
        dg : int
            data group number
        channelNameList : bool
            Flag to reads only channel blocks for listChannels4 method
        """
        if self['DGBlock'][dg]['dg_cg_first']:
            cg = 0
            self['CNBlock'][dg] = {}
            self['CNBlock'][dg][cg] = {}
            self['CCBlock'][dg] = {}
            self['CCBlock'][dg][cg] = {}
            self['CGBlock'][dg] = {}
            self['CGBlock'][dg][cg] = {}
            self['CGBlock'][dg][cg].update(CGBlock(fid, self['DGBlock'][dg]['dg_cg_first']))
            VLSDCGBlock = []

            if not channelNameList:
                # reads Source Information Block
                self['CGBlock'][dg][cg]['SIBlock'] = SIBlock(fid, self['CGBlock'][dg][cg]['cg_si_acq_source'])

                # reads Sample Reduction Block
                self['CGBlock'][dg][cg]['SRBlock'] = self.readSRBlock(fid, self['CGBlock'][dg][cg]['cg_sr_first'])

            if not self['CGBlock'][dg][cg]['cg_flags'] & 0b1:  # if not a VLSD channel group
                # reads Channel Block
                self.readCNBlock(fid, dg, cg, channelNameList)
            else:
                VLSDCGBlock.append(cg)

            while self['CGBlock'][dg][cg]['cg_cg_next']:
                cg += 1
                self['CGBlock'][dg][cg] = {}
                self['CGBlock'][dg][cg].update(CGBlock(fid, self['CGBlock'][dg][cg - 1]['cg_cg_next']))
                self['CNBlock'][dg][cg] = {}
                self['CCBlock'][dg][cg] = {}
                if not channelNameList:
                    # reads Source Information Block
                    self['CGBlock'][dg][cg]['SIBlock'] = SIBlock(fid, self['CGBlock'][dg][cg]['cg_si_acq_source'])

                    # reads Sample Reduction Block
                    self['CGBlock'][dg][cg]['SRBlock'] = self.readSRBlock(fid, self['CGBlock'][dg][cg]['cg_sr_first'])

                if not self['CGBlock'][dg][cg]['cg_flags'] & 0b1:  # if not a VLSD channel group
                    # reads Channel Block
                    self.readCNBlock(fid, dg, cg, channelNameList)
                else:
                    VLSDCGBlock.append(cg)

            if VLSDCGBlock:  # VLSD CG Block exiting
                self['VLSD_CG'] = {}
            # Matching VLSD CGBlock with corresponding channel
            for VLSDcg in VLSDCGBlock:
                VLSDCGBlockAdress = self['CGBlock'][dg][VLSDcg]['pointer']
                for cg in self['CGBlock'][dg]:
                    if cg not in VLSDCGBlock:
                        for cn in self['CNBlock'][dg][cg]:
                            if VLSDCGBlockAdress == self['CNBlock'][dg][cg][cn]['cn_data']:
                                # found matching channel with VLSD CGBlock
                                temp = {}
                                temp['cg_cn'] = (cg, cn)
                                self['VLSD_CG'][self['CGBlock'][dg][VLSDcg]['cg_record_id']] = temp
                                break

            # reorder channel blocks and related blocks(CC, SI, AT, CA) based on byte offset
            # this reorder is meant to improve performance while parsing records using core.records.fromfile
            # as it will not use cn_byte_offset
            # first, calculate new mapping/order
            nChannel = len(self['CNBlock'][dg][cg])
            Map = zeros(shape=nChannel, dtype=[('index', 'u4'), ('bit_offset', 'u4')])
            for cn in range(nChannel):
                Map[cn] = (cn, self['CNBlock'][dg][cg][cn]['cn_byte_offset'] * 8 + self['CNBlock'][dg][cg][cn]['cn_bit_offset'])
            orderedMap = sort(Map, order='bit_offset')

            toChangeIndex = Map == orderedMap
            for cn in range(nChannel):
                if not toChangeIndex[cn]:
                    # offset all indexes of indexes to be moved
                    self['CNBlock'][dg][cg][cn + nChannel] = self['CNBlock'][dg][cg].pop(cn)
                    self['CCBlock'][dg][cg][cn + nChannel] = self['CCBlock'][dg][cg].pop(cn)
            for cn in range(nChannel):
                if not toChangeIndex[cn]:
                    # change to ordered index
                    self['CNBlock'][dg][cg][cn] = self['CNBlock'][dg][cg].pop(orderedMap[cn][0] + nChannel)
                    self['CCBlock'][dg][cg][cn] = self['CCBlock'][dg][cg].pop(orderedMap[cn][0] + nChannel)

    def readCNBlock(self, fid, dg, cg, channelNameList=False):
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
        """
        cn = 0
        self['CNBlock'][dg][cg][cn] = {}
        self['CCBlock'][dg][cg][cn] = {}
        self['CNBlock'][dg][cg][cn] = CNBlock(fid, self['CGBlock'][dg][cg]['cg_cn_first'])
        MLSDChannels = []
        # check for MLSD
        if self['CNBlock'][dg][cg][cn]['cn_type'] == 5:
            MLSDChannels.append(cn)
        # check if already existing channel name
        if self['CNBlock'][dg][cg][cn]['name'] in self['ChannelNamesByDG'][dg]:
            self['CNBlock'][dg][cg][cn]['name'] = self['CNBlock'][dg][cg][cn]['name'] + str(cg) + '_'  + str(cn)
        self['ChannelNamesByDG'][dg].add(self['CNBlock'][dg][cg][cn]['name'])

        if self['CGBlock'][dg][cg]['cg_cn_first']:  # Can be NIL for VLSD
            # reads Channel Conversion Block
            self['CCBlock'][dg][cg][cn] = CCBlock(fid, self['CNBlock'][dg][cg][cn]['cn_cc_conversion'])
            if not channelNameList:
                # reads Channel Source Information
                self['CNBlock'][dg][cg][cn]['SIBlock'] = SIBlock(fid, self['CNBlock'][dg][cg][cn]['cn_si_source'])

                # reads Channel Array Block
                if self['CNBlock'][dg][cg][cn]['cn_composition']:  # composition but can be either structure of channels or array
                    fid.seek(self['CNBlock'][dg][cg][cn]['cn_composition'])
                    id = fid.read(4)
                    if id in ('##CA', b'##CA'):
                        self['CNBlock'][dg][cg][cn]['CABlock'] = CABlock(fid, self['CNBlock'][dg][cg][cn]['cn_composition'])
                    elif id in ('##CN', b'##CN'):
                        self['CNBlock'][dg][cg][cn]['CNBlock'] = CNBlock(fid, self['CNBlock'][dg][cg][cn]['cn_composition'])
                    else:
                        raise('unknown channel composition')

                # reads Attachment Block
                if self['CNBlock'][dg][cg][cn]['cn_attachment_count'] > 1:
                    for at in range(self['CNBlock'][dg][cg][cn]['cn_attachment_count']):
                        self['CNBlock'][dg][cg][cn]['attachment'][at].update(self.readATBlock(fid, self['CNBlock'][dg][cg][cn]['cn_at_reference'][at]))
                elif self['CNBlock'][dg][cg][cn]['cn_attachment_count'] == 1:
                    self['CNBlock'][dg][cg][cn]['attachment'][0].update(self.readATBlock(fid, self['CNBlock'][dg][cg][cn]['cn_at_reference']))

            while self['CNBlock'][dg][cg][cn]['cn_cn_next']:
                cn = cn + 1
                self['CNBlock'][dg][cg][cn] = CNBlock(fid, self['CNBlock'][dg][cg][cn - 1]['cn_cn_next'])
                # check for MLSD
                if self['CNBlock'][dg][cg][cn]['cn_type'] == 5:
                    MLSDChannels.append(cn)
                # reads Channel Conversion Block
                self['CCBlock'][dg][cg][cn] = CCBlock(fid, self['CNBlock'][dg][cg][cn]['cn_cc_conversion'])
                if not channelNameList:
                    # reads Channel Source Information
                    self['CNBlock'][dg][cg][cn]['SIBlock'] = SIBlock(fid, self['CNBlock'][dg][cg][cn]['cn_si_source'])

                    # check if already existing channel name
                    if self['CNBlock'][dg][cg][cn]['name'] in self['ChannelNamesByDG'][dg]:
                        self['CNBlock'][dg][cg][cn]['name'] = self['CNBlock'][dg][cg][cn]['name'] + str(cg) + '_'  + str(cn)
                    self['ChannelNamesByDG'][dg].add(self['CNBlock'][dg][cg][cn]['name'])

                    # reads Channel Array Block
                    if self['CNBlock'][dg][cg][cn]['cn_composition']:  # composition but can be either structure of channels or array
                        fid.seek(self['CNBlock'][dg][cg][cn]['cn_composition'])
                        id = fid.read(4)
                        if id in ('##CA', b'##CA'):
                            self['CNBlock'][dg][cg][cn]['CABlock'] = CABlock(fid, self['CNBlock'][dg][cg][cn]['cn_composition'])
                        elif id in ('##CN', b'##CN'):
                            self['CNBlock'][dg][cg][cn]['CNBlock'] = CNBlock(fid, self['CNBlock'][dg][cg][cn]['cn_composition'])
                        else:
                            raise('unknown channel composition')

                    # reads Attachment Block
                    if self['CNBlock'][dg][cg][cn]['cn_attachment_count'] > 1:
                        for at in range(self['CNBlock'][dg][cg][cn]['cn_attachment_count']):
                            print(self['CNBlock'][dg][cg][cn]['cn_at_reference'][at], file=stderr)
                            self['CNBlock'][dg][cg][cn]['attachment'][at].update(self.readATBlock(fid, self['CNBlock'][dg][cg][cn]['cn_at_reference'][at]))
                    elif self['CNBlock'][dg][cg][cn]['cn_attachment_count'] == 1:
                        self['CNBlock'][dg][cg][cn]['attachment'][0].update(self.readATBlock(fid, self['CNBlock'][dg][cg][cn]['cn_at_reference']))

        MLSDChannels = self.readComposition(fid, dg, cg, MLSDChannels, channelNameList=False)

        if MLSDChannels:
            self['MLSD'] = {}
            self['MLSD'][dg] = {}
            self['MLSD'][dg][cg] = {}
        for MLSDcn in MLSDChannels:
            for cn in self['CNBlock'][dg][cg]:
                if self['CNBlock'][dg][cg][cn]['pointer'] == self['CNBlock'][dg][cg][MLSDcn]['cn_data']:
                    self['MLSD'][dg][cg][MLSDcn] = cn
                    break

    def readComposition(self, fid, dg, cg, MLSDChannels, channelNameList=False):
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
        chan = list(self['CNBlock'][dg][cg].keys())[-1] + 1
        for cn in list(self['CNBlock'][dg][cg].keys()):
            if self['CNBlock'][dg][cg][cn]['cn_composition']:
                fid.seek(self['CNBlock'][dg][cg][cn]['cn_composition'])
                ID = MDFBlock()
                ID = ID.mdfblockreadCHAR(fid, 4)
                if ID in ('##CN', b'##CN'):  # Structures
                    self['CNBlock'][dg][cg][chan] = CNBlock(fid, self['CNBlock'][dg][cg][cn]['cn_composition'])
                    self['CCBlock'][dg][cg][chan] = CCBlock(fid, self['CNBlock'][dg][cg][chan]['cn_cc_conversion'])
                    if self['CNBlock'][dg][cg][chan]['cn_type'] == 5:
                        MLSDChannels.append(chan)
                    while self['CNBlock'][dg][cg][chan]['cn_cn_next']:
                        chan += 1
                        self['CNBlock'][dg][cg][chan] = CNBlock(fid, self['CNBlock'][dg][cg][chan - 1]['cn_cn_next'])
                        self['CCBlock'][dg][cg][chan] = CCBlock(fid, self['CNBlock'][dg][cg][chan]['cn_cc_conversion'])
                        if self['CNBlock'][dg][cg][chan]['cn_type'] == 5:
                            MLSDChannels.append(chan)
                    self['CNBlock'][dg][cg][cn]['cn_type'] = 6  # makes the channel virtual
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
        self['HDBlock'].update(HDBlock(fid))

        # reads Data Group, channel groups and channel Blocks  recursively but not the other metadata block
        self.readDGBlock(fid, True)

        for dg in self['DGBlock']:
            for cg in self['CGBlock'][dg]:
                for cn in self['CNBlock'][dg][cg]:
                    channelNameList.append(self['CNBlock'][dg][cg][cn]['name'])

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
    for dg in info['DGBlock']:
            master = ''
            mastertype = 0
            for cg in info['CGBlock'][dg]:
                channelNameList = []
                for cn in info['CNBlock'][dg][cg]:
                    name = info['CNBlock'][dg][cg][cn]['name']
                    if name in allChannelList or \
                            info['CNBlock'][dg][cg][cn]['cn_type'] in (2, 3):
                        name += '_' + str(dg)
                    if channelList is None or name in channelList:
                        channelNameList.append(name)
                        allChannelList.add(name)
                        # create mdf channel
                        mdfdict[name] = {}
                        mdfdict[name][dataField] = None
                        if 'description' in info['CNBlock'][dg][cg][cn]:
                            mdfdict[name][descriptionField] = info['CNBlock'][dg][cg][cn]['description']
                        else:
                            mdfdict[name][descriptionField] = ''
                        if 'unit' in info['CNBlock'][dg][cg][cn]:
                            mdfdict[name][unitField] = info['CNBlock'][dg][cg][cn]['unit']
                        else:
                            mdfdict[name][unitField] = ''
                        mdfdict[name][masterField] = 0 # default is time
                        mdfdict[name][masterTypeField] = None
                    if info['CNBlock'][dg][cg][cn]['cn_sync_type']: # master channel of cg
                        master = name
                        mastertype = info['CNBlock'][dg][cg][cn]['cn_sync_type']
                for chan in channelNameList:
                    mdfdict[chan][masterField] = master
                    mdfdict[chan][masterTypeField] = mastertype
            MasterChannelList[master] = channelNameList
    return (MasterChannelList, mdfdict)
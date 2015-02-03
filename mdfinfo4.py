# -*- coding: utf-8 -*-
""" Measured Data Format blocks paser for version 4.x

Platform and python version
----------------------------------------
With Unix and Windows for python 2.6+ and 3.2+

Created on Sun Dec 15 12:57:28 2013

:Author: `Aymeric Rateau <http://code.google.com/p/mdfreader/>`__

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
from struct import calcsize, unpack, unpack_from
from sys import version_info
from numpy import sort, zeros
PythonVersion=version_info
PythonVersion=PythonVersion[0]

#datatypes
LINK='<Q' 
REAL='<d'
BOOL='<h'
UINT8='<B'
BYTE='<c'
INT16='<h'
UINT16='<H'
UINT32='<I'
INT32='<i'
UINT64='<Q'
INT64='<q'
CHAR='<c'

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
        self['id']={}
        self['reserved']={}
        self['length']={}
        self['link_count']={}
        self['pointer']=0
    
    def loadHeader(self, fid,  pointer):
        """ reads block's header and put in class dict
        
        Parameters
        ----------------
        fid : float
            file identifier
        pointer : int
            position of block in file
        """
        #All blocks have the same header
        if pointer != 0 and pointer is not None:
            fid.seek(pointer)
            self['id']=self.mdfblockreadCHAR(fid, 4)
            self['reserved']=self.mdfblockreadBYTE(fid, 4)
            self['length']=self.mdfblockread(fid, UINT64, 1)
            self['link_count']=self.mdfblockread(fid, UINT64, 1)
            self['pointer']=pointer
        
    @staticmethod
    def mdfblockread( fid, type,  count ):
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
        valueSize=calcsize(type)
        value=fid.read(valueSize*count)
        if value:
            if count==1:
                value=unpack( type,value)
                value=value[0]
            else:
                vect={}
                for i in range(count):
                    vect[i]=unpack_from( type,value, offset=valueSize*i)
                return vect
        else:
            value=None
        return value
        
    @staticmethod
    def mdfblockreadCHAR( fid,  count ):
        """ reads a character chain of length count encoded in latin. Removes trailing 0
        
        Parameters
        ----------------
        count : int
            number of characters to read
            
        Returns
        -----------
        a string of length count
        """
        value=fid.read( count )
        if PythonVersion<3:
            value.replace('\x00', '')
        else:
            value.decode('iso8859-1', 'replace').replace('\x00', '') 
        return value
        
    @staticmethod
    def mdfblockreadBYTE( fid, count ):
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
        value=fid.read( count )
        value=value.decode('UTF-8', 'ignore')
        if PythonVersion<3:
            value=value.replace(b'\x00', b'') 
        else:
            value=value.replace('\x00', '') 
        return value
        
class IDBlock(MDFBlock):
    """ reads ID Block and save in class dict
    """
    def __init__(self, fid):
        # reads IDBlock
        fid.seek(0)
        self['id_file']=self.mdfblockreadCHAR(fid, 8)
        self['id_vers']=self.mdfblockreadCHAR(fid, 8)
        self['id_prog']=self.mdfblockreadCHAR(fid, 8)
        self['id_reserved1']=self.mdfblockreadBYTE(fid, 4)
        self['id_ver']=self.mdfblockread(fid, UINT16, 1)
        self['id_reserved2']=self.mdfblockreadBYTE(fid, 34)

class HDBlock(MDFBlock):
    """ reads Header block and save in class dict
    """
    def __init__(self, fid,  pointer=64):
        # block header
        self.loadHeader(fid, pointer)
        # header block
        self['hd_dg_first']=self.mdfblockread(fid, LINK, 1)
        self['hd_fh_first']=self.mdfblockread(fid, LINK, 1)
        self['hd_ch_first']=self.mdfblockread(fid, LINK, 1)
        self['hd_at_first']=self.mdfblockread(fid, LINK, 1)
        self['hd_ev_first']=self.mdfblockread(fid, LINK, 1)
        self['hd_md_comment']=self.mdfblockread(fid, LINK, 1)
        self['hd_start_time_ns']=self.mdfblockread(fid, UINT64, 1)
        self['hd_tz_offset_min']=self.mdfblockread(fid, INT16, 1)
        self['hd_dst_offset_min']=self.mdfblockread(fid, INT16, 1)
        self['hd_time_flags']=self.mdfblockread(fid, UINT8, 1)
        self['hd_time_class']=self.mdfblockread(fid, UINT8, 1)
        self['hd_flags']=self.mdfblockread(fid, UINT8, 1)
        self['hd_reserved']=self.mdfblockreadBYTE(fid, 1)
        self['hd_start_angle_rad']=self.mdfblockread(fid, REAL, 1)
        self['hd_start_distance']=self.mdfblockread(fid, REAL, 1)
        if self['hd_md_comment']: # if comments exist
            self['Comment']=CommentBlock(fid, self['hd_md_comment'], 'HD')

class FHBlock(MDFBlock):
    """ reads File History block and save in class dict
    """
    def __init__(self, fid,  pointer):
        # block header
        self.loadHeader(fid, pointer)
        # history block
        self['fh_fh_next']=self.mdfblockread(fid, LINK, 1)
        self['fh_md_comment']=self.mdfblockread(fid, LINK, 1)
        self['fh_time_ns']=self.mdfblockread(fid, UINT64, 1)
        self['fh_tz_offset_min']=self.mdfblockread(fid, INT16, 1)
        self['fh_dst_offset_min']=self.mdfblockread(fid, INT16, 1)
        self['fh_time_flags']=self.mdfblockread(fid, UINT8, 1)
        self['fh_reserved']=self.mdfblockreadBYTE(fid, 3)
        if self['fh_md_comment']: # comments exist
            self['Comment']=CommentBlock(fid, self['fh_md_comment'], 'FH')
            
class CHBlock(MDFBlock):
    """ reads Channel Hierarchy block and saves in class dict
    """
    def __init__(self, fid,  pointer):
        # block header
        self.loadHeader(fid, pointer)
        # Channel hierarchy block
        self['ch_ch_next']=self.mdfblockread(fid, LINK, 1)
        self['ch_ch_first']=self.mdfblockread(fid, LINK, 1)
        self['ch_tx_name']=self.mdfblockread(fid, LINK, 1)
        self['ch_md_comment']=self.mdfblockread(fid, LINK, 1)
        self['ch_element']=self.mdfblockread(fid, LINK, (self['link_count']-4)/3)
        # data section
        self['ch_element_count']=self.mdfblockread(fid, LINK, 1)
        self['ch_type']=self.mdfblockread(fid, UINT8, 1)
        self['ch_reserved']=self.mdfblockreadBYTE(fid, 3)
        if self['ch_md_comment']: # comments exist
            self['Comment']=CommentBlock(fid, self['ch_md_comment'])
        if self['ch_tx_name']: # text block containing name of hierarchy level
            self['ch_name_level']=CommentBlock(fid, self['ch_tx_name'])

class CommentBlock(MDFBlock):
    """ reads Comment block and saves in class dict
    
    Notes
    --------
    Can read xml (MD metadata) or text (TX) comments from several kind of blocks
    """
    def __init__(self, fid,  pointer, MDType=None):
        self.namespace='{http://www.asam.net/mdf/v4}'
        if pointer>0:
            # block header
            self.loadHeader(fid, pointer)
            if self['id'] in ('##MD',b'##MD'):
                # Metadata block
                self['Comment']=self.mdfblockreadBYTE(fid, self['length']-24) #[:-1] # removes normal 0 at end
                import xml.etree.ElementTree as ET
                utf8parser = ET.XMLParser(encoding="utf-8")
                try:
                    self['xml_tree']=ET.XML(self['Comment'].encode('utf-8'), parser=utf8parser)
                    self['xml']=elementTreeToDict(self['xml_tree'])
                    # specific action per comment block type, extracts specific tags from xml
                    if MDType=='CN': # channel comment
                        self['description']=self.extractXmlField(self['xml_tree'], 'TX')
                        self['names']=self.extractXmlField(self['xml_tree'], 'names')
                        self['linker_name']=self.extractXmlField(self['xml_tree'], 'linker_name')
                        self['linker_address']=self.extractXmlField(self['xml_tree'], 'linker_address')
                        self['address']=self.extractXmlField(self['xml_tree'], 'address')
                        self['axis_monotony']=self.extractXmlField(self['xml_tree'], 'axis_monotony')
                        self['raster']=self.extractXmlField(self['xml_tree'], 'raster')
                        self['formula']=self.extractXmlField(self['xml_tree'], 'formula')
                    elif MDType=='unit': # channel comment
                        self['unit']=self.extractXmlField(self['xml_tree'], 'TX')
                    elif MDType=='HD': # header comment
                        self['TX']=self.extractXmlField(self['xml_tree'], 'TX')
                        tmp=self['xml_tree'].find(self.namespace+'common_properties')
                        if tmp is not None:
                            self[tmp[0].attrib['name']]=tmp[0].text # subject
                            self[tmp[1].attrib['name']]=tmp[1].text # project
                            self[tmp[2].attrib['name']]=tmp[2].text # department
                            self[tmp[3].attrib['name']]=tmp[3].text # author
                    elif MDType=='FH': # File History comment
                        self['TX']=self.extractXmlField(self['xml_tree'], 'TX')
                        self['tool_id']=self.extractXmlField(self['xml_tree'], 'tool_id')
                        self['tool_vendor']=self.extractXmlField(self['xml_tree'], 'tool_vendor')
                        self['tool_version']=self.extractXmlField(self['xml_tree'], 'tool_version')
                        self['user_name']=self.extractXmlField(self['xml_tree'], 'user_name')
                    elif MDType=='SI':
                        self['TX']=self.extractXmlField(self['xml_tree'], 'TX')
                        self['names']=self.extractXmlField(self['xml_tree'], 'names')
                        self['path']=self.extractXmlField(self['xml_tree'], 'path')
                        self['bus']=self.extractXmlField(self['xml_tree'], 'bus')
                        self['protocol']=self.extractXmlField(self['xml_tree'], 'protocol')
                    elif MDType=='CC':
                        self['TX']=self.extractXmlField(self['xml_tree'], 'TX')
                        self['names']=self.extractXmlField(self['xml_tree'], 'names')
                        self['ho:COMPU_METHOD']=self.extractXmlField(self['xml_tree'], 'ho:COMPU_METHOD')
                        self['formula']=self.extractXmlField(self['xml_tree'], 'formula')
                    else:
                        if MDType is not None:
                            print('No recognized MDType')
                            print(MDType)
                except:
                    print('problem parsing metadata for '+MDType+' block')
                    print(self['Comment'])

            elif self['id'] in ('##TX',b'##TX'):
                if MDType=='CN': # channel comment
                    self['name']=self.mdfblockreadBYTE(fid, self['length']-24)
                else:
                    self['Comment']=self.mdfblockreadBYTE(fid, self['length']-24)
    
    def extractXmlField(self,  xml_tree, field):
        """ Extract Xml field from a xml tree
        
        Parameters
        ----------------
        xml_tree : xml tree from xml.etree.ElementTree
        field : str
        
        Returns
        -----------
        field value in xml tree
        """
        return xml_tree.findtext(self.namespace+field)

def elementTreeToDict(element):
    """ converts xml tree into dictionnary
    
    Parameters
    ----------------
    element : xml tree from xml.etree.ElementTree
    
    Returns
    -----------
    dict of xml tree flattened
    """
    node = dict()

    text = getattr(element, 'text', None)
    if text is not None:
        node['text'] = text

    node.update(element.items()) # element's attributes

    child_nodes = {}
    for child in element: # element's children
        child_nodes.setdefault(child, []).append( elementTreeToDict(child) )

    # convert all single-element lists into non-lists
    for key, value in child_nodes.items():
        if len(value) == 1:
             child_nodes[key] = value[0]

    node.update(child_nodes.items())

    return node
    
class DGBlock(MDFBlock):
    """ reads Data Group block and saves in class dict
    """
    def __init__(self, fid,  pointer):
        # block header
        self.loadHeader(fid, pointer)
        # Data Group block : Links
        self['dg_dg_next']=self.mdfblockread(fid, LINK, 1)
        self['dg_cg_first']=self.mdfblockread(fid, LINK, 1)
        self['dg_data']=self.mdfblockread(fid, LINK, 1)
        self['dg_md_comment']=self.mdfblockread(fid, LINK, 1)
        # data section
        self['dg_rec_id_size']=self.mdfblockread(fid, UINT8, 1)
        self['dg_reserved']=self.mdfblockreadBYTE(fid, 7)
        if self['dg_md_comment']: # comments exist
            self['Comment']=CommentBlock(fid, self['dg_md_comment'])
        
class CGBlock(MDFBlock):
    """ reads Channel Group block and saves in class dict
    """
    def __init__(self, fid,  pointer):
        # block header
        self.loadHeader(fid, pointer)
        # Channel Group block : Links
        self['cg_cg_next']=self.mdfblockread(fid, LINK, 1)
        self['cg_cn_first']=self.mdfblockread(fid, LINK, 1)
        self['cg_tx_acq_name']=self.mdfblockread(fid, LINK, 1)
        self['cg_si_acq_source']=self.mdfblockread(fid, LINK, 1)
        self['cg_sr_first']=self.mdfblockread(fid, LINK, 1)
        self['cg_md_comment']=self.mdfblockread(fid, LINK, 1)
        # data section
        self['cg_record_id']=self.mdfblockread(fid, UINT64, 1)
        self['cg_cycle_count']=self.mdfblockread(fid, UINT64, 1)
        self['cg_flags']=self.mdfblockread(fid, UINT16, 1)
        self['cg_path_separator']=self.mdfblockread(fid, UINT16, 1)
        self['cg_reserved']=self.mdfblockreadBYTE(fid, 4)
        self['cg_data_bytes']=self.mdfblockread(fid, UINT32, 1)
        self['cg_invalid_bytes']=self.mdfblockread(fid, UINT32, 1)
        if self['cg_md_comment']: # comments exist
            self['Comment']=CommentBlock(fid, self['cg_md_comment'])
        if self['cg_tx_acq_name']: # comments exist
            self['acq_name']=CommentBlock(fid, self['cg_tx_acq_name'])
        
class CNBlock(MDFBlock):
    """ reads Channel block and saves in class dict
    """
    def __init__(self, fid,  pointer):
        if pointer != 0 and pointer is not None:
            fid.seek( pointer )
            self['pointer']=pointer
            self['id']=self.mdfblockreadCHAR(fid, 4)
            self['reserved']=self.mdfblockreadBYTE(fid, 4)
            self['length']=self.mdfblockread(fid, UINT64, 1)
            self['link_count']=self.mdfblockread(fid, UINT64, 1)
            # data section
            fid.seek(pointer+self['length']-72)
            self['cn_type']=self.mdfblockread(fid, UINT8, 1)
            self['cn_sync_type']=self.mdfblockread(fid, UINT8, 1)
            self['cn_data_type']=self.mdfblockread(fid, UINT8, 1)
            self['cn_bit_offset']=self.mdfblockread(fid, UINT8, 1)
            self['cn_byte_offset']=self.mdfblockread(fid, UINT32, 1)
            self['cn_bit_count']=self.mdfblockread(fid, UINT32, 1)
            self['cn_flags']=self.mdfblockread(fid, UINT32, 1)
            self['cn_invalid_bit_pos']=self.mdfblockread(fid, UINT32, 1)
            self['cn_precision']=self.mdfblockread(fid, UINT8, 1)
            self['cn_reserved']=self.mdfblockreadBYTE(fid, 1)
            self['cn_attachment_count']=self.mdfblockread(fid, UINT16, 1)
            self['cn_val_range_min']=self.mdfblockread(fid, REAL, 1)
            self['cn_val_range_max']=self.mdfblockread(fid, REAL, 1)
            self['cn_limit_min']=self.mdfblockread(fid, REAL, 1)
            self['cn_limit_max']=self.mdfblockread(fid, REAL, 1)
            self['cn_limit_ext_min']=self.mdfblockread(fid, REAL, 1)
            self['cn_limit_ext_max']=self.mdfblockread(fid, REAL, 1)
            # Channel Group block : Links
            fid.seek(pointer+24)
            self['cn_cn_next']=self.mdfblockread(fid, LINK, 1)
            self['cn_composition']=self.mdfblockread(fid, LINK, 1)
            self['cn_tx_name']=self.mdfblockread(fid, LINK, 1)
            self['cn_si_source']=self.mdfblockread(fid, LINK, 1)
            self['cn_cc_conversion']=self.mdfblockread(fid, LINK, 1)
            self['cn_data']=self.mdfblockread(fid, LINK, 1)
            self['cn_md_unit']=self.mdfblockread(fid, LINK, 1)
            self['cn_md_comment']=self.mdfblockread(fid, LINK, 1)
            if self['cn_attachment_count']>0:
                self['cn_at_reference']=self.mdfblockread(fid, LINK, self['cn_attachment_count'])
                self['attachment']={}
                if self['cn_attachment_count']>1:
                    for at in range(self['cn_attachment_count']):
                        self['attachment'][at]=ATBlock(fid, self['cn_at_reference'][at][0])
                else:
                    self['attachment'][0]=ATBlock(fid, self['cn_at_reference'])
            if self['link_count']>(8+self['cn_attachment_count']):
                self['cn_default_x']=self.mdfblockread(fid, LINK, 3)
            else:
                self['cn_default_x']=None
            if self['cn_md_comment']: # comments exist
                self['Comment']=CommentBlock(fid, self['cn_md_comment'], MDType='CN')
            if self['cn_md_unit']: # comments exist
                self['unit']=CommentBlock(fid, self['cn_md_unit'], 'unit')
            if self['cn_tx_name']: # comments exist
                self['name']=CommentBlock(fid, self['cn_tx_name'], MDType='CN')['name']

class CCBlock(MDFBlock):
    """ reads Channel Conversion block and saves in class dict
    """
    def __init__(self, fid,  pointer):
        # block header
        if pointer != 0 and not pointer == None:
            fid.seek( pointer )
            self['id']=self.mdfblockreadCHAR(fid, 4)
            self['reserved']=self.mdfblockreadBYTE(fid, 4)
            self['length']=self.mdfblockread(fid, UINT64, 1)
            self['link_count']=self.mdfblockread(fid, UINT64, 1)
            # Channel Conversion block : Links
            self['cc_tx_name']=self.mdfblockread(fid, LINK, 1)
            self['cc_md_unit']=self.mdfblockread(fid, LINK, 1)
            self['cc_md_comment']=self.mdfblockread(fid, LINK, 1)
            self['cc_cc_inverse']=self.mdfblockread(fid, LINK, 1)
            self['cc_ref']=self.mdfblockread(fid, LINK, self['link_count']-4)
            # data section
            self['cc_type']=self.mdfblockread(fid, UINT8, 1)
            self['cc_precision']=self.mdfblockread(fid, UINT8, 1)
            self['cc_flags']=self.mdfblockread(fid, UINT16, 1)
            self['cc_ref_count']=self.mdfblockread(fid, UINT16, 1)
            self['cc_val_count']=self.mdfblockread(fid, UINT16, 1)
            self['cc_phy_range_min']=self.mdfblockread(fid, REAL, 1)
            self['cc_phy_range_max']=self.mdfblockread(fid, REAL, 1)
            if self['cc_val_count']:
                self['cc_val']=self.mdfblockread(fid, REAL, self['cc_val_count'])
            if self['cc_type']==3: # reads Algebraic formula
                self['cc_ref']=CommentBlock(fid, self['cc_ref'])
            elif self['cc_type']in (7, 8, 9, 10): # text list
                for i in range(self['cc_ref_count']):
                    fid.seek( self['cc_ref'][i][0] )
                    ID=self.mdfblockreadCHAR(fid, 4) # find if TX/MD or another CCBlock
                    if ID in ('##TX', '##MD', b'##TX', b'##MD'): # for algebraic formulae
                        temp=CommentBlock(fid, self['cc_ref'][i][0])
                        self['cc_ref'][i]=temp['Comment']
                    elif ID in ('##CC', b'##CC'): # for table conversion
                        # much more complicated nesting conversions !!!
                        self['cc_ref'][i]=CCBlock(fid, self['cc_ref'][i][0])['name']['Comment']
            if self['cc_md_comment']: # comments exist
                self['Comment']=CommentBlock(fid, self['cc_md_comment'], MDType='CC')
            if self['cc_md_unit']: # comments exist
                self['unit']=CommentBlock(fid, self['cc_md_unit'])
            if self['cc_tx_name']: # comments exist
                self['name']=CommentBlock(fid, self['cc_tx_name'])
        else: # no conversion
            self['cc_type']=0

class CABlock(MDFBlock):
    """ reads Channel Array block and saves in class dict
    """
    def __init__(self, fid,  pointer):
        # block header
        if pointer != 0 and not pointer == None:
            fid.seek( pointer )
            self['id']=self.mdfblockreadCHAR(fid, 4)
            self['reserved']=self.mdfblockreadBYTE(fid, 4)
            self['length']=self.mdfblockread(fid, UINT64, 1)
            self['link_count']=self.mdfblockread(fid, UINT64, 1)
            # reads data section
            fid.seek( pointer +24 + self['link_count']*8)
            self['ca_type']=self.mdfblockread(fid, UINT8, 1) 
            self['ca_storage']=self.mdfblockread(fid, UINT8, 1) 
            self['ca_ndim']=self.mdfblockread(fid, UINT16, 1) 
            self['ca_flags']=int(self.mdfblockread(fid, UINT32, 1) )
            self['ca_byte_offset_base']=self.mdfblockread(fid, INT32, 1) 
            self['ca_invalid_bit_pos_base']=self.mdfblockread(fid, UINT32, 1)
            self['ca_dim_size']=self.mdfblockread(fid, UINT64, self['ca_ndim'])
            try: # more than one dimension, processing dict
                SNd=0
                PNd=1
                for x in self['ca_dim_size']:
                    SNd+=self['ca_dim_size'][x][0]
                    PNd*=self['ca_dim_size'][x][0]
            except: # only one dimension, processing int
                SNd=self['ca_dim_size']
                PNd=SNd
            if 1<<5 & self['ca_flags']: # bit5
                self['ca_axis_value']=self.mdfblockread(fid, REAL, SNd)
            if self['ca_storage']>=1:
                self['ca_cycle_count']=self.mdfblockread(fid, UINT64, PNd)
            # Channel Conversion block : Links
            fid.seek( pointer +24 )
            self['ca_composition']=self.mdfblockread(fid, LINK, 1) # point to CN for array of structures or CA for array of array
            if self['ca_storage']==2:
                self['ca_data']=self.mdfblockread(fid, LINK, PNd) 
            if 1<<0 &self['ca_flags']: # bit 0
                self['ca_dynamic_size']=self.mdfblockread(fid, LINK, self['ca_ndim']*3)
            if 1<<1 &self['ca_flags']: # bit 1
                self['ca_input_quantity']=self.mdfblockread(fid, LINK, self['ca_ndim']*3)
            if 1<<2 &self['ca_flags']: # bit 2
                self['ca_output_quantity']=self.mdfblockread(fid, LINK, 3)
            if 1<<3 &self['ca_flags']: # bit 3
                self['ca_comparison_quantity']=self.mdfblockread(fid, LINK, 3)
            if 1<<4 &self['ca_flags']: # bit 4
                self['ca_cc_axis_conversion']=self.mdfblockread(fid, LINK, self['ca_ndim'])
            if 1<<4 &self['ca_flags'] and not 1<<5 &self['ca_flags']: # bit 4 and 5
                self['ca_axis']=self.mdfblockread(fid, LINK, self['ca_ndim']*3)
            # nested arrays
            if self['ca_composition']:
                self['CABlock']=CABlock(fid, self['ca_composition'])

class ATBlock(MDFBlock):
    """ reads Attachment block and saves in class dict
    """
    def __init__(self, fid,  pointer):
        # block header
        if pointer != 0 and not pointer == None:
            fid.seek( pointer )
            self['id']=self.mdfblockreadCHAR(fid, 4)
            self['reserved']=self.mdfblockreadBYTE(fid, 4)
            self['length']=self.mdfblockread(fid, UINT64, 1)
            self['link_count']=self.mdfblockread(fid, UINT64, 1)
            # links
            self['at_at_next']=self.mdfblockread(fid, LINK, 1)
            self['at_tx_filename']=self.mdfblockread(fid, LINK, 1)
            self['at_tx_mimetype']=self.mdfblockread(fid, LINK, 1)
            self['at_md_comment']=self.mdfblockread(fid, LINK, 1)
            # data
            self['at_flags']=self.mdfblockread(fid, UINT16, 1)
            self['at_creator_index']=self.mdfblockread(fid, UINT16, 1)
            self['at_reserved']=self.mdfblockreadBYTE(fid, 4)
            self['at_md5_checksum']=self.mdfblockreadBYTE(fid, 16)
            self['at_original_size']=self.mdfblockread(fid, UINT64, 1)
            self['at_embedded_size']=self.mdfblockread(fid, UINT64, 1)
            if self['at_embedded_size']>0:
                self['at_embedded_data']=self.mdfblockreadBYTE(fid, self['at_embedded_size'])
            if self['at_md_comment']: # comments exist
                self['Comment']=CommentBlock(fid, self['at_md_comment'])
        
class EVBlock(MDFBlock):
    """ reads Event block and saves in class dict
    """
    def __init__(self, fid,  pointer):
        # block header
        if pointer != 0 and not pointer == None:
            fid.seek( pointer )
            self['id']=self.mdfblockreadCHAR(fid, 4)
            self['reserved']=self.mdfblockreadBYTE(fid, 4)
            self['length']=self.mdfblockread(fid, UINT64, 1)
            self['link_count']=self.mdfblockread(fid, UINT64, 1)
            #data section
            fid.seek( pointer +self['length']-32)
            self['ev_type']=self.mdfblockread(fid, UINT8, 1)
            self['ev_sync_type']=self.mdfblockread(fid, UINT8, 1)
            self['ev_range_type']=self.mdfblockread(fid, UINT8, 1)
            self['ev_cause']=self.mdfblockread(fid, UINT8, 1)
            if self['ev_cause']==0:
                self['ev_cause']='OTHER'
            elif self['ev_cause']==1:
                self['ev_cause']=='ERROR'
            elif self['ev_cause']==2:
                self['ev_cause']=='TOOL'
            elif self['ev_cause']==3:
                self['ev_cause']=='SCRIPT'
            elif self['ev_cause']==4:
                self['ev_cause']=='USER'
            self['ev_flags']=self.mdfblockread(fid, UINT8, 1)
            self['at_reserved']=self.mdfblockreadBYTE(fid, 3)
            self['ev_scope_count']=self.mdfblockread(fid, UINT32, 1)
            self['ev_attachment_count']=self.mdfblockread(fid, UINT16, 1)
            self['ev_creator_index']=self.mdfblockread(fid, UINT16, 1)
            self['ev_sync_base_value']=self.mdfblockread(fid, UINT64, 1)
            self['ev_sync_factor']=self.mdfblockread(fid, REAL, 1)
            # link section
            fid.seek( pointer +24)
            self['ev_ev_next']=self.mdfblockread(fid, LINK, 1)
            self['ev_ev_parent']=self.mdfblockread(fid, LINK, 1)
            self['ev_ev_range']=self.mdfblockread(fid, LINK, 1)
            self['ev_tx_name']=self.mdfblockread(fid, LINK, 1)
            self['ev_md_comment']=self.mdfblockread(fid, LINK, 1)
            self['ev_scope']=self.mdfblockread(fid, LINK, self['ev_scope_count'])
            # post treatment
            if self['ev_attachment_count']>0:
                self['ev_at_reference']=self.mdfblockread(fid, LINK, self['ev_attachment_count'])
                for at in range(self['ev_attachment_count']):
                    self['attachment'][at]=ATBlock(fid, self['ev_at_reference'][at])
            if self['ev_md_comment']: # comments exist
                self['Comment']=CommentBlock(fid, self['ev_md_comment'])
            if self['ev_tx_name']: # comments exist
                self['name']=CommentBlock(fid, self['ev_tx_name'])

class SRBlock(MDFBlock):
    """ reads Sample Reduction block and saves in class dict
    """
    def __init__(self, fid,  pointer):
        # block header
        if pointer != 0 and not pointer == None:
            fid.seek( pointer )
            self['id']=self.mdfblockreadCHAR(fid, 4)
            self['reserved']=self.mdfblockreadBYTE(fid, 4)
            self['length']=self.mdfblockread(fid, UINT64, 1)
            self['link_count']=self.mdfblockread(fid, UINT64, 1)
            # link section
            self['sr_sr_next']=self.mdfblockread(fid, LINK, 1)
            self['sr_data']=self.mdfblockread(fid, LINK, 1)
            # data section
            self['sr_cycle_count']=self.mdfblockread(fid, UINT64, 1)
            self['sr_interval']=self.mdfblockread(fid, REAL, 1)
            self['sr_sync_type']=self.mdfblockread(fid, UINT8, 1)
            self['sr_flags']=self.mdfblockread(fid, UINT8, 1)
            self['sr_reserved']=self.mdfblockreadBYTE(fid, 6)
            
class SIBlock(MDFBlock):
    """ reads Source Information block and saves in class dict
    """
    def __init__(self, fid,  pointer):
        # block header
        if pointer != 0 and not pointer == None:
            fid.seek( pointer )
            self['id']=self.mdfblockreadCHAR(fid, 4)
            self['reserved']=self.mdfblockreadBYTE(fid, 4)
            self['length']=self.mdfblockread(fid, UINT64, 1)
            self['link_count']=self.mdfblockread(fid, UINT64, 1)
            # link section
            self['si_tx_name']=self.mdfblockread(fid, LINK, 1)
            self['si_tx_path']=self.mdfblockread(fid, LINK, 1)
            self['si_md_comment']=self.mdfblockread(fid, LINK, 1)
            # data section
            self['si_type']=self.mdfblockread(fid, UINT8, 1)
            if self['si_type']==0:
                self['si_type']='OTHER' # unknown
            elif self['si_type']==1:
                self['si_type']='ECU'
            elif self['si_type']==2:
                self['si_type']='BUS'
            elif self['si_type']==3:
                self['si_type']='I/O'
            elif self['si_type']==4:
                self['si_type']='TOOL'
            elif self['si_type']==5:
                self['si_type']='USER'
            self['si_bus_type']=self.mdfblockread(fid, UINT8, 1)
            if self['si_bus_type']==0:
                self['si_bus_type']='NONE'
            elif self['si_bus_type']==1:
                self['si_bus_type']='OTHER'
            elif self['si_bus_type']==2:
                self['si_bus_type']='CAN'
            elif self['si_bus_type']==3:
                self['si_bus_type']='LIN'
            elif self['si_bus_type']==4:
                self['si_bus_type']='MOST'
            elif self['si_bus_type']==5:
                self['si_bus_type']='FLEXRAY'
            elif self['si_bus_type']==6:
                self['si_bus_type']='K_LINE'
            elif self['si_bus_type']==7:
                self['si_bus_type']='ETHERNET'
            elif self['si_bus_type']==8:
                self['si_bus_type']='USB'
            self['si_flags']=self.mdfblockread(fid, UINT8, 1)
            self['si_reserved']=self.mdfblockreadBYTE(fid, 5)
            # post treatment
            self['source_name']=CommentBlock(fid, self['si_tx_name'])
            self['source_path']=CommentBlock(fid, self['si_tx_path'])
            self['comment']=CommentBlock(fid, self['si_md_comment'], MDType='SI')

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
        self['IDBlock'] = {} # Identifier Block
        self['HDBlock'] = {} # Header Block
        self['FHBlock'] = {}
        self['CHBlock'] = {}
        self['DGBlock']= {} # Data Group Block
        self['CGBlock'] = {} # Channel Group Block
        self['CNBlock'] = {}# Channel Block
        self['CCBlock'] = {} # Conversion block
        self['ATBlock'] = {} # Attachment block
        self.fileName = fileName
        if fileName != None and fid==None:
            try:
                fid = open( self.fileName, 'rb' )
            except IOError:
                print('Can not find file '+self.fileName)
                raise
            self.readinfo( fid )
            # Close the file
            fid.close()
        elif fileName == None and fid!=None:
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
        
        #print('reads File History blocks, always exists')
        fh=0 # index of fh blocks
        self['FHBlock'][fh] = {}
        self['FHBlock'][fh] .update(FHBlock(fid, self['HDBlock']['hd_fh_first']))
        while self['FHBlock'][fh]['fh_fh_next']:
            self['FHBlock'][fh+1]={}
            self['FHBlock'][fh+1] .update(FHBlock(fid, self['FHBlock'][fh]['fh_fh_next']))
            fh+=1
        
        # print('reads Channel Hierarchy blocks')
        if self['HDBlock']['hd_ch_first']:
            ch=0
            self['CHBlock'][ch] = {}
            self['CHBlock'][ch] .update(CHBlock(fid, self['HDBlock']['hd_ch_first']))
            while self['CHBlock'][ch]['ch_ch_next']:
                self['CHBlock'][ch] .update(CHBlock(fid, self['CHBlock'][ch]['ch_ch_next']))
                ch+=1
        
        # reads Attachment block
        self['ATBlock']=self.readATBlock(fid, self['HDBlock']['hd_at_first'])
        
        # reads Event Block
        if self['HDBlock']['hd_ev_first']:
            ev=0
            self['EVBlock']={}
            self['EVBlock'][ev]=EVBlock(fid, self['HDBlock']['hd_ev_first'])
            while self['EVBlock'][ev]['ev_ev_next']:
                ev+=1
                self['EVBlock'][ev]=EVBlock(fid, self['EVBlock'][ev-1]['ev_ev_next'])
        
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
        if self['HDBlock']['hd_dg_first']:
            dg=0
            self['DGBlock'][dg] = {}
            self['DGBlock'][dg].update(DGBlock(fid,self['HDBlock']['hd_dg_first']))
            # reads Channel Group blocks
            self.readCGBlock(fid, dg, channelNameList)
            while self['DGBlock'][dg]['dg_dg_next']:
                dg+=1
                self['DGBlock'][dg]={}
                self['DGBlock'][dg].update(DGBlock(fid, self['DGBlock'][dg-1]['dg_dg_next']))
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
        if self['DGBlock'][dg] ['dg_cg_first']:
            cg=0
            self['CNBlock'][dg] = {}
            self['CNBlock'][dg][cg] = {}
            self['CCBlock'][dg] = {}
            self['CCBlock'][dg][cg] = {}
            self['CGBlock'][dg] = {}
            self['CGBlock'][dg] [cg]={}
            self['CGBlock'][dg][cg].update(CGBlock(fid, self['DGBlock'][dg]['dg_cg_first']))
            VLSDCGBlock=[]
            
            if not channelNameList:
                # reads Source Information Block
                self['CGBlock'][dg][cg]['SIBlock']=SIBlock(fid, self['CGBlock'][dg][cg]['cg_si_acq_source'])
                
                # reads Sample Reduction Block
                self['CGBlock'][dg][cg]['SRBlock']=self.readSRBlock(fid, self['CGBlock'][dg][cg]['cg_sr_first'])
            
            if not self['CGBlock'][dg][cg]['cg_flags'] & 0b1: # if not a VLSD channel group
                # reads Channel Block
                self.readCNBlock(fid, dg, cg, channelNameList)
            else:
                VLSDCGBlock.append(cg)

            while self['CGBlock'][dg][cg]['cg_cg_next']:
                cg+=1
                self['CGBlock'][dg][cg]={}
                self['CGBlock'][dg][cg].update(CGBlock(fid, self['CGBlock'][dg][cg-1]['cg_cg_next']))
                self['CNBlock'][dg][cg] = {}
                self['CCBlock'][dg][cg] = {}
                if not channelNameList:
                    # reads Source Information Block
                    self['CGBlock'][dg][cg]['SIBlock']=SIBlock(fid, self['CGBlock'][dg][cg]['cg_si_acq_source'])
                    
                    # reads Sample Reduction Block
                    self['CGBlock'][dg][cg]['SRBlock']=self.readSRBlock(fid, self['CGBlock'][dg][cg]['cg_sr_first'])
                    
                if not self['CGBlock'][dg][cg]['cg_flags'] & 0b1: # if not a VLSD channel group
                    # reads Channel Block
                    self.readCNBlock(fid, dg, cg, channelNameList)
                else:
                    VLSDCGBlock.append(cg)
            
            if VLSDCGBlock: #VLSD CG Block exiting
                self['VLSD_CG']={}
            # Matching VLSD CGBlock with corresponding channel
            for VLSDcg in VLSDCGBlock:
                VLSDCGBlockAdress=self['CGBlock'][dg][VLSDcg]['pointer']
                for cg in self['CGBlock'][dg]:
                    if cg not in VLSDCGBlock:
                        for cn in self['CNBlock'][dg][cg]:
                            if VLSDCGBlockAdress==self['CNBlock'][dg][cg][cn]['cn_data']:
                                # found matching channel with VLSD CGBlock
                                temp={}
                                temp['cg_cn']=(cg, cn)
                                self['VLSD_CG'][self['CGBlock'][dg][VLSDcg]['cg_record_id']]=temp
                                break
            
              # reorder channel blocks and related blocks(CC, SI, AT, CA) based on byte offset
              # this reorder is meant to improve performance while parsing records using core.records.fromfile
              # as it will not use cn_byte_offset
              # first, calculate new mapping/order
            nChannel=len(self['CNBlock'][dg][cg])
            Map = zeros(shape=len(self['CNBlock'][dg][cg]), dtype=[('index','u4'), ('byte_offset', 'u4')])
            for cn in range(nChannel):
                Map[cn]=(cn, self['CNBlock'][dg][cg][cn]['cn_byte_offset'])
            orderedMap=sort(Map, order='byte_offset')
            
            toChangeIndex=Map==orderedMap
            for cn in range(nChannel):
                if not toChangeIndex[cn]:
                    # offset all indexes of indexes to be moved
                    self['CNBlock'][dg][cg][cn+nChannel]=self['CNBlock'][dg][cg].pop(cn)
                    self['CCBlock'][dg][cg][cn+nChannel]=self['CCBlock'][dg][cg].pop(cn)
            for cn in range(nChannel):
                if not toChangeIndex[cn]:
                    # change to ordered index
                    self['CNBlock'][dg][cg][cn]=self['CNBlock'][dg][cg].pop(orderedMap[cn][0]+nChannel)
                    self['CCBlock'][dg][cg][cn]=self['CCBlock'][dg][cg].pop(orderedMap[cn][0]+nChannel)

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
        cn=0
        self['CNBlock'][dg][cg][cn] = {}
        self['CCBlock'][dg][cg][cn] = {}
        self['CNBlock'][dg][cg][cn]=CNBlock(fid, self['CGBlock'][dg][cg]['cg_cn_first'])
        MLSDChannels=[]
        # check for MLSD
        if self['CNBlock'][dg][cg][cn]['cn_type']==5:
            MLSDChannels.append(cn)
        # check if already existing channel name
        for cgp in self['CNBlock'][dg]:
            for chan in self['CNBlock'][dg][cgp]:
                if not chan==cn and self['CNBlock'][dg][cgp][chan]['name']==self['CNBlock'][dg][cg][cn]['name']: 
                    self['CNBlock'][dg][cg][cn]['name']=self['CNBlock'][dg][cg][cn]['name']+str(cg)
                    break
        
        if self['CGBlock'][dg][cg]['cg_cn_first']: # Can be NIL for VLSD
            if not channelNameList:
                # reads Channel Source Information
                self['CNBlock'][dg][cg][cn]['SIBlock']=SIBlock(fid, self['CNBlock'][dg][cg][cn]['cn_si_source'])
                
                # reads Channel Array Block
                if self['CNBlock'][dg][cg][cn]['cn_composition']: # composition but can be either structure of channels or array
                    fid.seek(self['CNBlock'][dg][cg][cn]['cn_composition'])
                    if fid.read(4) in ('##CA',b'##CA'):
                        self['CNBlock'][dg][cg][cn]['CABlock']=CABlock(fid, self['CNBlock'][dg][cg][cn]['cn_composition'])
                    elif self.mdfblockreadCHAR(fid, 4) in ('##CN',b'##CN'):
                        self['CNBlock'][dg][cg][cn]['CNBlock']=CNBlock(fid, self['CNBlock'][dg][cg][cn]['cn_composition'])
                    else:
                        raise('unknown channel composition')
                        
                # reads Attachment Block
                if self['CNBlock'][dg][cg][cn]['cn_attachment_count']>1:
                    for at in range(self['CNBlock'][dg][cg][cn]['cn_attachment_count']):
                        self['CNBlock'][dg][cg][cn]['attachment'][at].update(self.readATBlock(fid, self['CNBlock'][dg][cg][cn]['cn_at_reference'][at]))
                elif self['CNBlock'][dg][cg][cn]['cn_attachment_count']==1:
                    self['CNBlock'][dg][cg][cn]['attachment'][0].update(self.readATBlock(fid, self['CNBlock'][dg][cg][cn]['cn_at_reference']))
                
                # reads Channel Conversion Block
                self['CCBlock'][dg][cg][cn]=CCBlock(fid, self['CNBlock'][dg][cg][cn]['cn_cc_conversion'])
            while self['CNBlock'][dg][cg][cn]['cn_cn_next']:
                cn=cn+1
                self['CNBlock'][dg][cg][cn]=CNBlock(fid, self['CNBlock'][dg][cg][cn-1]['cn_cn_next'])
                # check for MLSD
                if self['CNBlock'][dg][cg][cn]['cn_type']==5:
                    MLSDChannels.append(cn)
                if not channelNameList:
                    # reads Channel Source Information
                    self['CNBlock'][dg][cg][cn]['SIBlock']=SIBlock(fid, self['CNBlock'][dg][cg][cn]['cn_si_source'])
                    
                    # check if already existing channel name
                    for cgp in self['CNBlock'][dg]:
                        for chan in self['CNBlock'][dg][cgp]:
                            if not chan==cn and self['CNBlock'][dg][cgp][chan]['name']==self['CNBlock'][dg][cg][cn]['name']: 
                                self['CNBlock'][dg][cg][cn]['name']=self['CNBlock'][dg][cg][cn]['name']+str(cg)
                                break
                                
                    # reads Channel Array Block
                    if self['CNBlock'][dg][cg][cn]['cn_composition']: # composition but can be either structure of channels or array
                        fid.seek(self['CNBlock'][dg][cg][cn]['cn_composition'])
                        id=fid.read(4)
                        if id in ('##CA',b'##CA'):
                            self['CNBlock'][dg][cg][cn]['CABlock']=CABlock(fid, self['CNBlock'][dg][cg][cn]['cn_composition'])
                        elif id in ('##CN',b'##CN'):
                            self['CNBlock'][dg][cg][cn]['CNBlock']=CNBlock(fid, self['CNBlock'][dg][cg][cn]['cn_composition'])
                        else:
                            raise('unknown channel composition')
                        
                    # reads Attachment Block
                    if self['CNBlock'][dg][cg][cn]['cn_attachment_count']>1:
                        for at in range(self['CNBlock'][dg][cg][cn]['cn_attachment_count']):
                            print(self['CNBlock'][dg][cg][cn]['cn_at_reference'][at])
                            self['CNBlock'][dg][cg][cn]['attachment'][at].update(self.readATBlock(fid, self['CNBlock'][dg][cg][cn]['cn_at_reference'][at][0]))
                    elif self['CNBlock'][dg][cg][cn]['cn_attachment_count']==1:
                        self['CNBlock'][dg][cg][cn]['attachment'][0].update(self.readATBlock(fid, self['CNBlock'][dg][cg][cn]['cn_at_reference']))

                    # reads Channel Conversion Block
                    self['CCBlock'][dg][cg][cn]=CCBlock(fid, self['CNBlock'][dg][cg][cn]['cn_cc_conversion'])
        
        MLSDChannels=self.readComposition(fid, dg, cg, MLSDChannels, channelNameList=False)

        if MLSDChannels:
            self['MLSD']={}
            self['MLSD'][dg]={}
            self['MLSD'][dg][cg]={}
        for MLSDcn in MLSDChannels:
            for cn in self['CNBlock'][dg][cg]:
                if self['CNBlock'][dg][cg][cn]['pointer']==self['CNBlock'][dg][cg][MLSDcn]['cn_data']:
                    self['MLSD'][dg][cg][MLSDcn]=cn
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
        chan=list(self['CNBlock'][dg][cg].keys())[-1]+1
        for cn in list(self['CNBlock'][dg][cg].keys()):
            if self['CNBlock'][dg][cg][cn]['cn_composition']:
                fid.seek(self['CNBlock'][dg][cg][cn]['cn_composition'])
                ID=MDFBlock()
                ID=ID.mdfblockreadCHAR(fid, 4)
                if ID in ('##CN',b'##CN'): # Structures
                    self['CNBlock'][dg][cg][chan]=CNBlock(fid, self['CNBlock'][dg][cg][cn]['cn_composition'])
                    self['CCBlock'][dg][cg][chan]=CCBlock(fid, self['CNBlock'][dg][cg][chan]['cn_cc_conversion'])
                    if self['CNBlock'][dg][cg][chan]['cn_type']==5:
                        MLSDChannels.append(chan)
                    while self['CNBlock'][dg][cg][chan]['cn_cn_next']:
                        chan+=1
                        self['CNBlock'][dg][cg][chan]=CNBlock(fid, self['CNBlock'][dg][cg][chan-1]['cn_cn_next'])
                        self['CCBlock'][dg][cg][chan]=CCBlock(fid, self['CNBlock'][dg][cg][chan]['cn_cc_conversion'])
                        if self['CNBlock'][dg][cg][chan]['cn_type']==5:
                            MLSDChannels.append(chan)
                    self['CNBlock'][dg][cg][cn]['cn_type']=6 # makes the channel virtual
                elif ID in ('##CA',b'##CA'): # arrays
                    pass
                else:
                    print('unknown channel composition')
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
        if pointer>0:
            sr=0
            srBlocks={}
            srBlocks[0]=SRBlock(fid, pointer)
            while srBlocks[sr]['sr_sr_next']>0:
                sr+=1
                srBlocks[sr]=SRBlock(fid, srBlocks[sr-1]['sr_sr_next'])
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
        if pointer >0:
            at=0
            atBlocks={}
            if type(pointer) in (tuple, list):
                pointer=pointer[0]
            atBlocks[0]=ATBlock(fid, pointer)
            while atBlocks[at]['at_at_next']>0:
                at+=1
                atBlocks[at]=(ATBlock(fid, atBlocks[at-1]['at_at_next']))
            return atBlocks
    
    def listChannels4( self, fileName = None ):
        """ Read MDF file and extract its complete structure
        
        Parameters
        ----------------
        fileName : str
            file name
        
        Returns
        -----------
        list of channel names contained in file
        """
        if not fileName == None:
            self.fileName = fileName
        # Open file
        fid = open( self.fileName, 'rb' )
        channelNameList=[]
        # reads Header HDBlock
        self['HDBlock'].update(HDBlock(fid))
            
        # reads Data Group, channel groups and channel Blocks  recursively but not the other metadata block
        self.readDGBlock(fid, True)
        
        for dg in list(self['DGBlock'].keys()):
            for cg in list(self['CGBlock'][dg].keys()):
                for cn in list(self['CNBlock'][dg][cg].keys()):
                    channelNameList.append(self['CNBlock'][dg][cg][cn]['name'])
        
        # CLose the file
        fid.close()
        return channelNameList

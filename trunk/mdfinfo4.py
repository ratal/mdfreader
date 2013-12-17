# -*- coding: utf-8 -*-
""" Measured Data Format blocks paser for version 4.x
Created on Sun Dec 15 12:57:28 2013

:Author: `Aymeric Rateau <http://code.google.com/p/mdfreader/>`__

"""

from struct import calcsize, unpack
from sys import version_info
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
    
    **Methods**
    
    *mdfblockread* : converts a byte array to a given data type
    
    **Attributes**
    
    none        
     
    """
    
    def __init__(self):
        pass
    
    def loadHeader(self, fid,  pointer):
        #All blocks have the same header
        if pointer != 0 and not pointer == None:
            fid.seek( pointer )
            self['id']=self.mdfblockreadCHAR(fid, 4)
            self['reserved']=self.mdfblockreadBYTE(fid, 4)
            self['length']=self.mdfblockread(fid, UINT64, 1)
            self['link_count']=self.mdfblockread(fid, UINT64, 1)
        
    @staticmethod
    def mdfblockread( fid, type,  count ):
        # reads value in file corresponding to format and size
        value=fid.read(calcsize(type)*count)
        if value:
            value=unpack( type,value)
            value=value[0]
        else:
            value=None
        return value
        
    @staticmethod
    def mdfblockreadCHAR( fid,  count ):
        # reads value in file corresponding to format and size
        value=fid.read( count )
        if PythonVersion<3:
            value.replace('\x00', '')
        else:
            value.decode('iso8859-1', 'replace').replace('\x00', '') 
        return value
        
    @staticmethod
    def mdfblockreadBYTE( fid, count ):
        # reads value in file corresponding to format and size
        value=fid.read( count )
        value.encode('UTF-8').replace('\x00', '') 
        return value
        
class IDBlock(MDFBlock):
    def __init__(self, fid):
    # reads IDBlock
        fid.seek(0)
        self['id_file']=self.mdfblockreadCHAR(fid, 8)
        self['id_vers']=self.mdfblockreadCHAR(fid, 8)
        self['id_prog']=self.mdfblockreadCHAR(fid, 8)
        self['id_reserved1']=self.mdfblockreadBYTE(fid, 4)
        self['id_ver']=self.mdfblockread(fid, UINT16, 1)
        self['id_reserved1']=self.mdfblockreadBYTE(fid, 34)

class HDBlock(MDFBlock):
    def __init__(self, fid,  pointer=64):
        # block header
        self.loadHeader(fid, pointer)
        # header block
        self['hd_dg_first']=self.mdfblockread(fid, LINK, 1)
        self['hd_fh_first']=self.mdfblockread(fid, LINK, 1)
        self['hd_ch_first']=self.mdfblockread(fid, LINK, 1)
        self['hd_at_first']=self.mdfblockread(fid, LINK, 1)
        self['hd_ev_first']=self.mdfblockread(fid, LINK, 1)
        self['hd_start_time_ns']=self.mdfblockread(fid, UINT64, 1)
        self['hd_tz_offset_min']=self.mdfblockread(fid, INT16, 1)
        self['hd_dst_offset_min']=self.mdfblockread(fid, INT16, 1)
        self['hd_time_flags']=self.mdfblockread(fid, UINT8, 1)
        self['hd_time_class']=self.mdfblockread(fid, UINT8, 1)
        self['hd_flags']=self.mdfblockread(fid, UINT8, 1)
        self['hd_reserved']=self.mdfblockreadBYTE(fid, 1)
        self['hd_start_angle_rad']=self.mdfblockread(fid, REAL, 1)
        self['hd_start_distance']=self.mdfblockread(fid, REAL, 1)

class MDBlock(MDFBlock):
    def __init__(self, fid,  pointer):
        # block header
        self.loadHeader(fid, pointer)
        # Metadata block
        self['md_data']=self.mdfblockreadBYTE(fid, self['length'])
        
class TXBlock(MDFBlock):
    def __init__(self, fid,  pointer):
        # block header
        self.loadHeader(fid, pointer)
        # Text block
        self['tx_data']=self.mdfblockreadBYTE(fid, self['length'])
        
class DGBlock(MDFBlock):
    def __init__(self, fid,  pointer):
        # block header
        self.loadHeader(fid, pointer)
        # Data Group block : Links
        self['dg_dg_next']=self.mdfblockread(fid, LINK, 1)
        self['dg_cg_next']=self.mdfblockread(fid, LINK, 1)
        self['dg_data']=self.mdfblockread(fid, LINK, 1)
        self['dg_md_comment']=self.mdfblockread(fid, LINK, 1)
        # data section
        self['dg_rec_id_size']=self.mdfblockread(fid, UINT8, 1)
        self['dg_reserved']=self.mdfblockreadBYTE(fid, 7)
        
class CGBlock(MDFBlock):
    def __init__(self, fid,  pointer):
        # block header
        self.loadHeader(fid, pointer)
        # Channel Group block : Links
        self['cg_cg_next']=self.mdfblockread(fid, LINK, 1)
        self['cg_cn_next']=self.mdfblockread(fid, LINK, 1)
        self['cg_tx_acq_name']=self.mdfblockread(fid, LINK, 1)
        self['cg_si_acq_source']=self.mdfblockread(fid, LINK, 1)
        self['cg_sr_first']=self.mdfblockread(fid, LINK, 1)
        self['dg_md_comment']=self.mdfblockread(fid, LINK, 1)
        # data section
        self['cg_record_id']=self.mdfblockread(fid, UINT64, 1)
        self['cg_cycle_count']=self.mdfblockread(fid, UINT64, 1)
        self['cg_flags']=self.mdfblockread(fid, UINT16, 1)
        self['cg_path_separator']=self.mdfblockread(fid, UINT16, 1)
        self['cg_reserved']=self.mdfblockreadBYTE(fid, 4)
        self['cg_data_bytes']=self.mdfblockread(fid, UINT32, 1)
        self['cg_invalid_bytes']=self.mdfblockread(fid, UINT32, 1)
        
class CNBlock(MDFBlock):
    def __init__(self, fid,  pointer,  attachmentNumber=None):
        # block header
        self.loadHeader(fid, pointer)
        # Channel Group block : Links
        self['cn_cn_next']=self.mdfblockread(fid, LINK, 1)
        self['cn_composition']=self.mdfblockread(fid, LINK, 1)
        self['cn_tx_name']=self.mdfblockread(fid, LINK, 1)
        self['cn_si_source']=self.mdfblockread(fid, LINK, 1)
        self['cn_cc_conversion']=self.mdfblockread(fid, LINK, 1)
        self['cn_data']=self.mdfblockread(fid, LINK, 1)
        self['cn_mdunit']=self.mdfblockread(fid, LINK, 1)
        self['cn_md_comment']=self.mdfblockread(fid, LINK, 1)
        self['cn_at_reference']=self.mdfblockread(fid, LINK, attachmentNumber)
        self['cn_default_x']=self.mdfblockread(fid, LINK, 3)
        # data section
        self['cn_type']=self.mdfblockread(fid, UINT8, 1)
        self['cn_sync_type']=self.mdfblockread(fid, UINT8, 1)
        self['cn_datatype']=self.mdfblockread(fid, UINT8, 1)
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

class CCBlock(MDFBlock):
    def __init__(self, fid,  pointer,  attachmentNumber=None):
        # block header
        self.loadHeader(fid, pointer)

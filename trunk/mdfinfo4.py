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
        self['id']={}
        self['reserved']={}
        self['length']={}
        self['link_count']={}
    
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
        # UTF-8 encoded bytes
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
        self['id_reserved2']=self.mdfblockreadBYTE(fid, 34)

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
            self['Comment']=MDBlock(fid, self['hd_md_comment'])

class FHBlock(MDFBlock):
    def __init__(self, fid,  pointer):
        # block header
        self.loadHeader(fid, pointer)
        # history block
        self['fh_fh_next']=self.mdfblockread(fid, LINK, 1)
        self['fh_md_comment']=self.mdfblockread(fid, LINK, 1)
        self['fh_time_ns']=self.mdfblockread(fid, UINT64, 1)
        self['fh_tz_offset_min']=self.mdfblockread(fid, UINT16, 1)
        self['fh_dst_offset_min']=self.mdfblockread(fid, UINT16, 1)
        self['fh_time_flags']=self.mdfblockread(fid, UINT8, 1)
        self['fh_reserved']=self.mdfblockreadBYTE(fid, 3)
        if self['fh_md_comment']: # comments exist
            self['Comment']=MDBlock(fid, self['fh_md_comment'])
            
class CHBlock(MDFBlock):
    def __init__(self, fid,  pointer):
        # block header
        self.loadHeader(fid, pointer)
        # Channel hierarchy block
        self['ch_ch_next']=self.mdfblockread(fid, LINK, 1)
        self['ch_ch_first']=self.mdfblockread(fid, LINK, 1)
        self['ch_tx_name']=self.mdfblockread(fid, LINK, 1)
        self['ch_md_comment']=self.mdfblockread(fid, LINK, 1)
        self['ch_element']=self.mdfblockread(fid, LINK, (self['link_count']-4)/3)
        self['ch_element_count']=self.mdfblockread(fid, LINK, 1)
        self['ch_type']=self.mdfblockread(fid, UINT8, 1)
        self['ch_reserved']=self.mdfblockreadBYTE(fid, 3)
        if self['ch_md_comment']: # comments exist
            self['Comment']=MDBlock(fid, self['ch_md_comment'])
        if self['ch_tx_name']: # text block containing name of hierarchy level
            self['ch_name_level']=TXBlock(fid, self['ch_tx_name'])

class MDBlock(MDFBlock):
    def __init__(self, fid,  pointer):
        # block header
        self.loadHeader(fid, pointer)
        # Metadata block
        self['md_data']=self.mdfblockreadBYTE(fid, self['length']-16)
        
class TXBlock(MDFBlock):
    def __init__(self, fid,  pointer):
        # block header
        self.loadHeader(fid, pointer)
        # Text block
        self['tx_data']=self.mdfblockreadBYTE(fid, self['length']-16)
        
class DGBlock(MDFBlock):
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
            self['Comment']=MDBlock(fid, self['dg_md_comment'])
        
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
            self['Comment']=MDBlock(fid, self['cg_md_comment'])
        if self['cg_tx_acq_name']: # comments exist
            self['acq_name']=TXBlock(fid, self['cg_tx_acq_name'])
        
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
        self['cn_md_unit']=self.mdfblockread(fid, LINK, 1)
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
        if self['cn_md_comment']: # comments exist
            self['Comment']=MDBlock(fid, self['cn_md_comment'])
        if self['cn_md_unit']: # comments exist
            self['unit']=MDBlock(fid, self['cn_md_unit'])
        if self['cn_tx_name']: # comments exist
            self['name']=TXBlock(fid, self['cn_tx_name'])

class CCBlock(MDFBlock):
    def __init__(self, fid,  pointer):
        # block header
        self.loadHeader(fid, pointer)
        # Channel Conversion block : Links
        self['cc_tx_name']=self.mdfblockread(fid, LINK, 1)
        self['cc_cn_next']=self.mdfblockread(fid, LINK, 1)
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
        if self['cc_md_comment']: # comments exist
            self['Comment']=MDBlock(fid, self['cc_md_comment'])
        if self['cc_md_unit']: # comments exist
            self['unit']=MDBlock(fid, self['cc_md_unit'])
        if self['cc_tx_name']: # comments exist
            self['name']=TXBlock(fid, self['cc_tx_name'])
        
class DTBlock(MDFBlock):
    def __init__(self, fid,  pointer):
        # block header
        self.loadHeader(fid, pointer)
        # Text block
        #self['dt_data']=self.mdfblockreadBYTE(fid, self['link_count']-16)
        
class info4(dict):
    def __init__(self, fileName=None, fid=None):
        self['IDBlock'] = {} # Identifier Block
        self['HDBlock'] = {} # Header Block
        self['FHBlock'] = {}
        self['CHBlock'] = {}
        self['DGBlock']= {} # Data Group Block
        self['CGBlock'] = {} # Channel Group Block
        self['CNBlock'] = {}# Channel Block
        self['CCBlock'] = {} # Conversion block
        self.fileName = fileName
        if fileName != None and fid==None:
            try:
                fid = open( self.fileName, 'rb' )
            except IOError:
                print('Can not find file'+self.fileName)
                raise
            self.readinfo( fid )
            # CLose the file
            fid.close()
        elif fileName == None and fid!=None:
            self.readinfo(fid)

    def readinfo(self, fid):
        # this function reads all blocks contained in mdf4 file
        
        # reads IDBlock
        self['IDBlock'].update(IDBlock(fid))
        
        # reads Header HDBlock
        self['HDBlock'].update(HDBlock(fid))        
        
        # reads File History blocks, always exists
        fh=0 # index of fh blocks
        self['FHBlock'][fh] = {}
        self['FHBlock'][fh] .update(FHBlock(fid, self['HDBlock']['hd_fh_first']))
        while self['FHBlock'][fh]['fh_fh_next']:
            self['FHBlock'][fh] .update(FHBlock(fid, self['FHBlock'][fh]['fh_fh_next']))
            fh+=1
        
        # reads Channel Hierarchy blocks
        if self['HDBlock']['hd_ch_first']:
            ch=0
            self['CHBlock'][ch] = {}
            self['CHBlock'][ch] .update(CHBlock(fid, self['HDBlock']['hd_ch_first']))
            while self['CHBlock'][ch]['ch_ch_next']:
                self['CHBlock'][ch] .update(CHBlock(fid, self['CHBlock'][ch]['ch_ch_next']))
                ch+=1
        
        # reads Attachment block
        
        # reads Event Block
        
        # reads Data Group Blocks
        if self['HDBlock']['hd_dg_first']:
            dg=0
            self['DGBlock'][dg] = {}
            self['DGBlock'][dg] .update(DGBlock(fid,self['HDBlock']['hd_dg_first']))
            while self['DGBlock'][dg]['dg_dg_next']:
                dg+=1
                self['DGBlock'][dg]={}
                self['DGBlock'][dg].update(DGBlock(fid, self['DGBlock'][dg-1]['dg_dg_next']))
                # reads Channel Group blocks
                if self['DGBlock'][dg] ['dg_cg_first']:
                    cg=0
                    self['CGBlock'][dg] = {}
                    self['CGBlock'][dg] [cg]={}
                    self['CGBlock'][dg][cg].update(CGBlock(fid, self['DGBlock'][dg]['dg_cg_first']))
                    while self['CGBlock'][dg][cg]['cg_cg_next']:
                        cg+=1
                        self['CGBlock'][dg][cg]={}
                        self['CGBlock'][dg][cg].update(CGBlock(fid, self['CGBlock'][dg][cg-1]['cg_cg_next']))
                        # reads Source Information Block
                        
                        # reads Sample Reduction Block
                        
                        # reads Channel Block
                        cn=0
                        self['CNBlock'][dg] = {}
                        self['CNBlock'][dg][cg] = {}
                        self['CNBlock'][dg][cg][cn] = {}
                        self['CNBlock'][dg][cg][cn].update(CNBlock(fid, self['CGBlock'][dg][cg]['cg_cn_first']))
                        while self['CNBlock'][dg][cg][cn]['cn_cn_next']:
                            cn+=1
                            self['CNBlock'][dg][cg][cn]={}
                            self['CNBlock'][dg][cg][cn].update(CNBlock(fid, self['CNBlock'][dg][cg][cn-1]['cn_cn_next']))
                            # reads Channel Source Information
                            
                            # reads Channel Array Block
                            
                            # reads Attachment Block
                            
                            # reads Channel Conversion Block
                            self['CCBlock'][dg][cg][cn]=CNBlock(fid, self['CNBlock'][dg][cg][cn]['cn_cc_conversion'])
        
    def listChannels( self, fileName = None ):
        # Read MDF file and extract its complete structure
        if not fileName == None:
            self.fileName = fileName
        # Open file
        fid = open( self.fileName, 'rb' )
        # CLose the file
        fid.close()
        #return channelNameList

# -*- coding: utf-8 -*-
''' Veusz data plugin for mdf format

Created on 18 nov. 2010

@author: Aymeric Rateau ; aymeric.rateau@gmail.com

Veusz version confirmed : from 1.16

Installation instructions :
----------------------------------
1. copy or link mdfForVeuszPlugin.py into the plugins directory of Veusz
2. link mdfreader.py in root or plugin directory of Veusz
3. Go in Veusz ; Edit / Preferences / Plugins tab and add mdfForVeuszPlugins.py
'''
from veusz.plugins import *
from veusz.plugins.datasetplugin import Dataset1D as ImportDataset1D
from veusz.plugins.field import FieldFloat as ImportFieldFloat
try:
    from mdfreader import mdf
    from mdfreader import mdfinfo
except ImportError:
    try:
        from veusz.plugins.mdfreader import mdf
        from veusz.plugins.mdfreader import mdfinfo
    except ImportError:
        import sys
        import os
        import os.path
        sys.path.append( os.getcwdu() )
        print( os.path.join( os.getcwdu(), 'plugins' ) )
        from mdfreader import mdf
        from mdfreader import mdfinfo

class ImportPlugin( mdfinfo ):
    """Define a plugin to read data in a particular format.
    
    override doImport and optionally getPreview to define a new plugin
    register the class by adding to the importpluginregistry list
    """

    name = 'Import plugin'
    author = 'Aymeric Rateau'
    description = 'Import MDF files'
    promote_tab = 'MDF'
    file_extensions = set(['.dat', '.mf4', '.mdf'])

    def __init__( self ):
        """Override this to declare a list of input fields if required."""
        # a list of ImportField objects to display
        self.fields = []

    def getPreview( self, params ):
        """Get data to show in a text box to show a preview.
        params is a ImportPluginParams object.
        Returns (text, okaytoimport)
        """
        
        info = mdfinfo( fileName=params.filename )
        
        if info.mdfversion < 400:
            f = ''
            f += 'Time: ' + info['HDBlock']['Date'] + ' '
            f += info['HDBlock']['Time'] + '\n'
            f += 'Author: ' + info['HDBlock']['Author'] + '\n'
            f += 'Organisation: ' + info['HDBlock']['Organization' ] + '\n'
            f += 'Project Name: ' + info['HDBlock']['ProjectName'] + '\n'
            f += 'Subject: ' + info['HDBlock']['Subject'] + '\n' + 'Channel List:\n'
        else:
            from time import gmtime, strftime
            fileDateTime=gmtime(info['HDBlock']['hd_start_time_ns']/1000000000)
            date=strftime('%Y-%m-%d', fileDateTime)
            time=strftime('%H:%M:%S', fileDateTime)
            f = ''
            f += 'Date Time: ' + date + '  ' + time + '\n'
            if 'Comment' in info['HDBlock']:
                Comment = info['HDBlock']['Comment']
                if 'author' in Comment:
                    f += 'Author: ' + Comment['author'] + '\n'
                if 'department' in Comment:
                    f += 'Organisation: ' + Comment['department' ] + '\n'
                if 'project' in Comment:
                    f += 'Project Name: ' + Comment['project'] + '\n'
                if 'subject' in Comment:
                    f += 'Subject: ' + Comment['subject'] + '\n' + 'Channel List:\n'
        for channelName in info.listChannels():
            f += '   ' + channelName + '\n'
        return f, True

    def doImport( self, params ):
        """Actually import data
        params is a ImportPluginParams object.
        Return a list of ImportDataset1D, ImportDataset2D objects
        """
        return []

class MdfImportPlugin( ImportPlugin, mdf ):
    """Plugin to import mdf (Mostly ETAS INCA or CANape files)"""

    name = "MDFImport plugin"
    author = "Aymeric Rateau"
    description = "Reads MDF files from INCA or CANAPE"

    def __init__( self ):
        ImportPlugin.__init__( self )
        self.fields = [ImportFieldFloat( "mult", descr = "Sampling", default = 0.1 )]

    def doImport( self, params ):
        """Actually import data
        params is a ImportPluginParams object.
        Return a list of ImportDataset1D, ImportDataset2D objects
        """

        data = mdf( params.filename )
        data.resample( samplingTime = params.field_results['mult'] )
        List = []
        for channelName in list(data.keys()):
            if len( data[channelName]['data'] ) > 0 and not data[channelName]['data'].dtype.kind in ['S', 'U']:
                #print( data[channelName]['data'].dtype )
                List.append( ImportDataset1D( channelName, data[channelName]['data'] ) )
        return List

importpluginregistry.append( MdfImportPlugin() )

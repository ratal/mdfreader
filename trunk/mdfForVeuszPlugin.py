'''
Created on 18 nov. 2010

@author: Aymeric Rateau ; aymeric.rateau@gmail.com
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

    def __init__( self ):
        """Override this to declare a list of input fields if required."""
        # a list of ImportField objects to display
        self.fields = []

    def getPreview( self, params ):
        """Get data to show in a text box to show a preview.
        params is a ImportPluginParams object.
        Returns (text, okaytoimport)
        """

        info = mdfinfo( params.filename )
        f = ''
        f += 'Time: ' + info.HDBlock['Date'] + ' '
        f += info.HDBlock['Time'] + '\n'
        f += 'Author: ' + info.HDBlock['Author'] + '\n'
        f += 'Organisation: ' + info.HDBlock['Organization' ] + '\n'
        f += 'Project Name: ' + info.HDBlock['ProjectName'] + '\n'
        f += 'Vehicle: ' + info.HDBlock['Vehicle'] + '\n' + 'Channel List:\n'
        info.channelNameList.sort()
        for channelName in info.channelNameList:
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
    description = "Reads MDF files from INCA"

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
        for channelName in data.keys():
            if len( data[channelName]['data'] ) > 0 and data[channelName]['data'].dtype != '|S1':
                print( data[channelName]['data'].dtype )
                List.append( ImportDataset1D( channelName, data[channelName]['data'] ) )
        return List

importpluginregistry.append( MdfImportPlugin() )

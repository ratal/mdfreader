# -*- coding: utf-8 -*-
"""Veusz import plugin for MDF files (mdfreader >= 4.0, Veusz >= 1.16).

Installation instructions:
  1. Go to Edit / Preferences / Plugins and add this file.
  2. Ensure mdfreader is installed: pip install mdfreader
"""

import veusz.plugins as plugins
from veusz.plugins.datasetplugin import Dataset1D
from veusz.plugins.field import FieldFloat

try:
    from mdfreader import Mdf, MdfInfo
except ImportError:
    raise ImportError('mdfreader is required. Install with: pip install mdfreader')


class MdfImportPlugin(plugins.ImportPlugin):
    """Import plugin for MDF files (ETAS INCA / CANape)."""

    name = 'MDF Import'
    author = 'Aymeric Rateau'
    description = 'Reads MDF (.mdf / .mf4 / .dat) files'
    promote_tab = 'MDF'
    file_extensions = {'.dat', '.mf4', '.mdf'}

    def __init__(self):
        self.fields = [
            FieldFloat('mult', descr='Resampling period (s)', default=0.1),
        ]

    def getPreview(self, params):
        """Return a text preview of the MDF file metadata."""
        info = MdfInfo(params.filename)

        lines = []
        hd = info.get('HDBlock', {})

        if info.mdfversion < 400:
            lines.append(f"Date/Time : {hd.get('Date', '')} {hd.get('Time', '')}")
            lines.append(f"Author    : {hd.get('Author', '')}")
            lines.append(f"Org       : {hd.get('Organization', '')}")
            lines.append(f"Project   : {hd.get('ProjectName', '')}")
            lines.append(f"Subject   : {hd.get('Subject', '')}")
        else:
            from time import gmtime, strftime
            ns = hd.get('hd_start_time_ns', 0)
            t = gmtime(ns / 1e9)
            lines.append(f"Date/Time : {strftime('%Y-%m-%d %H:%M:%S', t)}")
            comment = hd.get('Comment', {})
            for key, label in [('author', 'Author'), ('department', 'Org'),
                                ('project', 'Project'), ('subject', 'Subject')]:
                if key in comment:
                    lines.append(f"{label:9s}: {comment[key]}")

        lines.append('\nChannels:')
        for ch in info.list_channels():
            lines.append(f'  {ch}')

        return '\n'.join(lines), True

    def doImport(self, params):
        """Load and resample MDF file; return list of Dataset1D objects."""
        data = Mdf(params.filename)
        data.resample(sampling_time=params.field_results['mult'])

        datasets = []
        for ch in data.keys():
            arr = data.get_channel_data(ch)
            if arr is not None and len(arr) > 0 and arr.dtype.kind not in ('S', 'U', 'V'):
                datasets.append(Dataset1D(ch, arr))
        return datasets


plugins.importpluginregistry.append(MdfImportPlugin)  # class, not instance

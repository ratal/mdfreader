# ----------------------------------------------------------------------
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 3 as
# published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see http://www.gnu.org/licenses.
#
# ----------------------------------------------------------------------

__author__ = 'Aymeric Rateau (aymeric.rateau@gmail.com)'
__copyright__ = 'Copyright (c) 2017 Aymeric Rateau'
__license__ = 'GPLV3'
__version__ = "2.7.2"

from sys import path
from os.path import dirname, abspath
root = dirname(abspath(__file__))
path.append(root)
# if it's run as a script or imported within python, this happens
if __name__ == 'mdfreader':
    try:
        from mdfreader.mdfreader import mdf, mdfinfo
    except ImportError:  # python 2-3 differences, not understood
        from mdfreader import mdf, mdfinfo
    from mdf import mdf_skeleton
    from mdf3reader import mdf3
    from mdf4reader import mdf4
    from mdfinfo3 import info3
    from mdfinfo4 import info4, ATBlock  # ,MDFBlock
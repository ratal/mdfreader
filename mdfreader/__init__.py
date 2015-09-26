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
__copyright__ = 'Copyright (c) 2015 Aymeric Rateau'
__license__ = 'GPLV3'
__version__ = "0.2.0"

#if it's run as a script or imported within python, this happens
if __name__ == 'mdfreader':
    from .mdf import mdf_skeleton
    from .mdf3reader import mdf3
    from .mdf4reader import mdf4
    from .mdfinfo3 import info3
    from .mdfinfo4 import info4, MDFBlock, ATBlock
    from .mdfreader import mdf,mdfinfo

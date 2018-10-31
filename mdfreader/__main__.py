from __future__ import print_function
from warnings import warn
import sys
import mdfreader
warn(sys.argv[1:][0])
mdfreader.Mdf(sys.argv[1:][0])

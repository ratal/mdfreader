from __future__ import print_function

import sys
print(sys.argv[1:][0], file=sys.stderr)
import mdfreader
mdfreader.mdf(sys.argv[1:][0])

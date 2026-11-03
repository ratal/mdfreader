import sys
print(sys.argv[1:][0])
import mdfreader
from mdfreader import MdfInfo
from mdfreader.mdf4reader import Mdf4
import cProfile
import pstats
yopinfo = MdfInfo(sys.argv[1:][0])
yop = mdfreader
cProfile.run('yop.Mdf(sys.argv[1:][0])', 'profileStatsFile.prof')
p = pstats.Stats('profileStatsFile.prof')
p.sort_stats('cumulative', 'time').print_stats(400)
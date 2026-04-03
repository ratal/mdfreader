import sys
import mdfreader

if len(sys.argv) < 2:
    print('Usage: python -m mdfreader <mdf_file>', file=sys.stderr)
    sys.exit(1)

mdfreader.Mdf(sys.argv[1])

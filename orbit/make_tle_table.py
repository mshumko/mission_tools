# Make TLE table
import numpy as np
import glob
import argparse
import os

import read_tle

tleDir = '/home/mike/research/firebird/tle'

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("sc_id", help="Spacecraft id", type=str)
parser.add_argument("--tleDir", help='TLE directory, currently set to '
    '{}'.format(tleDir), default=tleDir)
args = parser.parse_args()

# Load in the TLE blacklist here...


# Now read in all of the TLEs
tlePaths = sorted(glob.glob(os.path.join(args.tleDir, '*')))
for path in tlePaths:
    try:
        l1, l2 = read_tle.read_tle(path, args.sc_id)
    except NameError as err:
        print(err)
        continue
    # Now figure out the epoch
    epoch = l1.split()[3] # Look at TLE documentation on the format.
    print(epoch, l1, l2)
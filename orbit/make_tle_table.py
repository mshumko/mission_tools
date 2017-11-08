# Make TLE table
#import numpy as np
import glob
import argparse
import os
import sys
from datetime import datetime, timedelta

import read_tle
sys.path.insert(0, '../misc/')
import dates_in_filenames

tleDir = '/home/mike/research/firebird/tle'

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("sc_id", help="Spacecraft id", type=str)
parser.add_argument("--tleDir", help='TLE directory, currently set to '
    '{}'.format(tleDir), default=tleDir)
args = parser.parse_args()

# Load in the TLE blacklist here...
blackListDates = []
with open('{}_tle_blacklist.txt'.format(args.sc_id.upper()), 'r') as f:
    for line in f:
        if line[0] != '#':
            blackListDates.append(line.strip('\n'))

# Open file for writing the tle table
with open('{}_tle_table.txt'.format(args.sc_id.upper()), 'w') as f:
    f.write('# {} TLE list generated on {}\n# epoch, line1, line2 \n'.format(
        args.sc_id, datetime.now()))
    # Now read in all of the TLEs, and save to list
    tlePaths = sorted(glob.glob(os.path.join(args.tleDir, '*')))
    for path in tlePaths:
        # Get date from current TLE path, and check if has been blacklisted
        basename = os.path.basename(path)
        date = os.path.splitext(basename)[0].split('_')[2:5]
        if '-'.join(date) in blackListDates: # Ingnore if file is blacklisted.
            continue
        try:
            l1, l2 = read_tle.read_tle(path, args.sc_id)
        except NameError as err:
            print(err)
            continue
        # Now figure out the epoch
        rawEpoch = l1.split()[3] # Look at TLE documentation on the format.
        y = 2000 + int(rawEpoch[0:2])
        doy = float(rawEpoch[2:])
        epoch = datetime(y, 1, 1) + timedelta(days=doy-1)
        f.write('{}, {}, {}\n'.format(epoch, l1, l2))
# This script reads in a tle dump file from space-track.org and 
# splits the file into daily files which are saved to saveDir.
import csv
from datetime import datetime, timedelta
import os
import numpy as np

# FOR AC6
#saveDir = '/home/mike/Desktop/ac6_tles'
#header = 'AEROCUBE 6A' # Shows up before first line as an identifier.
## Path to space-track TLE dump file
#dumpPath = '/home/mike/research/ac6/tle/ac6_tle_dump.txt' 

# FOR ELFIN
# saveDir = '/home/mike/research/elfin/tle'
# header = 'ELFIN B' # Shows up before first line as an identifier.
## Path to space-track TLE dump file
# fname = 'elfinb'
# dumpPath = '/home/mike/research/elfin/tle/{}_tle_dump.txt'.format(fname)

# FOR DSX
saveDir = '/home/mike/research/dsx/tle'
header = 'DSX' # Shows up before first line as an identifier.
# Path to space-track TLE dump file
dumpPath = '/home/mike/Desktop/dsx_tle_dump.txt'

with open(dumpPath) as f:
    r = csv.reader(f)
    r = list(r)
    for l1, l2 in zip(r[0::2], r[1::2]): # Loop over each line
        # Calculate the epoch
        rawEpoch = l1[0].split()[3] # Look at TLE documentation on the format.
        y = 2000 + int(rawEpoch[0:2])
        doy = float(rawEpoch[2:])
        epoch = datetime(y, 1, 1) + timedelta(days=doy-1)
        
        # Write each TLE to file.
        with open(os.path.join(saveDir, '{}_tle_{}.txt'.format(
                    header.lower(), epoch.date())), 'a') as s:
            s.write(header + '\n' + l1[0] + '\n' + l2[0] + '\n')
# This script reads in a tle dump file from space-track.org and 
# splits the file into daily files which are saved to saveDir.
import csv
from datetime import datetime, timedelta
import os
import numpy as np

#saveDir = '/home/mike/Desktop/ac6_tles'
#sc_id = 'AEROCUBE 6A' # Shows up before first line as an identifier.
#Path tp space-track TLE dump file
#dumpPath = '/home/mike/research/ac6/tle/ac6_tle_dump.txt' 

saveDir = '/home/mike/research/elfin/tles'
sc_id = 'ELFIN B' # Shows up before first line as an identifier.
#Path tp space-track TLE dump file
dumpPath = '/home/mike/research/elfin/tles/elfin_b_tle_dump.txt' 

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
        with open(os.path.join(saveDir, 'elfinb_tle_{}.txt'.format(epoch.date())), 'a') as s:
            s.write(sc_id + '\n' + l1[0] + '\n' + l2[0] + '\n')

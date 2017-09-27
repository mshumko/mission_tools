from datetime import datetime
import dateutil.parser 
import numpy as np
import sys, os, string
sys.path.insert(0, '/home/mike/Dropbox/0_firebird_research/usefull_functions')
import dates_in_filenames

tleDir = '/home/mike/Dropbox/0_firebird_research/orbit_propagator/TLE/'

def match_tle(tleDir, t):
    """
    NAME: match_tle(tleDir, t)
    USE:  Finds a TLE file in directory tleDir that matches the
          date in t. t can be a datetime object or string.
    RETURNS: full path to the matched TLE.
    MOD:   2017-04-16
    """
    # Convert t to a date if not already.
    if isinstance(t, datetime):
        t = t.date().isoformat()
    else:
        t = dateutil.parser.parse(t).date().isoformat()

    # Now load in the TLE file names and match the input date
    # with the corresponding TLE.
    paths, dates = dates_in_filenames.find_dates(tleDir,
    dateTimeConvert = True)
    dates = list(map(lambda x: x.date().isoformat(), np.array(dates)[:, 0]))
    matchedFileInd = np.where(t == np.array(dates))[0][0]
    # Return the full path to the matched TLE.
    return paths[matchedFileInd]
    

def read_tle(tleFPath, sc_id):
    """
    NAME: read_tle(tleFPath, sc_id)
    USE: Reads in a TLE file from the theFPath input, and 
         looks for two line elements matching the sc_id.    
    RETURNS: Two line element strings.
    MOD:   2017-04-16
    """
    line1, line2 = None, None
    scanInd = 0

    # Now read in the file
    with open(tleFPath, 'r') as tle:
        # Now read the TLE file and find the matching spacecraft.
        for line in tle:
            if scanInd == 1:
                scanInd += 1
                line1 = line
            elif scanInd == 2:
                line2 = line
                return line1, line2

            # If found a matching spacecraft, prepare to
            # save the next two lines. 
            if not line[0].isdigit() and sc_id in line:
                scanInd += 1
    return -9999

tleFDir = match_tle(tleDir, '2017-04-16')
print(read_tle(tleFDir, 'FU3'))

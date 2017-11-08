from datetime import datetime
import dateutil.parser 
import numpy as np
import sys, os, string
sys.path.insert(0, '../misc/')
import dates_in_filenames

tleDir = '/home/mike/research/firebird/tle'

def get_tle(tleDir, date, sc_id):
    # Wrapper to get the TLE.
    tlePath = find_tle(tleDir, date)
    return read_tle(tlePath, sc_id)


def find_tle(tleDir, t):
    """
    NAME: find_tle(tleDir, t)
    USE:  Finds a TLE file in directory tleDir that matches the
          date in t. t can be a datetime object or string.
    RETURNS: full path to the matched TLE.
    MOD:   2017-11-08
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
    try:
        matchedFileInd = np.where(t == np.array(dates))[0][0]
    except IndexError as err:
        print('*********** No file found on {} **************'.format(t))
        raise
    # Return the full path to the matched TLE.
    return paths[matchedFileInd]
    

def read_tle(tleFPath, sc_id):
    """
    NAME: read_tle(tleFPath, sc_id)
    USE: Reads in a TLE file from the theFPath input, and 
         looks for two line elements with a matching sc_id
         in the line before.    
    RETURNS: A tuple of two line elements.
    MOD:   2017-11-08
    """
    # Now read in the file
    with open(tleFPath, 'r') as tle:
        for line in tle:
            if sc_id.upper() in line: # If sc_id is found in line. 
                line1 = next(tle) # Look up iterrator next() command for file I/O.
                line2 = next(tle)
                return line1.strip('\n'), line2.strip('\n')
        raise NameError('Spacecraft {} not found in {}.'.format(sc_id, tleFPath))

if __name__ == '__main__':
    print(get_tle(tleDir, '2017-11-08', 'FU3'))
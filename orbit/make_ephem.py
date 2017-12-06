from datetime import datetime, timedelta
import numpy as np
import dateutil.parser
import glob
import csv
import os

from sgp4.earth_gravity import wgs72 # SGP4 specific libraries
from sgp4.io import twoline2rv # https://pypi.python.org/pypi/sgp4/

from eci_conversion import eci2gdz # For coordinate transformations

class Make_ephem:
    def __init__(self, sc_id, tBounds, dt):
        """
        NAME: Make_ephem(self, sc_id, tBounds, dt)
        USE: This is the front end of a set of classes that 
             uses Two Line Elements (TLE) to propagate orbits.
        ARGS:
            REQUIRED:
                sc_id: Spacecraft id
                tBounds: Time bounds (a 2 element array of datetimes)
                dt: Output ephemeris time step size in seconds.
            OPTIONAL:
                None
        RETURNS: Nothing (Class returns file at the end.)
        EXAMPLE:
            startTime = datetime.now()
            sc_id = 'FU4'
            tRange = [datetime(2017, 6, 28), datetime(2017, 7, 31)]
            dT = 60
            ephem = Make_ephem(sc_id, tRange, dT)
            ephem.loadTleTable()
            ephem.propagateOrbit()
            ephem.saveEphem()
            print('Run time: {}'.format(datetime.now() - startTime))
        AUTHOR: Mykhaylo Shumko
        MOD:   2017-11-15
        """
        self.sc_id = sc_id
        self.tBounds = tBounds
        self.dt = dt
        return

    def loadTleTable(self, tlePath=None):
        """
        NAME:    loadTleTable(self, tlePath=None)
        USE:     This function loads the TLE table generated by
                 the Make_TLE_table class.
        INPUT:   REQUIRED:
                    None
                 OPTIONAL:
                    tlePath=None: The path where TLE table is 
                        located. Looks for file with the name
                        <sc_id>_tle_table.txt
        AUTHOR:  Mykhaylo Shumko
        RETURNS: Nothing (instance of self.tleTable).
                 Also self.startTimes and self.endTimes, arrays of
                 datetimes for each TLE's start and end time. 
                 (Last TLE's end epoch is datetime.max)
        MOD:     2017-11-15
        """
        if tlePath is None:
            tlePath = '{}_tle_table.txt'.format(self.sc_id.upper())
        # Read contents of the file
        self.tleData = np.nan*np.ones((0, 3), dtype=object)
        # Read in the file, and ignore the header.
        with open(tlePath, 'r') as f:
            for line in f:
                if line[0] == '#':
                    continue
                splitLine = np.array(line.strip('\n').split(', '))

                # Create self.tleData which contains the epoch and TLE.
                self.tleData = np.row_stack(
                    (self.tleData[:, 0:3], splitLine))
        self._getTleEpochBounds() # Get time TLE bounds
        self._findTLEidx() # Get valid TLE indicies
        return

    def propagateOrbit(self):
        """
        NAME:    propagateOrbit(self)
        USE:     Propagates orbits with the sgp4 library. It intellegently
                 switches between TLEs at their epoch times.
        INPUT:   REQUIRED:
                    None
                 OPTIONAL:
                    None                    
        AUTHOR:  Mykhaylo Shumko
        RETURNS: Nothing (self.time, self.alt, self.lat, self.lon instances 
                 of arrays)
        MOD:     2017-11-15
        """
        # Create time array and empty position arrays
        self._dataSkeletons() 
        # Loop over valid TLE's and generate ephemeris
        for i in range(self.startIdx, self.endIdx+1):
            # Load TLE
            sat = twoline2rv(self.tleData[i, 1], 
                self.tleData[i, 2], wgs72)
            # Find time indicies to loop over
            idt = np.where((self.times >= self.startTimes[i])
                 & (self.times <= self.endTimes[i]))[0]

            for t in idt: # Loop over times.
                x, v = sat.propagate(self.times[t].year, 
                    self.times[t].month, self.times[t].day,
                    self.times[t].hour, self.times[t].minute,
                    self.times[t].second)
                self._checkError(sat) # Check if propagation was sucessfull
                XGDZ = eci2gdz(x, self.times[t]) # ECI -> DGZ transform
                self.lat[t] = XGDZ['Lat']
                self.lon[t] = XGDZ['Lon']
                self.alt[t] = XGDZ['Alt']
                self.v[t] = np.linalg.norm(v) # Save magntiude of velocity
        # Now filter the data to user times.
        validIdt = np.where((self.times >= self.tBounds[0]) & 
            (self.times <= self.tBounds[1]))[0]
        self.times = self.times[validIdt]
        self.lon = self.lon[validIdt]
        self.lat = self.lat[validIdt]
        self.alt = self.alt[validIdt]
        self.v = self.v[validIdt]
        return

    def saveEphem(self, fPath=None):
        """
        NAME:    saveEphem(self, fPath=None)
        USE:     Saves the (self.time, self.alt, self.lat, self.lon)
                 instances to a csv file given by fPath
        INPUT:   REQUIRED:
                    None
                 OPTIONAL:
                    fPath=None: If None, will save to the current 
                        directory with the name 
                        <sc_id>_<start_date>_<end_date>_LLA_ephemeris.txt
                        If not none, will save to the full path specified.
        AUTHOR:  Mykhaylo Shumko
        RETURNS: Nothing (instance of self.tleTable)
        MOD:     2017-11-15
        """
        if fPath is None:
            fPath = '{}_{}_{}_LLA_ephemeris.csv'.format(self.sc_id.upper(),
                self.tBounds[0].date(), self.tBounds[1].date())
        with open(fPath, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Time (ISO), Lat (deg), Lon (deg), '
                'Alt (km), Vel (km/s)']) # Header
            for row in zip(self.times, self.lat, self.lon, self.alt, self.v):
                writer.writerow(row)

    def _getTleEpochBounds(self):
        """
        Convert the TLE table epoch from string to datetime objects.
        Also find each TLE's end time by indexing from [1:] and filling
        in the end with a large dummy year (datetime.max).
        """
        self.startTimes = np.array(list(map(
            dateutil.parser.parse, self.tleData[:, 0])))
        self.endTimes = np.append(self.startTimes[1:], datetime.max)

    def _findTLEidx(self):
        """
        Find the index of start and end TLEs contain tRange time stamps.
        """
        self.startIdx = np.where((self.tBounds[0] > self.startTimes) & 
            (self.tBounds[0] < self.endTimes))[0]
        self.endIdx = np.where((self.tBounds[1] > self.startTimes) & 
            (self.tBounds[1] < self.endTimes))[0]
        assert len(self.startIdx) and len(self.endIdx), 'No TLEs found for the times specified!'
        self.startIdx = self.startIdx[0]
        self.endIdx = self.endIdx[0]
        return

    def _dataSkeletons(self):
        """
        This function creates empty arrays of correct size for propagation.
        """
        # Figure out when to start time.
        if self.tBounds[0] > self.startTimes[0]:
            startTime = self.tBounds[0]
        else:
            startTime = self.startTimes[self.startIdx].replace(second=0,
                minute=(self.startTimes[self.startIdx].minute+1)%60, 
                microsecond=0)
        if self.tBounds[1] > self.startTimes[-1]:
            # If the user wanted an end time after the epoch of last TLE. (This avoids propagating to year 9999)
            endTime = self.tBounds[1]
        else:
            endTime = self.endTimes[self.endIdx].replace(second=0,
                minute=(self.endTimes[self.endIdx].minute+1)%60,
                microsecond=0)
        
        totSteps = int((endTime - startTime).total_seconds()/self.dt)
        self.times = np.array([startTime + timedelta(seconds=self.dt*i)
            for i in range(totSteps)])

        self.lat = np.nan*np.ones(totSteps, dtype=float)
        self.lon = np.nan*np.ones(totSteps, dtype=float)
        self.alt = np.nan*np.ones(totSteps, dtype=float)
        self.v = np.nan*np.ones(totSteps, dtype=float)
        return

    def _checkError(self, satellite):
        """
        Checks if the SGP4 algorithm sucessfully ran.
        """
        failedInx = np.where(satellite.error != 0)[0]
        if len(failedInx) > 0:
            raise ValueError('SGP4 failed!')    


class Make_TLE_table:
    def __init__(self, sc_id, tleDir='/home/mike/research/firebird/tle', blacklistFname=None):
        """
        NAME: Make_TLE_table(self, sc_id)
        USE:  Makes a csv table of the TLEs for sc_id 
              spacecraft.
        ARGS:
            REQUIRED:
                sc_id: Spacecraft id
            OPTIONAL:
                tleDir='/home/mike/research/firebird/tle': The 
                    directory containing the spacecraft 
                    TLEs. If propagating for FIREBIRD, download 
                    them off of Europa.
                blacklistFname=None: The blacklist filename. Will
                    try to load <<sc_id>>_tle_blacklist.txt file 
                    if none. If no such file exists, one will be
                    created.
        RETURNS: Nothing (Class returns file at the end.)
        EXAMPLE:
            # A file FU3_tle_blacklist.txt will be created
            # and FU3_tle_table.txt
            obj = Make_TLE_table('FU3')
            obj.createTable()
        AUTHOR: Mykhaylo Shumko
        MOD:   2017-11-15
        """
        self.sc_id = sc_id
        self.tleDir = tleDir
        self._loadTLEblacklist(blacklistFname) # Load TLE blacklist file.
        return

    def createTable(self, fname=None):
        """
        NAME: createTable(self, fname=None)
        USE:  Makes a csv table of the TLEs for sc_id 
              spacecraft.
        ARGS:
            REQUIRED:
                none
            OPTIONAL:
                fname=None: The output TLE table filename.
                    If None, will create one with format
                    <<sc_id>>_tle_table.txt
        RETURNS: Nothing (Class returns file at the end.)
        AUTHOR: Mykhaylo Shumko
        MOD:   2017-11-15
        """
        if fname is None:
            fname = '{}_tle_table.txt'.format(self.sc_id.upper())

        # Open file for writing the tle table
        with open(fname, 'w') as f:
            print('Created TLE table file', fname)
            f.write('# {} TLE list generated on {}\n# epoch, line1, line2 \n'.format(
                self.sc_id.upper(), datetime.now())) # Header
            # Now read in all of the TLEs, and save to list
            tlePaths = sorted(glob.glob(os.path.join(self.tleDir, '*')))
            for path in tlePaths:
                # Get date from current TLE path, and check if has been blacklisted
                basename = os.path.basename(path)
                date = os.path.splitext(basename)[0].split('_')[2:5]
                if '-'.join(date) in self.blackListDates: # Ignore if file is blacklisted.
                    continue
                try:
                    l1, l2 = read_tle(path, self.sc_id)
                except NameError as err:
                    print(err)
                    continue
                # Now figure out the epoch
                rawEpoch = l1.split()[3] # Look at TLE documentation on the format.
                y = 2000 + int(rawEpoch[0:2])
                doy = float(rawEpoch[2:])
                epoch = datetime(y, 1, 1) + timedelta(days=doy-1)
                f.write('{}, {}, {}\n'.format(epoch, l1, l2))
        return

    def _loadTLEblacklist(self, fname=None):
        """
        Load in the TLE blacklist file...
        """
        if fname is None:
            fname = '{}_tle_blacklist.txt'.format(self.sc_id.upper())

        # If blacklist file has not been created
        if not os.path.exists(fname):
            with open(fname, 'w') as f:
                f.write("# This is the black list of TLE's"
                    " not to process for {}.\n".format(self.sc_id))
                f.write('# Insert dates for TLEs you want ignored '
                    'in YYYY-MM-DD format.')
                print('Created blacklist file named:', fname)
            
        self.blackListDates = []
        with open(fname, 'r') as f:
            for line in f:
                if line[0] != '#':
                    self.blackListDates.append(line.strip('\n'))
        return


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
        raise NameError('Spacecraft {} not found in {}.'.format(sc_id.upper(), tleFPath))

if __name__ == '__main__':
    for sc_id in ['FU4']:
        tableObj = Make_TLE_table(sc_id)
        tableObj.createTable()
        tBounds = [datetime(2015, 3, 28), datetime.now()]
        dT = 300
        ephemObj = Make_ephem(sc_id, tBounds, dT)
        ephemObj.loadTleTable()
        ephemObj.propagateOrbit()
        ephemObj.saveEphem()

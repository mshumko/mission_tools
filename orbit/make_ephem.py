# This program uses the file output from make_tle_table to make FIREBIRD ephemeris
# for a specified spacecraft, time range, and step size
from datetime import datetime, timedelta
import numpy as np
import dateutil.parser
#import multiprocessing

# SGP4 Specific libraries
# https://pypi.python.org/pypi/sgp4/
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv

from eci_conversion import eci2gdz

class make_ephem:
    def __init__(self, sc_id, tBounds, dt):
        """

        """
        self.sc_id = sc_id
        self.tBounds = tBounds
        self.dt = dt
        # Make the time array
        #ds = (tRange[1] - tRange[0]).total_seconds()
        #self.ephemTime = np.array([tRange[0] + timedelta(seconds=s*dt) 
        #    for s in range(int(ds/dt))])
        return

    def loadTleTable(self, tlePath=None):
        if tlePath is None:
            tlePath = '{}_tle_table.txt'.format(self.sc_id.upper())
        # Read contents of the file
        self.tleData = np.nan*np.ones((0, 3), dtype=object)

        with open(tlePath, 'r') as f:
            for line in f:
                if line[0] == '#':
                    continue
                splitLine = np.array(line.strip('\n').split(','))

                # Create self.tleData which contains the epoch and TLE.
                self.tleData = np.row_stack(
                    (self.tleData[:, 0:3], splitLine))
        self._getTleEpochBounds() # Get time TLE bounds
        self._findTLEidx() # Get valid TLE indicies
        return

    def propagateOrbit(self):
        # Loop over valid TLE's and generate ephemeris
        return

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
        assert len(self.startIdx) > 0 and len(self.endIdx), 'No TLEs found for the times specified!'
        return

                

if __name__ == '__main__':
    startTime = datetime.now()
    sc_id = 'fu3'
    tRange = [datetime(2016, 5, 19), datetime(2016, 7, 1)]
    dT = 60
    ephem = make_ephem(sc_id, tRange, dT)
    ephem.loadTleTable()
    print('Run time: {}'.format(datetime.now() - startTime))

    # for i in zip(ephem.startTimes, ephem.endTimes):
    #     print(i)
    
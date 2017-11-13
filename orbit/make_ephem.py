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
    def __init__(self, sc_id, tRange, dt):
        """

        """
        self.sc_id = sc_id
        # Make the time array
        ds = (tRange[1] - tRange[0]).total_seconds()
        self.ephemTime = np.array([tRange[0] + timedelta(seconds=s*dt) 
            for s in range(int(ds/dt))])
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
                # Get TLE bounds
                self._getTleEpochBounds()
        return

    def propagateOrbit(self)

    def _getTleEpochBounds(self):
        # Convert the TLE table epoch from string to datetime objects
        self.startTimes = np.array(list(map(
            dateutil.parser.parse, self.tleData[:, 0])))
        self.endTimes = np.append(self.startTimes[1:], datetime.max)

                

if __name__ == '__main__':
    startTime = datetime.now()
    sc_id = 'fu3'
    tRange = [datetime(2015, 1, 29), datetime(2015, 2, 23)]
    dT = 60
    ephem = make_ephem(sc_id, tRange, dT)
    ephem.loadTleTable()
    print('Run time: {}'.format(datetime.now() - startTime))

    for i in zip(ephem.startTimes, ephem.endTimes):
        print(i)
    
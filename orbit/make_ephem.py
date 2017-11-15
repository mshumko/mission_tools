# This program uses the file output from make_tle_table to make FIREBIRD ephemeris
# for a specified spacecraft, time range, and step size
from datetime import datetime, timedelta
import numpy as np
import dateutil.parser
import csv

import matplotlib.pyplot as plt # REMOVE WHEN DONE TESTING
import os # REMOVE WHEN DONE TESTING

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
                splitLine = np.array(line.strip('\n').split(', '))

                # Create self.tleData which contains the epoch and TLE.
                self.tleData = np.row_stack(
                    (self.tleData[:, 0:3], splitLine))
        self._getTleEpochBounds() # Get time TLE bounds
        self._findTLEidx() # Get valid TLE indicies
        return

    def propagateOrbit(self):
        """

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
                x, _ = sat.propagate(self.times[t].year, 
                    self.times[t].month, self.times[t].day,
                    self.times[t].hour, self.times[t].minute,
                    self.times[t].second)
                self._checkError(sat) # Check if propagation was sucessfull
                XGDZ = eci2gdz(x, self.times[t]) # ECI -> DGZ transform
                self.lat[t] = XGDZ['Lat']
                self.lon[t] = XGDZ['Lon']
                self.alt[t] = XGDZ['Alt']

        # Now filter the data to user times.
        validIdt = np.where((self.times >= self.tBounds[0]) & 
            (self.times <= self.tBounds[1]))[0]
        self.times = self.times[validIdt]
        self.lon = self.lon[validIdt]
        self.lat = self.lat[validIdt]
        self.alt = self.alt[validIdt]
        return

    def saveEphem(self, fPath=None):
        # Save the ephemeris to csv file.
        if fPath is None:
            fPath = '{}_{}_{}_LLA_ephemeris.txt'.format(self.sc_id.upper(),
                self.tBounds[0].date(), self.tBounds[1].date())
        with open(fPath, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Time (ISO), Lat (deg), Lon (deg), Alt (km)']) # Header
            for row in zip(self.times, self.lat, self.lon, self.alt):
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
        assert len(self.startIdx) > 0 and len(self.endIdx), 'No TLEs found for the times specified!'
        self.startIdx = self.startIdx[0]
        self.endIdx = self.endIdx[0]
        return

    def _dataSkeletons(self):
        startTime = self.startTimes[self.startIdx].replace(second=0,
            minute=self.startTimes[self.startIdx].minute+1, 
            microsecond=0)
        endTime = self.endTimes[self.endIdx].replace(second=0,
            minute=self.endTimes[self.endIdx].minute+1,
            microsecond=0)
        totSteps = int((endTime - startTime).total_seconds()/self.dt)
        self.times = np.array([startTime + timedelta(seconds=self.dt*i)
            for i in range(totSteps)])

        self.lat = np.nan*np.ones(totSteps, dtype=float)
        self.lon = np.nan*np.ones(totSteps, dtype=float)
        self.alt = np.nan*np.ones(totSteps, dtype=float)
        return

    def _checkError(self, satellite):
        # Checks if the SGP4 algorithm sucessfully ran.
        failedInx = np.where(satellite.error != 0)[0]
        if len(failedInx) > 0:
            raise ValueError('SGP4 failed!')


if __name__ == '__main__':
    # Propagate ephemeris
    startTime = datetime.now()
    sc_id = 'fu3'
    tRange = [datetime(2016, 5, 19), datetime(2016, 7, 1)]
    dT = 60
    ephem = make_ephem(sc_id, tRange, dT)
    ephem.loadTleTable()
    ephem.propagateOrbit()
    ephem.saveEphem()
    print('Run time: {}'.format(datetime.now() - startTime))

    # # Load existing ephemeris
    # fDir = '/home/mike/research/firebird/Datafiles/FU_{}/ephem'.format(
    #     sc_id[-1])
    # fname = 'FU3_LLA_camp08_2016-05-19_2016-07-01.csv'
    # N = int((datetime(2016, 7, 1) - datetime(2016, 5, 19)).total_seconds()/60)+1
    # import csv
    # stkAlt = np.ones(N, dtype=float)
    # stkLat = np.ones(N, dtype=float)
    # stkLon = np.ones(N, dtype=float)
    # stkTime = np.ones(N, dtype=object)

    # with open(os.path.join(fDir, fname)) as f:
    #     reader = csv.DictReader(f, 
    #         fieldnames=['Time', 'Lat', 'Lon', 'Alt'])
    #     reader.__next__()
    #     for idx, line in enumerate(reader):
    #         stkTime[idx] = dateutil.parser.parse(line['Time'])
    #         stkAlt[idx] = line['Alt']
    #         stkLat[idx] = line['Lat']
    #         stkLon[idx] = line['Lon']

    # # ### PLOT VALIDATION ###
    # fig, ax = plt.subplots(3, sharex=True)

    # ax[0].plot(stkTime, stkAlt, label='STK alt')
    # ax[1].plot(stkTime, stkLat, label='STK lat')
    # ax[2].plot(stkTime, stkLon, label='STK lon')

    # ax[0].plot(ephem.times, ephem.alt, label='Mikes Alt')
    # ax[1].plot(ephem.times, ephem.lat, label='Mikes Lat')
    # ax[2].plot(ephem.times, ephem.lon, label='Mikes Lon')

    # ax[0].set_ylabel('Alt (km)')
    # ax[0].set_title('STK vs Mikes make_ephem validation')
    # ax[1].set_ylabel('Lat')
    # ax[2].set_ylabel('Lon')
    # ax[0].legend()
    # plt.show()
    
# This file validates the make_ephem.py program
import matplotlib.pyplot as plt
from datetime import datetime
import dateutil.parser
import numpy as np
import csv
import os

import make_ephem

# Propagate ephemeris
startTime = datetime.now()
sc_id = 'FU4'
tRange = [datetime(2016, 5, 19), datetime(2016, 7, 1)]
dT = 60
ephem = make_ephem.Make_ephem(sc_id, tRange, dT)
ephem.loadTleTable()
ephem.propagateOrbit()
ephem.saveEphem()
print('Run time: {}'.format(datetime.now() - startTime))

# Load existing ephemeris
fDir = '/home/mike/research/firebird/Datafiles/FU_{}/ephem'.format(
    sc_id[-1])
fname = 'FU{}_LLA_camp08_2016-05-19_2016-07-01.csv'.format(sc_id[-1])
N = int((datetime(2016, 7, 1) - datetime(2016, 5, 19)).total_seconds()/60)+1
stkAlt = np.ones(N, dtype=float)
stkLat = np.ones(N, dtype=float)
stkLon = np.ones(N, dtype=float)
stkTime = np.ones(N, dtype=object)

with open(os.path.join(fDir, fname)) as f:
    reader = csv.DictReader(f, 
        fieldnames=['Time', 'Lat', 'Lon', 'Alt'])
    reader.__next__()
    for idx, line in enumerate(reader):
        stkTime[idx] = dateutil.parser.parse(line['Time'])
        stkAlt[idx] = line['Alt']
        stkLat[idx] = line['Lat']
        stkLon[idx] = line['Lon']

# ### PLOT VALIDATION ###
fig, ax = plt.subplots(3, sharex=True)

ax[0].plot(stkTime, stkAlt, label='STK alt')
ax[1].plot(stkTime, stkLat, label='STK lat')
ax[2].plot(stkTime, stkLon, label='STK lon')

ax[0].plot(ephem.times, ephem.alt, label='Mikes Alt')
ax[1].plot(ephem.times, ephem.lat, label='Mikes Lat')
ax[2].plot(ephem.times, ephem.lon, label='Mikes Lon')

ax[0].set_ylabel('Alt (km)')
ax[0].set_title('STK vs Mikes make_ephem validation')
ax[1].set_ylabel('Lat')
ax[2].set_ylabel('Lon')
ax[0].legend()

plt.show()
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 16:26:41 2017

@author: mike
"""
import spacepy.datamodel as dm
import spacepy.time as spt
import os
import dateutil.parser
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import datetime
import generate_ephem

import sys
sys.path.insert(0, '/home/mike/FIREBIRD/data_processing/magnetic_ephemeris')

from generate_magnetic_ephemeris import loadEphemeris

# Generate the ephemeris
sc_id = 3
    
if sc_id == 3:
    line1 = ('1 40377U 15003B   17065.77165762  .00002287  00000-0  10232-3 0  9990')
    line2 = ('2 40377  99.1141 229.3954 0137748 341.9237  17.7117 15.14671553115505')
elif sc_id == 4:
    line1 = ('1 40378U 15003C   17066.16165030  .00002541  00000-0  11295-3 0  9996')
    line2 = ('2 40378  99.1143 229.8923 0137929 340.5500  19.0493 15.14695458115566')
      
times = [datetime.datetime(2017, 3, 6, 0, 0, 1) + datetime.timedelta(minutes = i) \
for i in range(1440)]

p = generate_ephem.SGP4_ephemeris(line1, line2, times)
X, errCode = p.propagate()
failedInx = p.checkError()

if len(failedInx) == 0:
    p.write_to_file('/home/mike/Dropbox/0_firebird_research/' + 
    'orbit_propagator/validation', 'FU' + str(sc_id) + '_test.txt')

# generate_ephem validation
fDir = '/home/mike/Dropbox/0_firebird_research/orbit_propagator/validation'
STKname = 'FU3_STK_LLA_2017-03-06_to_2017-07-06_gen_w_2017-03-06T18UT_TLE.csv'
myFname = 'FU3_test.txt'

if 'STKdata' not in globals():
    STKdata = loadEphemeris(os.path.join(fDir, STKname))   
    STKtimes = [dateutil.parser.parse(i[0].decode()) for i in STKdata]
    
mydata = dm.readJSONheadedASCII(os.path.join(fDir, myFname))
myTimes = list(map(dateutil.parser.parse, mydata['dateTime']))
#myTimes = spt.Ticktock(mydata['dateTime']).UTC
startInd = np.where(np.array(STKtimes) == myTimes[0])[0][0]
endInd = np.where(np.array(STKtimes) == myTimes[-1])[0][0]
    
fig = plt.figure(figsize=(15, 10), dpi=80, facecolor = 'grey')
gs = gridspec.GridSpec(3, 1)
altPlt = fig.add_subplot(gs[0, 0])
latPlt = fig.add_subplot(gs[1, 0], sharex = altPlt)
lonPlt = fig.add_subplot(gs[2, 0], sharex = altPlt)

altPlt.plot(myTimes, mydata['Alt'], 'r', label = 'My implemenation')
altPlt.plot(STKtimes[startInd:endInd], np.array(STKdata['Alt_km'])[startInd:endInd], 'b', label = 'STK implemenation')
altPlt.set_ylabel('Altitude (km)')
altPlt.legend()

latPlt.plot(myTimes, mydata['Lat'], 'r', label = 'My implemenation')
latPlt.plot(STKtimes[startInd:endInd], np.array(STKdata['Lat_deg'])[startInd:endInd], 'b', label = 'STK implemenation')
latPlt.set_ylabel('Latitude (deg)')
#latPlt.legend()

lonPlt.plot(myTimes, mydata['Lon'], 'r', label = 'My implemenation')
lonPlt.plot(STKtimes[startInd:endInd], np.array(STKdata['Lon_deg'])[startInd:endInd], 'b', label = 'STK implemenation')
lonPlt.set_ylabel('Longitude (deg)')
#lonPlt.legend()
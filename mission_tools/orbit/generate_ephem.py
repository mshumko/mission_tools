#SGP4 library located at
#https://pypi.python.org/pypi/sgp4/

# SGP4 Specific libraries
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv

# General libraries
import spacepy.datamodel as dm
from datetime import datetime, timedelta
import dateutil.parser
import os, sys
import numpy as np
import time

# My libraries.
sys.path.insert(0, '/home/mike/Dropbox/0_firebird_research/usefull_functions')
import dates_in_filenames
from eci_conversion import eci2gdz

class SGP4_ephemeris:
    def __init__(self, t, **kwargs):
        """
        NAME: SGP4_ephemeris
        USE:   
        RETURNS:
        MOD:   2017-04-16
        """
        self.line1 = kwargs.get('line1', None)
        self.line2 = kwargs.get('line2', None)
        self.sc_id = kwargs.get('sc_id', None)
        self.tleDir = kwargs.get('tleDir', '/home/mike/' + 
        'Dropbox/0_firebird_research/orbit_propagator/TLE/')
        self.sphere_shape = kwargs.get('sphere_shape', True)

        # If did not supply datetime objects, convert!
        if type(t[0]) is not datetime:
            self.times = [dateutil.parser.parse(i) for i in t]
        else:
            self.times = t

		# Execute this comment block if the user wants the
		# program to load in the TLE from the date of the first
		# datetime in the t input.
        if (self.line1 is None) and (self.line2 is None):
            if sc_id is not None:
                # Read TLE, "intelegently" and save to line1
                # and line2 class instances.
                self.getTwoLineElements(self.sc_id, t[0])
            else:
                raise ValueError("Either specify the " +
                "spacecraft id, ie 'FU3', 'FU4', or supply" + 
                " the two-line elements!")
        self.sat = twoline2rv(self.line1, self.line2, wgs72)
        return

    def getTwoLineElements(self, sc_id, t):
        fPath = self._match_tle(self.tleDir, t)
        self.line1, self.line2 = self._read_tle(fPath, sc_id)
        return self.line1, self.line2
   
    def propagate(self):
        self.X = {}
        self.X['Alt'], self.X['Lat'], self.X['Lon'], self.errCode = \
        (-9999*np.ones_like(self.times) for i in range(4))
        
        for t in range(len(self.times)):
            # Propagate
            position, velocity = \
            self.sat.propagate(self.times[t].year, 
            self.times[t].month, self.times[t].day, 
            self.times[t].hour, self.times[t].minute,
            self.times[t].second)
            
            # Convert TEME coordinates to Lat Lon Alt.
            XGDZ = eci2gdz(position, self.times[t],
            sphere_shape = self.sphere_shape)
            self.X['Alt'][t] = XGDZ['Alt']
            self.X['Lat'][t] = XGDZ['Lat']
            self.X['Lon'][t] = XGDZ['Lon']
            
            self.errCode[t] = self.sat.error
        return self.X, self.errCode
        
    def checkError(self):
        failedInx = np.where(self.errCode != 0)[0]
        
        if len(failedInx) > 0:
            print('Error! SGP4 failed for ', len(failedInx), 
            ' points')
        return failedInx
        
    def write_to_file(self, fDir, fName):
        file = dm.SpaceData()
        if self.sphere_shape:
            file.attrs['Global'] = ('Sphere conversion gravity' 
            + ' model')
        else:
            file.attrs['Global'] = ('Oblate spheroid ' + 
            'conversion gravity model')
            
        order = ['dateTime', 'Lat', 'Lon', 'Alt']
        dateTimeAttrs = {'FORMAT':"ISO", "UNIT":'UTC', \
        'TITLE':'ISO 8601 standard date time'}
        altAttrs = {'TITLE':'Altitude', 'UNIT':'km'}
        latAttrs = {'TITLE':'Latitude', 'UNIT':'degrees'}
        lonAttrs = {'TITLE':'Longitude', 'UNIT':'degrees'}
        timeStrings = [i.isoformat() for i in self.times]
        
        file['dateTime'] = dm.dmarray(timeStrings, attrs=dateTimeAttrs)
        file['Alt'] = dm.dmarray(self.X['Alt'], attrs=altAttrs)
        file['Lat'] = dm.dmarray(self.X['Lat'], attrs=latAttrs)
        file['Lon'] = dm.dmarray(self.X['Lon'], attrs=lonAttrs)
        dm.toJSONheadedASCII(os.path.join(fDir, fName), file, \
        depend0='dateTime', order=order)
        pass

    def _match_tle(self, tleDir, t):
        """
        NAME: match_tle(tleDir, t)
        USE:  Finds a TLE file in directory tleDir that matches
              the date in t. t can be a datetime object or
              string.
        RETURNS: full path to the matched TLE.
        MOD:   2017-04-16
        """
        # Convert t to a date if not already.
        if isinstance(t, datetime):
            t = t.date().isoformat()
        else:
            t = dateutil.parser.parse(t).date().isoformat()

        # Now load in the TLE file names and match the input
        # date with the corresponding TLE.
        paths, dates = dates_in_filenames.find_dates(tleDir,
        dateTimeConvert = True)
        dates = list(map(lambda x: x.date().isoformat(),
        np.array(dates)[:, 0]))
        try:
            matchedFileInd = np.where(t == np.array(dates)
            )[0][0]
        except IndexError:
            print('Error, no matching times found!')
            sys.exit()

        # Return the full path to the matched TLE.
        return paths[matchedFileInd]

    def _read_tle(self, tleFPath, sc_id):
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


if __name__ == '__main__':
    start_time = time.time()
    sc_id = 4

    if sc_id == 3:
        line1 = ('1 40377U 15003B   17316.87243881 +.00003378 +00000-0 +14390-3 0  9990')
        line2 = ('2 40377 099.1116 170.5446 0132965 229.1728 129.7910 15.16264918153522')
    elif sc_id == 4:
        line1 = ('1 40378U 15003C   17317.12687510  .00003231  00000-0  13770-3 0  9990')
        line2 = ('2 40378  99.1117 170.8967 0133183 228.2288 130.7496 15.16281094153579')
              
    days = 40
    times = [datetime(2017, 11, 13, 0, 0, 0) + timedelta(
    minutes=i) for i in range(1440*days)]
    
#    times = [datetime(2017, 11, 8, 0, 0, 0) + timedelta(
#    seconds = i) for i in range(86401)]
        
    p = SGP4_ephemeris(times, line1 = line1, line2 = line2)
    X, errCode = p.propagate()
    failedInx = p.checkError()
    
    if len(failedInx) == 0:
        p.write_to_file('/home/mike/research/' 
        'mission-tools/orbit/data/', 
        'FU' + str(sc_id) + '_SGP4_LLA_' + 
        times[0].date().isoformat() + '_to_' + 
        times[-1].date().isoformat() + 
        '_gen_w_2017-11-13_TLE.txt')
    else:
        print('Some values failed to propegate!')
    print("Ephemeris execution time %s seconds" % (time.time() - start_time))

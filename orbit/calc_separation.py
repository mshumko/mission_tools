# This script makes a file that calculates the cross-spacecraft separation and 
# in-track lag.
import dateutil.parser 
import numpy as np
import csv

Re = 6371 # km.

class CalcSep:
    def __init__(self, fPathA, fPathB, convertTimes=False):
        """
        Convert times will convert timestamps to datetime objects (slow!)
        """
        self.convertTimes = convertTimes
        self.dataA = self._readFile(fPathA)
        self.dataB = self._readFile(fPathB)
        return

    def calcDistance(self):
        """ 
        This function will blindly calculate the spacecraft separation 
        assuming the time stamps all line up.
        """
        # Combine the two data sets
        lat = np.stack((self.dataA['lat'], self.dataB['lat']))
        lon = np.stack((self.dataA['lon'], self.dataB['lon']))
        # Calculate the average altitude for both spacecraft
        alt = np.mean(np.stack((self.dataA['lat'], self.dataB['lat'])), axis=0)
        self.crowFliesDist(lat, lon, alt)
        return

    def crowFliesDist(self, lat, lon, alt):
        """ 
        This function uses the haverside foruma to calculate the great-circle
        distance between two points.
        
        lat/lon have to be dimentions of 2 or 2*nT. alt dimentions should be 1 or nT
        """
        # Convert degrees to radians
        lat = np.deg2rad(lat)
        lon = np.deg2rad(lon)
        dlat = np.subtract(lat[0], lat[1])
        dlon = np.subtract(lon[0], lon[1])
        # Calculate distance
        a = np.sin(dlat/2)**2 + np.cos(lat[0])*np.cos(lat[1])*np.sin(dlon/2)**2
        c = 2*np.arctan2(np.sqrt(a), np.sqrt(1-a))
        self.d = (Re+alt)*c
        return self.d

    def saveData(self, savePath):
        """
        This function will save data in a CSV format with time and distance 
        columns.
        """
        with open(savePath, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['dateTime', 'separation (km)'])
            for line in zip(self.dataA['dateTime'], self.d):
                writer.writerow(line)

    def _readFile(self, fPath):
        """ This function uses the csv module to read in the ephem file """
        with open(fPath, 'r') as f:
            reader = csv.reader(f)
            next(reader) # Skip header
            d = np.array(list(reader))

        data = {}
        if self.convertTimes:
            data['dateTime'] = np.array([dateutil.parser.parse(i[0]) 
                for i in d]) 
        else:
            data['dateTime'] = np.array([i[0] for i in d])
        data['lat'] = np.array([i[1] for i in d], dtype=float)
        data['lon'] = np.array([i[2] for i in d], dtype=float)
        data['alt'] = np.array([i[3] for i in d], dtype=float)
        data['vel'] = np.array([i[4] for i in d], dtype=float)
        return data
    
# Testing
if __name__ == '__main__':
    fPathA = ('/home/mike/research/mission-tools/orbit/data/'
        'FU3_2015-03-28_2017-12-06_LLA_ephemeris.csv')
    fPathB = ('/home/mike/research/mission-tools/orbit/data/'
        'FU4_2015-03-28_2017-12-06_LLA_ephemeris.csv')
    savePath = '2015_2017_FB_separation.txt'
    dObj = CalcSep(fPathA, fPathB, convertTimes=True)
    dObj.calcDistance()
    #dObj.saveData(savePath)

    import matplotlib.pyplot as plt
    plt.plot(dObj.dataA['dateTime'], dObj.dataA['alt'])
    plt.ylabel('FU3 Altitude (km)')
    plt.xlabel('Time')
    plt.show()    
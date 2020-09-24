import numpy as np
import os
import spacepy.datamodel as dm
import datetime

import downloadKp

class downloadAE:
    """
    Class that downloads the AE index for a specific year. Possible kwargs 
    are: url, outputDir, and outputFName. The year is a required argument. 
    If RUN_SAVE_SCRIPT kwargs is True, the data will be downloaded and saved
    to whatever directory and filename specified. 
    """
    def __init__(self, year, **kwargs):
        self.url = kwargs.get('url', \
                'ftp://ftp.ngdc.noaa.gov/STP/GEOMAGNETIC_DATA/INDICES/' + \
                'AURORAL_ELECTROJET/ONE_MINUTE/')
        self.outputDir = kwargs.get('outputDir', \
                '/home/mike/Desktop/')
        self.outputFName = kwargs.get('outputFName', \
                str(year) + '_minute_ae.txt')
        self.iType = kwargs.get('iType', 'AE')
        self.year = year
        assert self.year > 10, 'Error, input year must be in the full year '+\
        'format (ex 2003).'
        
        if kwargs.get('RUN_SAVE_SCRIPT', False):
            self.downloadData()
            self.saveToFile()
        
    def downloadData(self):
        """
        NAME: downloadData(self)
        USE: Given a year, figures out which subdirectories to pull from in
        the URL provided in __init__. 
        RETURNS: dateTime and index arrays which contain one of the A index 
            types
        MOD: 2017-02-14
        AUTHOR: Mykhaylo Shumko
        """ 
        
        assert self.iType in ['AE', 'AL', 'AO', 'AU'], 'Error, not correct index'+\
        ' type! Check "iType" keyword'
        
        urlDir = self.url + 'ae_' + str(self.year) + '_minute.txt'
        raw_data = downloadKp.getFTPdata(urlDir, max_tries = 5)
        
        # Find the indicies where there is the correct index type.
        dTypes = np.array([raw_data[i][21:23].decode()\
        for i in range(len(raw_data))])
        idx = np.where(dTypes == self.iType)[0]
        
        self.indexArr = -9999*np.ones(60*len(idx)) 
        self.dateTimeArr = 60*len(idx)*[None]
        
        for i in range(len(idx)):
            # Now decode and parse the data for the correct index type.
            line = raw_data[idx[i]].decode()
            if int(line[12:14]) < 75:
                year = int(str(20) + line[12:14])
            else:
                year = int(str(19) + line[12:14])
                
            month = int(line[14:16])
            day = int(line[16:18])
            hour = int(line[19:21])
            
            # Read and convert the 1 minute AE values to integers.
            indicies = list(map(int, line[34:394].split())) 
            # Now fill minute arrays
            self.indexArr[(i)*60:(i+1)*60] = indicies
            self.dateTimeArr[i*60:(i+1)*60] = \
            [datetime.datetime(year, month, day, hour) + \
            k*datetime.timedelta(minutes = 1) for k in range(60)]

        return  self.indexArr, self.dateTimeArr
        
    def saveToFile(self):
        """
        Saves the dateTime and index data to a JSON headed ASCII file in the 
        directory and filename specified by __init__().
        """
        saveObj = dm.SpaceData()
        outputFName = os.path.join(self.outputDir, self.outputFName)
                  
        dateTimeAttrs = {'FORMAT':"ISO", "UNIT":'UTC', \
        'DESCRIPTION':'ISO 8601 standard date time'}
        IndexAttrs = {'UNIT':'None?', 'DESCRIPTION':self.iType + ' geomagnetic index.'}
        saveObj['dateTime'] = dm.dmarray(self.dateTimeArr, attrs = dateTimeAttrs)
        saveObj[self.iType] = dm.dmarray(self.indexArr, attrs = IndexAttrs)
        order = ['dateTime', self.iType]
        dm.toJSONheadedASCII(outputFName, saveObj, order = order, depend0='dateTime') 
        
if __name__ == '__main__':
    years =  np.arange(1990, 2003)
    for year in years:
        d = downloadAE(year, RUN_SAVE_SCRIPT = True, \
        outputDir = '/home/mike/research/indicies/AE')

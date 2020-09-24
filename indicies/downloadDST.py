# DST downloader class

import numpy as np
#import sys
import os
import spacepy.datamodel as dm
from datetime import date
import downloadKp

class downloadDST:
    """
    Class that downloads the DST index for a specific year. Possible kwargs 
    are: url, outputDir, and outputFName. The year is a required argument. 
    
    Example: This creates a file called 2015_quicklook_DST.txt on mike's 
    Desktop directory. 
    >>>data = downloadDST(2015)
    >>>data.downloadData()
    >>>data.processDST()
    >>>data.writeToJSONfile()
    """
    def __init__(self, year, **kwargs):
        self.url = kwargs.get('url', \
                'http://wdc.kugi.kyoto-u.ac.jp/dst_realtime/')
        self.outputDir = kwargs.get('outputDir', \
                '/home/mike/Desktop/')
        self.outputFName = kwargs.get('outputFName', \
                str(year) + '_quicklook_DST.txt')
        #self.times = np.array([])
        #self.DST = np.array([])
        self.year = year
        
    def loadRawDST(self, url):
        """
        NAME: loadRawDST(self, url)
        USE: Downloads and line splits raw html from url provided.
        RETURNS: Raw HTML with DST
        MOD: 2016-10-31
        AUTHOR: Mykhaylo Shumko
        """
        self.raw_data = downloadKp.getFTPdata(url)
        
    def downloadData(self):
        """
        NAME: downloadData(self)
        USE: Given a year, figures out which subdirectories to pull from in
        the URL mprovided in __init__. Takes the raw data from 
        downloadKp.getFTPdata and appends it to self.rawLineData with the 
        first three columns as year, month, day, and the other 24 columns are
        the hour values of the DST.
        RETURNS: Parsed lines of int DST values from the raw HTML
        (self.rawLineData)
        MOD: 2016-10-31
        AUTHOR: Mykhaylo Shumko
        """ 
        self.rawLineData = []
        dateTimeToday = date.today()
        if self.year == dateTimeToday.year:
            monthsNum = dateTimeToday.month
        else:
            monthsNum = 12
        for month in range(1, monthsNum+1):
            urlDir = self.url + str(self.year) + "%02d" % month + '/index.html'
            raw_data = downloadKp.getFTPdata(urlDir)
            
            i = 40
            while i < len(raw_data):
                try:
                    if (len(raw_data[i]) != 0 and int(raw_data[i][0:2]) in range(1, 32)):
                        tempLine = [int(numeric_string) for numeric_string in raw_data[i].split()]
                        tempLine.insert(0, month)
                        tempLine.insert(0, self.year)
                        self.rawLineData.append(tempLine)
                    i+=1
                except ValueError:
                    i+=1
            i = 0
        return
        
    def processDST(self):
        """
        NAME: processDST(self)
        USE: Processes the line by line paresed DST index, and creates a
        corresponding time and DST arrays at a one hour cadence. The error flag
        is caught, and not included into the arrays.
        RETURNS: Dictionary with 'dateTime':self.times, 'DST':self.DST
        MOD: 2016-11-01
        AUTHOR: Mykhaylo Shumko
        """ 
        self.timeswErrors = 24*len(self.rawLineData)*[99999999999999999999999999999999]
        #self.timeswErrors = np.empty(24*len(self.rawLineData), dtype = str)
        self.DSTwErrors = 99999999999999999999999999999999*np.ones(24*len(self.rawLineData))
        
        for md in range(len(self.rawLineData)):
            for hr in range(24):
                if max(self.rawLineData[md]) < 9999:
                    self.timeswErrors[md*24+hr] = '{}-{}-{}T{}:00:00'.format(
                        self.rawLineData[md][0], "%02d" % self.rawLineData[md][1], 
                        "%02d" % self.rawLineData[md][2], "%02d" % hr)
                    self.DSTwErrors[24*md+hr] = self.rawLineData[md][3+hr]      
        self.timeswErrors = np.array(self.timeswErrors)
        self.times = self.timeswErrors[np.where(self.DSTwErrors != 1e+32)[0]]
        self.DST = self.DSTwErrors[np.where(self.DSTwErrors != 1e+32)[0]]
                    
        return {'dateTime':self.times, 'DST':self.DST}
        
    def writeToJSONfile(self):
        """
        NAME: writeToJSONfile(self)
        USE: Given the outputDir and outputFName kwargs, will save the DST and
        corresponsing times into a JSON headed ASCII. If kwargs not provided,
        will save in same directory that exists on Mike's computer. 
        RETURNS: JSON headed ASCII file with dateTimes and DST.
        MOD: 2016-11-01
        AUTHOR: Mykhaylo Shumko
        """
        self.output = dm.SpaceData()
        outputFName = os.path.join(self.outputDir, self.outputFName)
                      
        dateTimeAttrs = {'FORMAT':"ISO", "UNIT":'UTC', 'DESCRIPTION':'ISO 8601 standard date time'}
        DSTAttrs = {'UNIT':'none', 'DESCRIPTION':'Disturbance Storm Time index. ' + 
            'The Dst index is an index of magnetic activity derived from a '+
            'network of near-equatorial geomagnetic observatories that measures' + 
            ' the intensity of the globally symmetrical equatorial electrojet ' + 
            '(the "ring current"). Dst is maintained at NGDC and is available'+
            ' via FTP from 1957 to the present.'}
        self.output['dateTime'] = dm.dmarray(self.times, attrs = dateTimeAttrs)
        self.output['DST'] = dm.dmarray(self.DST, attrs = DSTAttrs)
        order = ['dateTime', 'DST']
        dm.toJSONheadedASCII(outputFName, self.output, order = order, depend0='dateTime') 
            
#        
if __name__ == '__main__':
    data = downloadDST(2017)
    data.downloadData()
    data.processDST()
    data.writeToJSONfile()

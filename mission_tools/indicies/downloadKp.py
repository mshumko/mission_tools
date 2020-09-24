# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 14:00:04 2016

@author: mike
"""

import numpy as np
import sys, os
import spacepy.datamodel as dm
from datetime import date

def downloadKp(year, url = ''):
    """
    NAME: downloadKpIndex(year, url = '')
    USE:  Downloads the Kp index for the supplied year argument from an optional 
    URL that defaults to ftp://ftp.ngdc.noaa.gov/STP/GEOMAGNETIC_DATA/INDICES/KP_AP/ + year
    
    
    **** Kp index formatting ****    
    The National Geophysical Data Center provides the observations here with no re-
    strictions on their use.  Please contact us at the address below with your com-
    ments and questions about the form and the content of this information product
    or about the measurements themselves.
    
    NATIONAL GEOPHYSICAL DATA CENTER
    325 Broadway     Mail Code E/GC2
    Boulder, Colorado 80303-3328 USA
    
    Telephone: (303) 497-6135   Telex: 592811 NOAA MASC BDR
    FAX:       (303) 497-6513
    
    
    -------------------------------------------------------------------------------
    SELECTED GEOMAGNETIC AND SOLAR ACTIVITY INDICES
    -------------------------------------------------------------------------------
    
    
    WORLDWIDE INDICES--------------------------------------------------------------
    The subscript "p" means planetary and designates a global magnetic activity in-
    dex.  The following 13 observatories, which lie between 46 and 63 degrees north
    and south geomagnetic latitude, now contribute to the planetary indices: Ler-
    wick (UK), Eskdalemuir (UK), Hartland (UK), Ottawa (Canada), Fredericksburg
    (USA), Meannook (Canada), Sitka (USA), Eyrewell (New Zealand), Canberra (Aus-
    tralia), Lovo (Sweden), Brorfelde (Denmark), Wingst (Germany), and Niemegk
    (Germany).
    
    
    THREE-HOUR-RANGE INDEX K-------------------------------------------------------
    K indices isolate solar particle effects on the earth's magnetic field; over a
    3-hour period, they classify into disturbance levels the range of variation of
    the more unsettled horizontal field component.  Each activity level relates al-
    most logarithmically to its corresponding disturbance amplitude.  Three-hour
    indices discriminate conservatively between true magnetic field perturbations
    and the quiet-day variations produced by ionospheric currents.
    
    K indices range in 28 steps from 0 (quiet) to 9 (greatly disturbed) with frac-
    tional parts expressed in thirds of a unit.  A K-value equal to 27, for exam-
    ple, means 2 and 2/3 or 3-; a K-value equal to 30 means 3 and 0/3 or 3 exactly;
    and a K-value equal to 33 means 3 and 1/3 or 3+.  The arithmetic mean of the K
    values scaled at the 13 observatories listed above gives Kp.
    ------------------------------------------------------------------------------
    RETURNS: ISO formatted UTC and Kp for each 3 hour interval. 
    AUTHOR: Mykhaylo Shumko
    MOD:     2016-09-27
    """  
    
    if url == '':
        url = 'ftp://ftp.ngdc.noaa.gov/STP/GEOMAGNETIC_DATA/INDICES/KP_AP/' + str(year)
    
    
    # Sometimes the NOAA website gets overwheled, so try to connect to the url a few times.
    data = getFTPdata(url)
            
    outputTimes = [None] * (8*len(data))
    outputKp = np.zeros(8*len(data))
    
    # Date parsing assumes that the year is two numbers. 
    # Iterate over each line (each day)
    for dayData in range(len(data)):
        if str(year)[0:2] == '20':
            year = 2000 + int(data[dayData][0:2])
        elif str(year)[0:2] == '19':
            year = 1900 + int(data[dayData][0:2])
        else:
            print('Year out of range')
            return -1
        month = int(data[dayData][2:4])
        day = int(data[dayData][4:6])
        
        # Iterate over each 3-hour block, starting at 0000 - 0300 UT, and ending at 2100 - 2400 UT.
        for threeHourBlock in range(8):
            outputTimes[dayData*8 + threeHourBlock] = str("%04d" % year) + '-' + \
            str("%02d" % month) + '-' + str("%02d" % day) + 'T' + "%02d" % (3*threeHourBlock) + ':00:00'
            outputKp[dayData*8 + threeHourBlock] = int(data[dayData][12 + 2*threeHourBlock : 14 + 2*threeHourBlock])
        
    return {"dateTime":outputTimes, 'kp':outputKp}
    
def downloadCurrentKp(year, url = ''):
    if url == '':
        url = 'ftp://ftp.swpc.noaa.gov/pub/indices/old_indices/'
        
    # Determine if we can just download the yearly data, or have to download it
    # in quarters.
    year = int(year)
    dateTimeToday = date.today()
    if year != dateTimeToday.year:
        url = url + str(year) + '_DGD.txt'
        data = getFTPdata(url)
    else:
        # Now figure out which quarters to download.
        if(dateTimeToday.month <= 3):
            quarter = [str(year) + 'Q1_DGD.txt']
        elif(dateTimeToday.month > 3 and dateTimeToday.month <= 6):
            quarter = [str(year) + 'Q1_DGD.txt', str(year) + 'Q2_DGD.txt']
        elif(dateTimeToday.month > 6 and dateTimeToday.month <= 9):
            quarter = [str(year) + 'Q1_DGD.txt', str(year) + 'Q2_DGD.txt', 
                str(year) + 'Q3_DGD.txt']
        elif(dateTimeToday.month > 9 and dateTimeToday.month <= 12):
            quarter = [str(year) + 'Q1_DGD.txt', str(year) + 'Q2_DGD.txt', 
                str(year) + 'Q3_DGD.txt', str(year) + 'Q4_DGD.txt']
        
        # Now download quarters. 
        quarterData = []
        for i in range(len(quarter)):
            quarterData.append(getFTPdata(url + quarter[i]))
        data = np.hstack((quarterData))    
    
    # If unable to download data.
    if len(data) == 0:
        print('Could not download data from ' + url)
        return {"dateTime":-99, 'kp':-99}
        
    # Decode the data for Python 3!
    if sys.version_info[0] == 3:
        data = decodeData(data)
        
    # Now parse the data! 
    # Calculate array length without the metadata.
    data = np.array(data)[np.where(np.array(data) != '')[0]]
    firstChar =[x[0] for x in data]
    validInd = np.where((np.array(firstChar) != '#') & (np.array(firstChar) != ':'))[0]
        
    outputTimes = [None] * 8*len(validInd)
    outputKp = -1*np.ones(8*len(validInd), dtype = int)

    # Date parsing assumes that the year is two numbers. 
    # Iterate over each line (each day)
    for dayData in range(len(validInd)):

        # Parse the data, by spacing out the minus signes to then parse it.
        # .replace creates a space in front of the minus sign. 
        # .split splits up the string line with the space delimiter.
        # list and map function convers the string numbers into ints.
            
        parsedLine = list(map(int, data[validInd[dayData]].replace('-', ' -').split()))
        year = parsedLine[0]
        month = parsedLine[1]
        day = parsedLine[2]
        
        # Iterate over each 3-hour block, starting at 0000 - 0300 UT, and ending at 2100 - 2400 UT.
        for threeHourBlock in range(8):         
            outputTimes[dayData*8 + threeHourBlock] = str("%04d" % year) + '-' + \
            str("%02d" % month) + '-' + str("%02d" % day) + 'T' + "%02d" % (3*threeHourBlock) + ':00:00'
            outputKp[dayData*8 + threeHourBlock] = 10*parsedLine[-8 + threeHourBlock]
    #return data, validInd
    return {"dateTime":outputTimes, 'kp':outputKp}
        
def getFTPdata(url, max_tries = 5):
    """
    Works with Python 2.7 and 3.5!
    """
    
    tried = 0 
    connected = False

    # Attemp to download data from url
    if sys.version_info[0] == 2:
        import urllib2
        while connected == False:
            try:
                response = urllib2.urlopen(url) 
                connected = True
                html = response.read()
                data = html.splitlines()
            except urllib2.URLError:
                print("Could not connect to: " + url + ". Trying again")
                tried += 1
                if tried == max_tries:
                    #print("Tried to connect to " + url + " 5 times. Aborting")
                    raise IOError("Tried to connect to " + url + " 5 times. Aborting")
                    #return -1

        
    elif sys.version_info[0] == 3:
        import urllib.request
        while connected == False:
            try:
                with urllib.request.urlopen(url) as response:
                    html = response.read()
                connected = True
                data = html.splitlines()
            except urllib.error.URLError:
                print("Could not connect to: " + url + ". Trying again")
                tried += 1
                if tried == max_tries:
                    raise IOError("Tried to connect to " + url + " 5 times. Aborting")
                    #print("Tried to connect to " + url + " 5 times. Aborting")
                    #return -1
    return data    

def decodeData(data):
    """
    Decode data for Python 3.5
    """
    decodedData = [None]*len(data)
    for i in range(len(data)):
        decodedData[i] = data[i].decode()
    return decodedData    
    
def writeToJSON(year, fDir, fName, newDataFlag = 0):
    if newDataFlag == 0:
        kpData = downloadKp(year)
    elif newDataFlag == 1:
        kpData = downloadCurrentKp(year)

    if kpData == -1:
        return

    kp = dm.SpaceData()
    outputFName = os.path.join(fDir, fName)
                  
    dateTimeAttrs = {'FORMAT':"ISO", "UNIT":'UTC', 'DESCRIPTION':'ISO 8601 standard date time'}
    kpAttrs = {'UNIT':'10*Kp', 'DESCRIPTION':'Global geomagnetic storm index'}
    kp['dateTime'] = dm.dmarray(kpData['dateTime'], attrs = dateTimeAttrs)
    kp['kp'] = dm.dmarray(kpData['kp'], attrs = kpAttrs)
    order = ['dateTime', 'kp']
    dm.toJSONheadedASCII(outputFName, kp, order = order, depend0='dateTime') 
    
if __name__ == '__main__':
    # Update Kp inidicies for 2015 and 2016.
    for i in [2018]:
        print(i)
        writeToJSON(i, './indicies/', str(i)+'_kp.txt', newDataFlag=1)

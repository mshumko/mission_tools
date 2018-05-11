import numpy as np
import dateutil.parser
import multiprocessing
from datetime import datetime
import itertools

import csv

# Libs for testing
import glob, os, sys
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
from matplotlib import dates

try:
    sys.path.insert(0, '/home/mike/Dropbox/0_grad_work')
    import directories
except SystemError:
    print('Could not import the directories.py file, '
        'please supply data directories manualy!')

def read_ac_data(filePath, dType=None, verbose=False):
    """
    This function reads in AC6 data products in CSV formatted files. If dType is not
    specified, it will try to find it automatically.
    """
    if dType is None:
        availTypes = ['10Hz', 'survey', 'coords', 'att']
        for s in availTypes:
            if s in os.path.basename(filePath):
                dType = s
    assert dType is not None, 'Could not identify correct data type!'
    
    with open(filePath, 'r') as f:
        reader = csv.reader(f)
        keys = next(reader)
        rawData = np.array(list(reader))
        assert rawData.shape[0] > 1, 'File is empty!'

    data = {}
    for col, key in enumerate(keys):
        data[key] = rawData[:, col].astype(float)

    # Convert decimal part of seconds to microseconds.
    data['musec'] = np.array(list(map(
        lambda x: int(1E6*(x - int(x))), 
        data['second'])))
    # Prepare input for starmap
    inputTuple = zip(data['year'].astype(int), data['month'].astype(int), 
        data['day'].astype(int), data['hour'].astype(int), data['minute'].astype(int), 
        data['second'].astype(int), data['musec']) 
    # get a dateTime array
    data['dateTime'] = np.array(list(itertools.starmap(datetime, inputTuple)))
    return data

def read_ac_data_wrapper(sc_id, date, dType='10Hz', tRange=None, plot=False):
    """
    This is a plotting wrapper for AC6 data. This plots the dos1 and dos2 rates
    in the top panel, lat/lon in the middle, and OPQ MLT, L shell, and particle
    traped flag on the bottom pannel. 
    """
    ### Load in the data ###
    #path = '/home/ms30715/ssd_data/ac6/ac6{}/ascii/level2'.format(sc_id.lower())
    path = directories.ac6_dir(sc_id)
    splitDate = date.date().isoformat().split('-')
    dateJoined = ''.join(splitDate)
    files = glob.glob(os.path.join(path, 
            'AC6-{}_{}_L2_{}_V03.csv'.format(sc_id.upper(), dateJoined, dType)))
    assert len(files) == 1, 'None or > 1 AC6 files found in {}'.format(path)
    fPath = files[0]
    data = read_ac_data(fPath, dType = dType)

    ### Filter all of the data by time ###
    if tRange is not None:
        timeInd = np.where((data['dateTime'] > tRange[0]) &
            (data['dateTime'] < tRange[1]))[0]
        assert len(timeInd) > 0, 'ERROR: No time filetered data found!'
        
        for key in data.keys():
            data[key] = data[key][timeInd]
    
    if plot:
        ### PLOT DATA ###
        fig = plt.figure(figsize = (10, 8), dpi = 80)
        gs = gridspec.GridSpec(3, 1)
        dosPlt = plt.subplot(gs[0, 0])
        posPlt = plt.subplot(gs[1, 0], sharex = dosPlt)
        magPosPlt = plt.subplot(gs[2, 0], sharex = dosPlt)
        
        # Plot the dos rates
        dos1Valid = np.where(data['dos1rate'] != -1E31)[0]
        dos2Valid = np.where(data['dos2rate'] != -1E31)[0]
        dos3Valid = np.where(data['dos3rate'] != -1E31)[0]
        validL = np.where(data['Lm_OPQ'] != -1E31)[0]

        dosPlt.plot(data['dateTime'][dos1Valid], data['dos1rate'][dos1Valid], 
            label = 'dos1rate') # > 35 keV electron channel
        dosPlt.plot(data['dateTime'][dos2Valid], data['dos2rate'][dos2Valid], 
            label = 'dos2rate') # > 35 keV electron channel
        # On AC6A, dos3 responds to > 1 MeV electrons and > 20 MeV protons.
        # On AC6B, dos3 measures mainly > 20 MeV
        dosPlt.plot(data['dateTime'][dos3Valid], data['dos3rate'][dos3Valid],
            label = 'dos3rate') 
        dosPlt.legend(loc = 1)
        dosPlt.set(ylabel = 'Dose', title = os.path.basename(fPath))
        dosPlt.axes.get_xaxis().set_visible(False)
        
        # Plot the geographic position
        posPlt.plot(data['dateTime'], data['lat'], label = 'lat')
        posPlt.plot(data['dateTime'], data['lon'], label = 'lon')
        posPlt.legend(loc = 1)
        posPlt.set(ylabel = 'Degrees', xlabel = 'UTC')
        posPlt.axes.get_xaxis().set_visible(False)
        
        # Plot the magnetic position
        magPosPlt.plot(data['dateTime'], data['flag'], label = 'Flag')
        magPosPlt.plot(data['dateTime'], data['MLT_OPQ'], label = 'MLT (OPQ)')
        magPosPlt.plot(data['dateTime'][validL], data['Lm_OPQ'][validL], 
            label = r'Lm (OPQ)')
        magPosPlt.plot(data['dateTime'][validL], data['Loss_Cone_Type'][validL], 
            label = 'Loss_Cone_Type')
        magPosPlt.set(ylabel = 'Unitless')
        magPosPlt.legend(loc = 1)

        # Tick formatting that does not quite work yet!
        #hfmt = dates.DateFormatter('%H:%M:%S')
        #magPosPlt.xaxis.set_major_locator(dates.MinuteLocator())
        #magPosPlt.xaxis.set_major_formatter(hfmt)
        #plt.setp(magPosPlt.xaxis.get_majorticklabels(), rotation=30)
        dosPlt.axes.get_xaxis().set_visible(False)
        gs.tight_layout(fig)
        plt.show()
    return data

if __name__ == '__main__':
    import time
    startTime = time.time()
    sc_id = 'A'
    date = datetime(2018, 4, 21)
    #tRange = [datetime(2018, 2, 27, 10, 30), datetime(2018, 2, 27, 10, 50)]
    data = read_ac_data_wrapper(sc_id, date, dType = 'survey', tRange=None, 
        plot=True)
    print('Run time {}'.format(time.time() - startTime))

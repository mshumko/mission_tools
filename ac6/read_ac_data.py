import numpy as np
import dateutil.parser
import multiprocessing
from datetime import datetime

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

#ac6path = directories.ac6_dir

def read_ac_data(filePath, dType = None, verbose = False):
    """
    This function reads in AC6 data products in CSV formatted files. If dType is not
    specified, it will try to find it automatically.
    """
    # Checks to catch improper function use.
    # If not specified, scan the filename for a dType.
    if dType is None:
        baseName, extention = os.path.splitext(filePath)
        baseNameSplit = filePath.split('_')
        dType = baseNameSplit[4]
        assert extention in ['.csv', '.txt'], 'Error, not a typical file exention!'
    assert dType in ['coords', 'survey', 'att', 'atts', '10Hz'],  'Specify a valid dType! ({})'.format(dType)
    
    # Use the correct list of keys for each type of data.
    if dType == 'coords':
        dataKeys = ['alt', 'lat', 'lon', 'dos1rate', 'dos2rate', 'dos3rate',
            'flag', 'Lm_IGRF', 'Bmag_IGRF', 'MLT_IGRF', 'InvLat_IGRF', 'Lm_OPQ',
            'Bmag_OPQ', 'MLT_OPQ', 'InvLat_OPQ', 'Loss_Cone_Type', 'Bx_GEO',
            'By_GEO', 'Bz_GEO', 'Beq', 'I', 'K', 'K_Z', 'Lstar', 'Lstar_Z', 'hmin', 'hmin_Z', 
            'Loss_Cone_Near', 'Loss_Cone_Far', 'B100N', 'LAT100N', 'LON100N', 
            'B100S', 'LAT100S', 'LON100S']
    elif dType == 'survey':
       dataKeys = ['alt', 'lat', 'lon', 'X_GEO', 'Y_GEO', 'Z_GEO', 'dos1l', 'dos1m', 
            'dos1rate', 'dos2l', 'dos2m', 'dos2rate', 'dos3l', 'dos3m', 'dos3rate', 'flag', 
            'Sample_Rate', 'Lm_OPQ', 'Bmag_OPQ', 'MLT_OPQ', 'InvLat_OPQ', 
            'Loss_Cone_Type', 'Lstar', 'hmin', 'Alpha', 'Alpha_Eq', 'Beta', 'Phi_B',
            'Dist_In_Track', 'Lag_In_Track', 'Dist_Cross_Track_Horiz', 
            'Dist_Cross_Track_Vert', 'Dist_Total'] 
    elif (dType == 'att') or (dType == 'atts'):
        dataKeys = ['alt', 'lat', 'lon', 'dos1rate', 'dos2rate', 'dos3rate', 'flag', ' Alpha',
            'Alpha_X', 'Alpha_Y', 'Alpha_Eq', 'Beta', 'Beta_X', 'Beta_Y', ' Phi_B', 
            'OmegaX_GEO', 'OmegaY_GEO', 'OmegaZ_GEO', 'B_Spin', 'Spin_Sun']
    elif dType == '10Hz': # Burst mode
        dataKeys = ['alt', 'lat', 'lon', 'dos1l', 'dos1m', 'dos1rate', 'dos2l', 'dos2m',
            'dos2rate', 'dos3l', 'dos3m', 'dos3rate', 'flag', 'Subcom', 'Lm_OPQ', 
            'Bmag_OPQ', 'MLT_OPQ', 'InvLat_OPQ', 'Loss_Cone_Type', 'K_Z', 'Lstar_Z',
            'hmin_Z', 'Alpha', 'Beta', 'Dist_In_Track', 'Lag_In_Track', 
            'Dist_Cross_Track_Horiz', 'Dist_Cross_Track_Vert', 'Dist_Total']
        pass
    if verbose: print('Reading in data of type: {}'.format(dType))
    # Read in the raw data
    rawData = np.genfromtxt(filePath, delimiter = ',', skip_header = 1)
    assert len(rawData.shape) == 2, 'Error, the data is not 2D (Empty file?)'
    #for i in rawData:
    #    print(i[:6])
    # Format time, last term is to format the fractional second to microseconds
    dateTime = np.array([datetime(int(i[0]), int(i[1]), int(i[2]), int(i[3]), 
        int(i[4]), int(i[5]), int((i[5] - int(i[5]))*1e6)) for i in rawData])

#    dateTime = np.array([])
#    for i in rawData:
#        #print(i[:6])
#        try:
#            dateTime = np.append(dateTime, datetime(int(i[0]), int(i[1]), 
#            int(i[2]), int(i[3]), int(i[4]), int(i[5]), int((i[5] - int(i[5]))*1e6)))
#        except ValueError:
#            print(i[:6])
#            raise
        
    # Populate the rest of the arrays
    if verbose: print('Populating dictionary with data keys.')
    data = {}
    data['dateTime'] = dateTime
    dOffset = 6 # Data offset because the date/time is stored in sperate columns
    for i in range(len(dataKeys)):
        data[dataKeys[i]] = rawData[:, i + dOffset]
    return data

def read_ac_data_wrapper(sc_id, date, dType = '10Hz', tRange = None, plot = False):
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
    sc_id = 'A'
    date = datetime(2017, 1, 19)
    #tRange = [datetime(2016, 8, 30, 22, 47), datetime(2016, 8, 30, 22, 50)]
    data = read_ac_data_wrapper(sc_id, date, dType = '10Hz', tRange = None, 
        plot = True)
    #for i in data['dateTime']:
    #    print(i)

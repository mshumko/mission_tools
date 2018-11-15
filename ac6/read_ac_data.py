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
import matplotlib.dates as mdates
from matplotlib.ticker import FuncFormatter, MaxNLocator

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

def read_ac_data_wrapper(sc_id, date, dType='10Hz', tRange=None):
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
    return data

def plot_data(data, sc_id, dtype):
    """
    Plot AC6 data
    """
    # Get plot labels. This is very bad coding practice!
    # Too lazy to rewrite it all.
    # label_ouput =  _plotLabels(data)
    # global time = label_ouput[0]
    # global numTimes = label_ouput[1]
    # global labels = label_ouput[2]

    _, ax = plt.subplots(figsize=(12, 8))

    date = data['dateTime'][0].date()
    idx1 = np.where(data['dos1rate'] != -1E31)[0]
    idx2 = np.where(data['dos2rate'] != -1E31)[0]
    idx3 = np.where(data['dos3rate'] != -1E31)[0]

    # plot dos 1-3 rate channels
    ax.plot(data['dateTime'][idx1], data['dos1rate'][idx1], 
        label = 'dos1rate') # > 35 keV electron channel
    ax.plot(data['dateTime'][idx2], data['dos2rate'][idx2], 
        label = 'dos2rate') # > 35 keV electron channel
    # On AC6A, dos3 responds to > 1 MeV electrons and > 20 MeV protons.
    # On AC6B, dos3 measures mainly > 20 MeV
    ax.plot(data['dateTime'][idx3], data['dos3rate'][idx3],
        label = 'dos3rate') 

    # Format axes
    ax.set_yscale('log')
    ax.set_title('AC6-{} {} {}'.format(sc_id.upper(), dtype, date))
    ax.set_ylabel('dos rate [counts/s]')
    ax.legend()

    # Start interactive session
    ax.format_coord = lambda x, y: '{}, {}'.format(
            dates.num2date(x).replace(tzinfo=None).isoformat(), round(y))
    ax.figure.canvas.mpl_connect('key_press_event', _plotMouseTime)

    # Format x-axis labels
    ax.xaxis.set_major_formatter(FuncFormatter(format_fn))
    ax.set_xlabel('time\nL\nMLT\nlat\nlon')
    ax.xaxis.set_label_coords(-0.1,-0.03)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)
    plt.show()
    return

def _plotMouseTime(event):
    """ 
    # Initialize interactive sesson
    When a user presses 't', this function will print the spacecraft info.
    """
    #print('Not implemented yet!')
    # if event.key == 't':
    #     time = mdates.num2date(event.xdata).replace(tzinfo=None).isoformat()
        # Lt = self.L[np.argmin(np.abs(event.xdata - mdates.date2num(self.time)))]
        # print(time, 'L =', Lt)
    # elif event.key == 'm':
    #     tRange = [t.replace(tzinfo=None) for t in mdates.num2date(ax.get_xlim())]
    #     self.plot_map(tRange)
    return

def _plotLabels(data, skip_n=5):
    ### FORMAT X-AXIS to show more information ###
    data['Lm_OPQ'] = np.round(data['Lm_OPQ'], decimals=1)
    L = np.copy(data['Lm_OPQ']).astype(object)
    L[L < 0] = ''
    time = data['dateTime']
    # This code is a nifty way to format the x-ticks to my liking.
    labels = ['{}\n{}\n{}\n{}\n{}'.format(
                t.replace(microsecond=0).time(),
                L, round(MLT,1), round(lat,1), round(lon,1)) for 
                (t, L, MLT, lat, lon) in zip(
                time[::skip_n], L[::skip_n], data['MLT_OPQ'][::skip_n], 
                data['lat'][::skip_n], data['lon'][::skip_n])]  
    numTimes = mdates.date2num(time[::skip_n])
    return time, numTimes, labels 

def format_fn(tick_val, tick_pos):
    """
    The tick magic happens here. pyplot gives it a tick time, and this function 
    returns the closest label to that time. Read docs for FuncFormatter().
    """
    idx = np.argmin(np.abs(numTimes-tick_val))
    return labels[idx]

if __name__ == '__main__':
    #### If running in interactive mode, use argparse ###
    import argparse

    parser = argparse.ArgumentParser(description=('This script plots the '
        'AC6 level 2 data.'))
    parser.add_argument('sc_id',  type=str,
        help=('This is the AC6 unit arument (A or B)'))
    parser.add_argument('date', nargs=3, type=int,
        help=('This is the data of the data to plot')) 
    parser.add_argument('-d', '--dtype', type=str, default='10Hz',
        help=('AC6 data type to plot (10Hz or survey)'))  
    parser.add_argument('-p', '--plot', type=bool, default=True,
        help=('Plot AC6 data'))
    args = parser.parse_args()

    date = datetime(*args.date)
    data = read_ac_data_wrapper(args.sc_id, date, dType=args.dtype, tRange=None)

    if args.plot:
        time, numTimes, labels = _plotLabels(data)
        plot_data(data, args.sc_id, args.dtype)
   

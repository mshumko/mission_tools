import numpy as np
import dateutil.parser
import multiprocessing
from datetime import datetime
import itertools
import csv
import pandas as pd
from pandas.plotting import register_matplotlib_converters
import matplotlib.colors as colors
import cartopy.crs as ccrs

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

register_matplotlib_converters()

def read_ac_data(filePath, dType=None, verbose=False, use_pandas=True):
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

    if use_pandas:
        data = pd.read_csv(filePath, na_values='-1e+31')
        data['dateTime'] = pd.to_datetime(
            data[['year', 'month', 'day', 'hour', 'minute', 'second']])

        assert data.shape[0] > 1, 'File is empty!'
    else:
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


def get_ac6_path(sc_id, date, dType):
    path = directories.ac6_dir(sc_id)
    splitDate = date.date().isoformat().split('-')
    dateJoined = ''.join(splitDate)
    files = glob.glob(os.path.join(path, 
            'AC6-{}_{}_L2_{}_V03.csv'.format(sc_id.upper(), dateJoined, dType)))
    assert len(files) == 1, 'None or > 1 AC6 files found in {}'.format(path)
    return files[0]

def read_ac_data_wrapper(sc_id, date, dType='10Hz', tRange=None, use_pandas=True):
    """
    This is a plotting wrapper for AC6 data. This plots the dos1 and dos2 rates
    in the top panel, lat/lon in the middle, and OPQ MLT, L shell, and particle
    traped flag on the bottom pannel. 
    """
    ### Load in the data ###
    #path = '/home/ms30715/ssd_data/ac6/ac6{}/ascii/level2'.format(sc_id.lower())
    fPath = get_ac6_path(sc_id, date, dType)
    data = read_ac_data(fPath, dType=dType, use_pandas=use_pandas)

    ### Filter all of the data by time ###
    if tRange is not None:
        timeInd = np.where((data['dateTime'] > tRange[0]) &
            (data['dateTime'] < tRange[1]))[0]
        assert len(timeInd) > 0, 'ERROR: No time filetered data found!'
        
        for key in data.keys():
            data[key] = data[key][timeInd]
    return data

class Plot_AC6:
    def __init__(self, data, sc_id, dtype):
        """
        Plot the AC6 time series data as well as AC6's location.
        """
        self.data = data
        self.time, self.numTimes, self.labels = self._plotLabels(data)
        self.sc_id = sc_id
        self.dtype = dtype
        return

    def plot_data(self):
        """
        Plot AC6 data
        """
        _, self.ax = plt.subplots(figsize=(12, 8))

        date = self.time[0].date()
        # idx1 = np.where(data['dos1rate'] != -1E31)[0]
        # idx2 = np.where(data['dos2rate'] != -1E31)[0]
        # idx3 = np.where(data['dos3rate'] != -1E31)[0]

        # plot dos 1-3 rate channels
        self.ax.plot(data['dateTime'], data['dos1rate'], 
            label = 'dos1rate') # > 35 keV electron channel
        self.ax.plot(data['dateTime'], data['dos2rate'], 
            label = 'dos2rate') # > 35 keV electron channel
        # On AC6A, dos3 responds to > 1 MeV electrons and > 20 MeV protons.
        # On AC6B, dos3 measures mainly > 20 MeV
        self.ax.plot(data['dateTime'], data['dos3rate'],
            label = 'dos3rate') 

        # Format axes
        self.ax.set_yscale('log')
        self.ax.set_title('AC6-{} {} {}'.format(self.sc_id.upper(), self.dtype, date))
        self.ax.set_ylabel('dos rate [counts/s]')
        self.ax.legend(loc=1)

        # Start interactive session
        self.ax.format_coord = lambda x, y: '{}, {}'.format(
                dates.num2date(x).replace(tzinfo=None).isoformat(), round(y))
        self.ax.figure.canvas.mpl_connect('key_press_event', self._plotMouseTime)

        # Format x-axis labels
        self.ax.xaxis.set_major_formatter(FuncFormatter(self.format_fn))
        self.ax.set_xlabel('time\nL\nMLT\nlat\nlon\nflag\nLCT')
        self.ax.xaxis.set_label_coords(-0.1,-0.02)
        plt.tight_layout()
        plt.subplots_adjust(bottom=0.16)
        plt.show()
        return

    def plot_map(self, tRange, channel='dos1rate', decimate=10):
        """
        This function plots the map of the orbit with the same time bounds
        as the time series.
        """
        data_flt = self.data[(self.data['dateTime'] > tRange[0]) & 
                             (self.data['dateTime'] < tRange[1])]
        # idx = np.where((self.data['dateTime'] > tRange[0]) & 
        #             (d['dateTime'] < tRange[1]))[0]
        fig = plt.figure(figsize=(12, 6))
        ax = plt.subplot(111, projection=ccrs.PlateCarree())
        ax.stock_img()
        #ax.coastlines()
        sc = ax.scatter(data_flt.lon[::decimate], data_flt.lat[::decimate], 
                    c=data_flt[channel][::decimate],
                    transform=ccrs.PlateCarree())
        # load and plot L shell data
        lons = np.load('/home/mike/research/mission_tools'
                        '/misc/irbem_l_lons.npy')
        lats = np.load('/home/mike/research/mission_tools'
                        '/misc/irbem_l_lats.npy')
        L = np.load('/home/mike/research/mission_tools'
                        '/misc/irbem_l_l.npy')
        levels = np.arange(2, 10, 2)
        CS = plt.contour(lons, lats, L, levels=levels, colors='k')
        plt.clabel(CS, inline=1, fontsize=10, fmt='%d')
        
        # Mark starting point with a red star
        ax.text(data_flt.lon.iat[0], data_flt.lat.iat[0], '*',
                ha='center', va='center', fontsize=20, color='red',
                transform=ccrs.PlateCarree())

        fig.suptitle('AC6{} ground track on {}\n{} to {}'.format(
            self.sc_id, tRange[0].date(), tRange[0].replace(microsecond=0).time(), 
            tRange[1].replace(microsecond=0).time()), fontsize=16)
        plt.tight_layout()
        # Plot colorbar
        fig.subplots_adjust(right=0.89, top=0.9, bottom=0.1)
        cbar_ax = fig.add_axes([0.9, 0.15, 0.05, 0.7])
        fig.colorbar(sc, orientation='vertical', cax=cbar_ax,
                label=channel)
        plt.show()
        return

    def _plotMouseTime(self, event):
        """ 
        When a user presses 'm', a map will be plotted of the spacecraft location.
        """
        if event.key == 'm':
            tRange = [t.replace(tzinfo=None) for t in mdates.num2date(self.ax.get_xlim())]
            self.plot_map(tRange)
        return

    def _plotLabels(self, data, skip_n=5):
        ### FORMAT X-AXIS to show more information ###
        data['Lm_OPQ'] = np.round(data['Lm_OPQ'], decimals=1)
        L = pd.DataFrame(data['Lm_OPQ'])#.astype(object))
        L = L.replace(np.nan, '', regex=True)
        time = data['dateTime']
        # This code is a nifty way to format the x-ticks to my liking.
        labels = ['{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(
                    t.replace(microsecond=0).time(),
                    L, round(MLT,1), round(lat,1), round(lon,1), flag, lct) for 
                    (t, L, MLT, lat, lon, flag, lct) in zip(
                    time[::skip_n], L.loc[::skip_n, 'Lm_OPQ'], data['MLT_OPQ'][::skip_n], 
                    data['lat'][::skip_n], data['lon'][::skip_n], 
                    data['flag'][::skip_n], data['Loss_Cone_Type'][::skip_n])]  
        numTimes = mdates.date2num(time[::skip_n])
        return time, numTimes, labels 

    def format_fn(self, tick_val, tick_pos):
        """
        The tick magic happens here. pyplot gives it a tick time, and this function 
        returns the closest label to that time. Read docs for FuncFormatter().
        """
        idx = np.argmin(np.abs(self.numTimes-tick_val))
        return self.labels[idx]

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
    import time
    t = time.time()
    data = read_ac_data_wrapper(args.sc_id, date, dType=args.dtype, 
            tRange=None)

    if args.plot:
        p = Plot_AC6(data, args.sc_id, args.dtype)
        p.plot_data()

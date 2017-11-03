# This class will open and analyze GOES electron and proton flux data.
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
import os, sys
import itertools

from netCDF4 import Dataset

try:
    sys.path.insert(0, '/home/mike/Dropbox/0_grad_work')
    import directories
except SystemError:
    print('Could not import the directories.py file, '
        'please supply data directories manualy!')

# Set global font size and color loop to be constistant
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 15}
plt.rc('font', **font)
colorArr = ['r', 'b', 'g', 'k']
colors = itertools.cycle(colorArr)

class GOES_reader:
    def __init__(self, sc_id, date, tRange=None, dType=None):
        """

        """
        self.sc_id = sc_id
        self.date = date
        splitDate = date.date().isoformat().split('-') 
        self.joinedDate = ''.join(splitDate)
        self.dType = dType
        return

    def read_GOES(self):
        """

        """
        fname = 'g{0}_{1}_e1ew_4s_{2}_{2}.nc'.format(self.sc_id, self.dType,
            self.joinedDate)
        filePath = os.path.join(directories.goes_dir(self.sc_id), fname)
        self.goesData = Dataset(filePath)
        self.get_datetimes()
        return

    def get_datetimes(self):
        self.time = [datetime(1970, 1, 1) + timedelta(milliseconds=i) for i 
            in self.goesData['time_tag'][:]] 

    def calc_magLoc(self):
        """

        """
        # Try to import IRBEM
        try:
            from IRBEM import MagFields
        except:
            print("Error: can't find IRBEM!")
            raise

        N = len(self.time)
        # Input variables
        alt = 35790*np.ones(N) # km (geostationary orbit)
        lon = -self.goesData['west_longitude'][0]*np.ones(N)
        lat = np.zeros(N)
        # Output variables
        self.MLT = np.nan*np.ones(N)  
        self.L = np.nan*np.ones(N)

        model = MagFields()
        for i in range(N):
            model.make_lstar({'dateTime':self.time[i], 'x1':alt[i], 
                'x2':lat[i], 'x3':lon[i]}, {'Kp':20})
            self.MLT[i] = model.make_lstar_output['MLT'][0]
            self.L[i] = model.make_lstar_output['Lm'][0]
        #X = {'dateTime': }
    
    def plot_electron_flux(self, ax=None, pltLabels=True, pltLegend=True, formatTimes=True):
        """

        """
        if ax is None:
            fig = plt.figure(figsize=(15, 10), dpi=80, facecolor = 'white')
            gs = gridspec.GridSpec(1,1)
            self.ax = fig.add_subplot(gs[0, 0])
        else:
            self.ax = ax

        if formatTimes:
            myFmt = mdates.DateFormatter('%H:%M:%S')
            self.ax.xaxis.set_major_formatter(myFmt)

        self.ax.plot(self.time, self.goesData['E1W_UNCOR_FLUX'][:], 
            c=next(colors), label='GOES{} {}'.format(self.sc_id, 
                self.goesData['E1W_UNCOR_FLUX'].long_label))
        if pltLabels:
            self.ax.set(xlabel='UTC', ylabel=r'Flux $(e/(cm^2 \ s \ sr)$', 
                title='GOES{} Electron Flux from {}'.format(self.sc_id, 
                    self.date.date().isoformat()))
        if pltLegend:
            self.ax.legend(loc=1)
        if ax is None:
            plt.show()
        return

if __name__ == '__main__':
    fig = plt.figure(figsize=(15, 10), dpi=80, facecolor = 'white')
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0, 0])
    bx = ax.twinx()#fig.add_subplot(gs[1, 0], sharex=ax)

    d = GOES_reader(13, datetime(2017, 3, 31), dType='epead')
    d.read_GOES()
    #print(d.goesData.variables)
    d.calc_magLoc()
    bx.plot(d.time, d.MLT, '--', c=colorArr[0])
    d.plot_electron_flux(pltLegend=False, ax=ax)

    d2 = GOES_reader(15, datetime(2017, 3, 31), dType='epead')
    d2.read_GOES()
    d2.calc_magLoc()
    bx.plot(d2.time, d2.MLT, '--', c=colorArr[1])
    bx.set(ylabel='MLT')
    d2.plot_electron_flux(pltLegend=False, ax=ax)

    ax.legend(loc=1)
    fig.autofmt_xdate()
    gs.tight_layout(fig)
    plt.show()
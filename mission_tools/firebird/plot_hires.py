#!/usr/bin/env python3
"""
This script plots the HiRes data from the FIREBIRD-II mission.
It accepts sc_id (3 or 4) and date (yyyy mm dd) arguments.
"""
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FuncFormatter
import matplotlib
import matplotlib.dates as mdates
import matplotlib.lines as mlines
import matplotlib.colors as colors
import cartopy.crs as ccrs
import argparse
import os 

import spacepy.datamodel as dm
import spacepy.time as spt
import dateutil.parser
from datetime import datetime, timedelta, date
import numpy as np

plt.style.use('bmh')

class plotHiResData:
    def __init__(self, sc_id, date, level=2):
        self.date = date
        self.sc_id = sc_id
        # Load HIRes data
        self.get_hires(self.sc_id, date, level=level)
        self._labels() # Make plotting labels
        return
        
    def plotTimeSeries(self, **kwargs):
        """
        NAME:    plotHiResData(sc_num, level, hiresName)
        USE:     This function plots high resolution (HiRes) data from FIREBIRD-II collimated detectors. 
                 sc_num is the spacraft id, level is the data processing level. NOTE:
                 after running update_hires(), the HiRes data is at level 0. hiresName
                 is the name of the daata file to plot.
        RETURNS: The HiRes plot of log counts vs UTC. 
        MOD:     2016-08-22
        """    
        #level = kwargs.get('level', 2)
        self.dataKey = kwargs.get('dataKey', 'Col_counts')
        errors = kwargs.get('errors', False)
        yscale = kwargs.get('yscale', 'log')
        #cadence = kwargs.get('cadence', 18.75E-3)
        self.ax = kwargs.get('ax', None)
        G = kwargs.get('G', 9)
        c = ['r', 'g', 'b', 'c', 'm','k']
        cadence = 1000*float(self.hires.attrs['CADENCE'])

        if self.ax is None:
            f, self.ax = plt.subplots(figsize=(8, 6))
        
        assert self.dataKey in ['Col_flux', 'Col_counts', 'Sur_counts', 'Sur_flux'], \
        'Epecify the correct data key! (Sur_counts, Sur_flux, Col_flux or Col_counts)'

        # Find breaks in the times
        self.dT = np.array([(self.time[i+1] - 
            self.time[i]).total_seconds() for i in 
            range(len(self.time)-1)])
        self.tBreaks = np.where(np.abs(
            np.convolve(self.dT, [-0.5, 0.5])) > 1)[0]
        self.tBreaks = np.insert(
            self.tBreaks, 0, 0)
        self.tBreaks = np.append(self.tBreaks, 
            len(self.time))

        if errors:
            if 'counts' in self.dataKey.split('_'):
                yerr = np.sqrt(self.hires['Col_counts'])                
            else:
                raise NotImplemented('Flux uncertnaity not implemeted yet!')
                # yerr = [np.sqrt(self.hires['Col_counts'][:, i])/\
                # (cadence*G*self.energyWidths[i]) for i in range(5)]
    
        # Plot the N energy channels.
        if 'flux' in self.dataKey:
            N = 5 # Can't get an integral channel flux.
        else:
            N = 6
        labelHandles = [None]*N

        for i in range(N): 
            label = mlines.Line2D([], [], color=c[i], markersize=15, 
                label=self.energyRanges[int(i)])
            labelHandles[i] = label
        
            # Loop over the continous chunks of times
            for itt in range(len(self.tBreaks)-1): 
                pltInd = range(self.tBreaks[itt], self.tBreaks[itt+1])
                if errors is True:
                    self.ax.errorbar(self.time[pltInd], 
                        self.hires[self.dataKey][pltInd, int(i)], 
                        yerr = yerr[int(i)], fmt='o', ms = 2, 
                        c = c[i])
                else:
                    self.ax.plot(self.time[pltInd], 
                        self.hires[self.dataKey][pltInd, int(i)], c = c[i])
        
        self.ax.set_xlabel('UTC')
        #self.ax.set_xticks(rotation = 30)
        if self.dataKey == 'Col_flux':
            self.ax.set_ylabel(r'Collimated flux $(s \ cm^2 \ sr \ keV)^{-1}$')
        elif self.dataKey == 'Col_counts':
            self.ax.set_ylabel(f'Collimated counts/{cadence} ms')
        elif self.dataKey == 'Sur_counts':
            self.ax.set_ylabel(r'Surface counts')
        elif self.dataKey == 'Sur_flux':
            self.ax.set_ylabel(r'Surface flux $(s \ cm^2 \ sr \ keV)^{-1}$')
        self.ax.set_yscale(yscale)
        self.ax.set_title('FU{} HiRes from {}'.format(self.hiresName[2], self.hiresName[10:20]))
        self.ax.legend(handles=labelHandles, loc='best')

        self.ax.xaxis.set_major_formatter(FuncFormatter(self._format_xaxis))
        self.ax.set_xlabel('time\nL\nMLT\nlat\nlon', 
        fontsize=10)
        self.ax.xaxis.set_label_coords(-0.05,-0.01)
        plt.tight_layout()
        plt.subplots_adjust(bottom=0.20)
        # Start interactive session
        self.ax.format_coord = lambda x, y: '{}, {}'.format(
            mdates.num2date(x).replace(tzinfo=None).isoformat(), round(y))
        self.ax.figure.canvas.mpl_connect('key_press_event', self._mouseTime)
        return

    def plot_map(self, tRange, channel=0):
        """
        This function plots the map of the orbit +/- 10 minutes of the click time.
        """
        idx = np.where((self.time > tRange[0]) & 
                    (self.time < tRange[1]))[0]
        fig = plt.figure(figsize=(12, 6))
        self.mc = fig.canvas
        ax = plt.subplot(111, projection=ccrs.PlateCarree())
        ax.coastlines()
        #ax.stock_img()
        sc = ax.scatter(self.hires['Lon'][idx], self.hires['Lat'][idx], 
                        c=self.hires[self.dataKey][idx, channel],
                    transform=ccrs.PlateCarree(), norm=colors.LogNorm())

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
        ax.text(self.hires['Lon'][idx[0]], self.hires['Lat'][idx[0]], '*',
                ha='center', va='center', fontsize=20, color='red',
                transform=ccrs.PlateCarree())

        fig.suptitle('FU{} HiRes {}\n{} to {}'.format(
                    args.sc_id, tRange[0].date(), 
                    tRange[0].replace(microsecond=0).time(), 
                    tRange[1].replace(microsecond=0).time()), 
                    fontsize=16)
        plt.tight_layout()
        # Plot colorbar
        fig.subplots_adjust(right=0.89, top=0.9, bottom=0.1)
        cbar_ax = fig.add_axes([0.9, 0.15, 0.05, 0.7])
        fig.colorbar(sc, orientation='vertical', cax=cbar_ax,
                label='channel {} [counts/bin]'.format(channel+1))
        return

    def get_hires(self, sc_id, date, level=2):
        # Find the filename of the HiRes file
        self.hiresName = 'FU{}_Hires_{}_L{}.txt'.format(sc_id, date, level)
        print('Loading HiRes file:', self.hiresName)
        # Change the directory, depending where you store your data. 
        if self.hiresName[2] == '3':
            directory = '/home/mike/research/firebird/Datafiles/FU_3/hires/level' + str(level) + '/'
        elif self.hiresName[2] == '4':
            directory = '/home/mike/research/firebird/Datafiles/FU_4/hires/level' + str(level) + '/'
        else:
            raise LookupError('Not a valid spacecraft number!')
        
        # Read in the files.    
        self.hires = dm.readJSONheadedASCII(
                    os.path.join(directory + self.hiresName))
        self.time = spt.Ticktock(self.hires["Time"]).UTC
        self.energyCenters = self.hires['Col_counts'].attrs['ELEMENT_LABELS']
        self.energyRanges = self.hires['Col_counts'].attrs['ENERGY_RANGES']
        return

    
    def _mouseTime(self, event):
        """ 
        # Initialize interactive sesson
        When a user presses 't', this function will print the spacecraft info.
        """
        t = mdates.num2date(ax.get_xlim()[0]).replace(tzinfo=None).replace(microsecond=0)
        save_datetime = t.strftime('%Y%m%d_%H%M')
        canvas.get_default_filename = lambda: '{}_FU{}_hires.png'.format(save_datetime, args.sc_id)

        if event.key == 't':
            time = t.isoformat()
            Lt = self.L[np.argmin(np.abs(event.xdata - mdates.date2num(self.time)))]
            print(time, 'L =', Lt)
        elif event.key == 'm':
            tRange = [t.replace(tzinfo=None) for t in mdates.num2date(self.ax.get_xlim())]
            self.plot_map(tRange)
            self.mc.get_default_filename = lambda: '{}_FU{}_hires_map.png'.format(save_datetime, args.sc_id)
            plt.show()
        return

    def _labels(self):
        ### FORMAT X-AXIS to show more information ###
        self.hires['McIlwainL'] = np.round(self.hires['McIlwainL'], decimals=1)
        self.L = np.copy(self.hires['McIlwainL']).astype(object)
        self.L[self.L > 100] = ''
        # This code is a nifty way to format the x-ticks to my liking.
        self.Labels = ['{}\n{}\n{}\n{}\n{}'.format(
                    t.replace(microsecond=0).time(),
                    L, round(MLT,1), round(lat,1), round(lon,1)) for 
                    (t, L, MLT, lat, lon) in zip(
                    self.time[::10], self.L[::10], self.hires['MLT'][::10], 
                    self.hires['Lat'][::10], self.hires['Lon'][::10])]  
        self.numTimes = mdates.date2num(self.time[::10])
        return

    def _format_xaxis(self, tick_val, tick_pos):
        """
        The tick magic happens here. pyplot gives it a tick time, and this function 
        returns the closest label to that time. Read docs for FuncFormatter().
        """
        dt = self.numTimes-tick_val
        # If time difference between matplotlib's tick and HiRes time stamp 
        # is larger than 30 minutes, skip that tick label.
        if np.min(np.abs(dt)) > 30/1440:
            return ''
        else:
            idx = np.argmin(np.abs(dt))
            return self.Labels[idx]
        
if __name__ == '__main__':
    ### Parse args ###
    parser = argparse.ArgumentParser(description=('This script plots the '
        'FIREBIRD-II HiRes level 2 data. No time correction is applied.'))
    parser.add_argument('sc_id',  type=int,
        help=('This is the FIREBIRD-II unit arument (3 or 4)'))
    parser.add_argument('date', nargs=3, type=int,
        help=('This is the FIREBIRD-II campaign number')) 
    parser.add_argument('-k', '--key', type=str, default='Col_counts',
        help=('This is the detector key (Col_counts, Col_flux, Sur_counts, and Sur_flux)'))     
    args = parser.parse_args()

    # Tap into the magical matplotlib backend to make magic happen.
    fig, ax = plt.subplots(figsize=(8, 6))
    canvas = fig.canvas 
    saveDir = '/home/mike/Dropbox/0_firebird_research/ops/'
    matplotlib.rcParams["savefig.directory"] = saveDir
    
    ### Plot data
    date = date(*args.date)
    hr = plotHiResData(args.sc_id, date)   
    hr.plotTimeSeries(ax=ax, dataKey=args.key)
    plt.show()

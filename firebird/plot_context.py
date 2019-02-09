#!/usr/bin/env python3
"""
This script plots the context data from the FIREBIRD-II mission.
It accepts sc_id (3 or 4) and campaign number arguments, and 
looks for the most recent data in the firebird/Datafiles
firectory.
"""
import spacepy.datamodel
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import dates
import matplotlib.colors as colors
import cartopy.crs as ccrs
from matplotlib.ticker import FuncFormatter, MaxNLocator
import dateutil.parser
import matplotlib
from datetime import datetime, timedelta
import glob
import os
import argparse

plt.style.use('bmh')

### ARGPARSE COMMANDS
parser = argparse.ArgumentParser(description=('This script plots the '
        'FIREBIRD-II Context level 2 data. No time correction is applied.'))
parser.add_argument('sc_id',  type=int,
    help=('This is the FIREBIRD-II unit arument (3 or 4)'))
parser.add_argument('camp', type=int,
    help=('This is the FIREBIRD-II campaign number'))     
args = parser.parse_args()

### FIND DATA ###
dataDir = ('/home/mike/research/firebird/Datafiles/FU_{}/context/'
            'level2'.format(args.sc_id))
fNames = glob.glob(os.path.join(dataDir, 
                    'FU{}_Context_camp{:02}_*.txt'.format(
                    args.sc_id, args.camp)))
# Get the latest campaign file if there are multiple files.
fname = sorted(fNames)[-1]

### READ DATA ###
d = spacepy.datamodel.readJSONheadedASCII(
    os.path.join(dataDir, fname))
d['Time'] = np.array([dateutil.parser.parse(i) for i in d['Time']])

### PLOT CONTEXT ###
fig, ax = plt.subplots(figsize=(12, 8))
# Tap into the magical matplotlib backend to make magic happen.
canvas = fig.canvas 
saveDir = '/home/mike/Dropbox/0_firebird_research/ops/camp_{}'.format(args.camp)
matplotlib.rcParams["savefig.directory"] = saveDir
if not os.path.exists(saveDir): 
    os.makedirs(saveDir)
    print('Made new directory:', saveDir)

ax.plot(d['Time'], d['D0'], 'r', label='D0')
ax.plot(d['Time'], d['D1'], 'b', label='D1')
ax.legend(loc=1)
ax.set(yscale='log', ylabel='Counts/6 s',
    title='FU{} Context'.format(args.sc_id)
    )

### FORMAT X-AXIS to show more information ###
d['McIlwainL'] = np.round(d['McIlwainL'], decimals=1)
L = np.copy(d['McIlwainL']).astype(object)
L[L > 100] = ''
# This code is a nifty way to format the x-ticks to my liking.
Labels = ['{}\n{}\n{}\n{}\n{}\n{}'.format(
            t.replace(microsecond=0).date(), t.replace(microsecond=0).time(),
            L, round(MLT,1), round(lat,1), round(lon,1)) for 
            (t, L, MLT, lat, lon) in zip(
            d['Time'][::10], L[::10], d['MLT'][::10], 
            d['Lat'][::10], d['Lon'][::10])]  
numTimes = dates.date2num(d['Time'][::10])

def format_fn(tick_val, tick_pos):
    """
    The tick magic happens here. pyplot gives it a tick time, and this function 
    returns the closest label to that time. Read docs for FuncFormatter().
    """
    dt = numTimes-tick_val
    # If time difference between matplotlib's tick and HiRes time stamp 
    # is larger than 10 minutes, skip that tick label.
    if np.min(np.abs(dt)) > 10/1440:
        return ''
    else:
        idx = np.argmin(np.abs(dt))
        return Labels[idx]

ax.xaxis.set_major_formatter(FuncFormatter(format_fn))
ax.set_xlabel('date\ntime\nL\nMLT\nlat\nlon')
ax.xaxis.set_label_coords(-0.1,-0.03)
plt.tight_layout()
plt.subplots_adjust(bottom=0.15)

# Initialize interactive sesson
def mouseTime(event):
    """ 
    When a user presses 't', this function will print the spacecraft info.
    """
    # Set default save filename
    t = dates.num2date(ax.get_xlim()[0]).replace(tzinfo=None).replace(microsecond=0)
    save_datetime = t.strftime('%Y%m%d_%H%M')
    canvas.get_default_filename = lambda: '{}_FU{}_context.png'.format(save_datetime, args.sc_id)

    if event.key == 't':
        time = t.isoformat()
        Lt = L[np.argmin(np.abs(event.xdata - dates.date2num(d['Time'])))]
        print(time, 'L =',Lt)
    elif event.key == 'm':
        tRange = [t.replace(tzinfo=None) for t in dates.num2date(ax.get_xlim())]
        plot_map(tRange)
        #time = dates.num2date(event.xdata).replace(tzinfo=None)
        # plot_map(time)
    return

def plot_map(tRange, channel='D0'):
    """
    This function plots the map of the orbit +/- 10 minutes of the click time.
    """
    idx = np.where((d['Time'] > tRange[0]) & 
                   (d['Time'] < tRange[1]))[0]
    fig = plt.figure(figsize=(12, 6))
    ax = plt.subplot(111, projection=ccrs.PlateCarree())
    ax.stock_img()
    sc = ax.scatter(d['Lon'][idx], d['Lat'][idx], c=d[channel][idx],
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
    ax.text(d['Lon'][idx[0]], d['Lat'][idx[0]], '*',
            ha='center', va='center', fontsize=20, color='red',
            transform=ccrs.PlateCarree())

    fig.suptitle('FU{} context {}\n{} to {}'.format(
        args.sc_id, tRange[0].date(), tRange[0].replace(microsecond=0).time(), 
        tRange[1].replace(microsecond=0).time()), fontsize=16)
    plt.tight_layout()
    # Plot colorbar
    fig.subplots_adjust(right=0.89, top=0.9, bottom=0.1)
    cbar_ax = fig.add_axes([0.9, 0.15, 0.05, 0.7])
    fig.colorbar(sc, orientation='vertical', cax=cbar_ax,
            label='{} [counts / 6s]'.format(channel))
    plt.show()
    return

ax.format_coord = lambda x, y: '{}, {}'.format(
            dates.num2date(x).replace(tzinfo=None).isoformat(), round(y))
ax.figure.canvas.mpl_connect('key_press_event', mouseTime)
plt.show()

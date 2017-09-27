# Plot EMFISIS data.
#import os
from spacepy import pycdf
import spacepy.datamodel
import matplotlib.pylab as plt
#from matplotlib import colors, ticker, cm, gridspec
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import matplotlib.ticker
import matplotlib.dates
#import matplotlib.gridspec as gridspec
import numpy as np
import os, sys, glob, copy
import dateutil.parser
from datetime import datetime#, timedelta

# Try to import the directorories module
try:
    sys.path.insert(0, '/home/mike/Dropbox/0_grad_work')
    import directories
except SystemError:
    print('Could not import the directories.py file, '
        'please supply data directories manualy!')

m_e = 9.1E-31 # kg
q_e = -1.6E-19 # coulombs

class EMFISISspectra:
    def __init__(self, sc_id, date, **kwargs):
        """
        NAME:    EMFISISspectra(sc_id, date, **kwargs)
        USE:     This function will plot emfisis data. There are helper 
                 functions that load in the spectra data, and the local 
                 (observed) magnetic field data or the minimum B from a magnetic
                 field model from a magnetic ephemeris data file.
        INPUT:   REQUIRED:
                    sc_id either 'A' or 'B', date is a 
                    datetime object of a date to plot. 
                 OPTIONAL:
                    ### FINISH WRITING THIS LATER ###
                    wfrDir      - WFR spectra directory
                    bDir        - Magnetic field measurments directory
                    magEphemDir - Magnetic ephemeris directory
                    tBounds     - Time bounds for spectra plotting

        AUTHOR:  Mykhaylo Shumko
        RETURNS: None
        MOD:     2017-06-27
        """
        self.sc_id = sc_id
        self.date = date
        
        # Create a merged date string for file lookup.
        splitDate = date.date().isoformat().split('-')
        self.joinedDate = ''.join(splitDate)

        # Can set the data directory parameters here
        self.wfrDir = kwargs.get('wfrDir', directories.emfisis_dir(sc_id))
        self.hfrDir = kwargs.get('hfrDir', directories.emfisis_dir(sc_id))
        self.bDir = kwargs.get('bDir', directories.emfisis_dir(sc_id))
        self.magEphemDir = kwargs.get('magEphemDir', directories.rbsp_magephem_dir(sc_id))
        self.tBounds = kwargs.get('tBounds', None)
        return

    def loadWFRSpectra(self):
        """
        NAME:    loadWFRSpectra(self)
        USE:     This function will load EMFISIS WFR data.
        INPUT:   None
        AUTHOR:  Mykhaylo Shumko
        RETURNS: None, self.wfrSpec with the magnitude of the diagonal 
                 magnetic field elements in array as 'spectra' key
        MOD:     2017-06-27
        """
        wfrName = sorted(glob.glob(os.path.join(self.wfrDir,
            ('*rbsp-' + self.sc_id.lower() + 
            '_WFR-spectral-matrix-diagonal_emfisis*' + 
            self.joinedDate + '*'))))

        assert len(wfrName) == 1, ('Error, none or multiple'+
        ' WFR files found! \n {} \n {}'.format(wfrName, self.wfrDir))
        wfrName = wfrName[0]
        self.spec = pycdf.CDF(wfrName).copy() # copy so that it can be modified.
        
        # Now find filtered data indicies if tBounds is specified.
        if self.tBounds is not None:
            wfrInd = np.where((self.spec['Epoch'][:] >= self.tBounds[0]) &
                (self.spec['Epoch'][:] <= self.tBounds[1]))[0]
            assert len(wfrInd) > 0, ('ERROR: no filtered spectra found in the time'
            ' range specified! Check tBounds keyword')
        else:
            wfrInd = range(len(self.spec['BuBu']))

        self.spec['Epoch'] = self.spec['Epoch'][wfrInd]
        self.spec['Spectra'] = np.sqrt(self.spec['BuBu'][wfrInd]**2 + 
            self.spec['BvBv'][wfrInd]**2 + self.spec['BwBw'][wfrInd]**2)
        return 
        
    def loadHFRSpectra(self):
        """
        NAME:    loadHFRSpectra(self)
        USE:     This function will load EMFISIS HFR data.
        INPUT:   None
        AUTHOR:  Mykhaylo Shumko
        RETURNS: None, self.wfrSpec with the magnitude of the diagonal 
                 magnetic field elements in array as 'spectra' key
        MOD:     2017-06-27
        """
        hfrName = sorted(glob.glob(os.path.join(self.hfrDir,
            ('*rbsp-' + self.sc_id.lower() + 
            '_HFR-spectra_emfisis*' + self.joinedDate + '*'))))

        assert len(hfrName) == 1, ('Error, none or multiple'+
        ' WFR files found! \n {} \n {}'.format(hfrName, self.hfrDir))
        hfrName = hfrName[0]
        self.spec = pycdf.CDF(hfrName).copy() # copy so that it can be modified.
        
        # Now find filtered data indicies if tBounds is specified.
        if self.tBounds is not None:
            hfrInd = np.where((self.spec['Epoch'][:] >= self.tBounds[0]) &
                (self.spec['Epoch'][:] <= self.tBounds[1]))[0]
            assert len(hfrInd) > 0, ('ERROR: no filtered spectra found in the time'
            ' range specified! Check tBounds keyword')
        else:
            hfrInd = range(len(self.spec['Epoch']))

        self.spec['Epoch'] = self.spec['Epoch'][hfrInd]
        self.spec['Spectra'] = self.spec['HFR_Spectra'][hfrInd, :]
        return     

    def loadMagEphem(self, Bmodel = 'TS04D'):
        """
        This function will load the magnetic ephemeris file and filter it by the
        time bounds, if necessary. This will also calculate the model gyrofrequency
        at the magnetic equator.
        """
        """
        NAME:    loadMagEphem(self, Bmodel = 'TS04D')
        USE:     This function will load the magnetic ephemeris file and 
                 filter it by the time bounds, if necessary. This will also
                 calculate the model equatorial gyrofrequency as self.feEq
       
        INPUT:   REQUIRED:
                    None
                 OPTIONAL:
                    Bmodel - Magnetic field model magEphem file to use
        AUTHOR:  Mykhaylo Shumko
        RETURNS: None
        MOD:     2017-06-27
        """
        magEphemName = sorted(glob.glob(os.path.join(self.magEphemDir,
            ('rbsp{}*MagEphem*{}*{}*.txt'.format(self.sc_id.lower(), 
            Bmodel, self.joinedDate)))))
        assert len(magEphemName) == 1, ('Error, none or multiple'
        ' magnetic ephemeris files found!')
        magEphemName = magEphemName[0]
        self.magEphem = spacepy.datamodel.readJSONheadedASCII(magEphemName)
        
        # Convert times
        self.magEphem['DateTime'] = np.array(list(map(
            lambda x: dateutil.parser.parse(x).replace(tzinfo=None),
            self.magEphem['DateTime'])))

        # Filter data by times
        if self.tBounds is not None:
            magInd = np.where((self.magEphem['DateTime'] >= self.tBounds[0]) &
                (self.magEphem['DateTime'] <= self.tBounds[1]))[0]
            assert len(magInd) > 0, ('ERROR: no filtered magEphem found in the time'
            ' range specified! Check tBounds keyword')
        else:
            magInd = range(len(self.magEphem['DateTime']))

        self.magEphem['DateTime'] = self.magEphem['DateTime'][magInd]
        self.magEphem['Bmin'] = self.magEphem['Bmin_gsm'][magInd, 3]
        self.magEphem['Lstar'] = self.magEphem['Lstar'][magInd, :]
        self.magEphem['EDMAG_MLT'] = self.magEphem['EDMAG_MLT'][magInd]
        self.magEphem['EDMAG_MLAT'] = self.magEphem['EDMAG_MLAT'][magInd]
        self.feEq = np.abs(q_e)*np.array(self.magEphem['Bmin'])/(
            np.power(10, 9)*2*np.pi*m_e)
        return

    def plotSpectra(self, ax = None, pltFe = True, **kwargs):
        """
        NAME:    EMFISISspectra(sc_id, date, **kwargs)
        USE:     This function will plot emfisis data. Optionally, it can 
                 superpose and the equatorial electron gyrofrequency from the 
                 magEphem data.
        INPUT:   REQUIRED:
                    None
                 OPTIONAL:
                    ax - Subplot object (Will create one otherwise)
                    pltFe - Plot equatorial electron gyrofrequency. (True)
                    plotCb - Plot colorbar ('vertical', 'v', 'horizontal',
                            'h', False)
        AUTHOR:  Mykhaylo Shumko
        RETURNS: None
        MOD:     2017-06-27
        """
        dcLevels = kwargs.get('dcLevels', 0.5)
        printTitle = kwargs.get('printTitle', True)
        cbX = kwargs.get('cbX', [-0.05, -0.2])
        plt.rcParams['font.size'] = kwargs.get('fontSize', 12)
        grid = kwargs.get('grid', True)
        plotCb = kwargs.get('plotCb', 'horizontal')
        printXlabel = kwargs.get('printXlabel', True)
        cAspect = kwargs.get('cAspect', 100)
        cax = kwargs.get('cax', None)
        legendLoc = kwargs.get('legendLoc', 1)
        instrument = kwargs.get('instrument', 'WFR')

        if ax is None:
            fig = plt.figure(figsize=(15, 10), dpi=80, facecolor = 'white')
            gs = gridspec.GridSpec(1,1)
            self.ax = fig.add_subplot(gs[0, 0], facecolor='k')
        else:
            self.ax = ax
            
        if instrument == 'WFR':
            freq = 'WFR_frequencies'
            lowF = kwargs.get('lowF', 100)
            highF = kwargs.get('highF', 10**4)
            vmin = kwargs.get('spectraMin', 10**-10)
            vmax = kwargs.get('spectraMax', 10**-5)
            levels = kwargs.get('cLevels', 
                1/np.power(10, np.arange(10, 5, -dcLevels)))
        else:
            freq = 'HFR_frequencies'
            lowF = kwargs.get('lowF', 10**4)
            highF = kwargs.get('highF', 5*10**5)
            vmin = kwargs.get('spectraMin', 10**-30)
            vmax = kwargs.get('spectraMax', 10**-8)
            levels = levels = kwargs.get('cLevels', 
                1/np.power(10, np.arange(30, 8, -dcLevels)))
            
        ### Create 2D mesh grid to visualize. ###
        wTT, wFF = np.meshgrid(range(len(self.spec['Epoch'])),
            np.array(self.spec[freq])[0, :])
        # Fill wTT with datetimes instead of indicies
        wTT = np.broadcast_to(self.spec['Epoch'][:], 
            (len(self.spec[freq][0, :]), len(self.spec['Epoch'][:]) ))

        # Set up x and y axis ticks for imshow keyword extent. Have to 
        # convert datetimes.
        ticks = [matplotlib.dates.date2num(self.spec['Epoch'][0]), 
                matplotlib.dates.date2num(self.spec['Epoch'][-1]),
                np.array(self.spec[freq])[0, 0], 
                np.array(self.spec[freq])[0, -1]]
        #self.ax.set_yscale('log')
        cs = self.ax.pcolormesh(wTT, wFF, np.transpose(self.spec['Spectra']),
             cmap = plt.get_cmap('gnuplot2'), norm=colors.LogNorm(),
             vmin = vmin, vmax = vmax)
        if vmin is not None:
            cs.set_clim(vmin = vmin)
        if vmax is not None:
            cs.set_clim(vmax = vmax)
    
        # Format time stamps
        self.ax.xaxis_date()
        fmtr = matplotlib.dates.DateFormatter('%H:%M:%S')
        self.ax.xaxis.set_major_formatter(fmtr)

        # Scientific notation on y-axis
        self.ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1e'))

        
        ### Plot electron gyrofrequencies ###
        if pltFe:
            self.ax.plot(self.magEphem['DateTime'], self.feEq, '--r', 
                label = r'$f_{e \ Eq}$', lw = 2)
            self.ax.plot(self.magEphem['DateTime'], self.feEq/2, '--g', 
                label = r'$f_{e \ Eq}/2$', lw = 2)
            self.ax.plot(self.magEphem['DateTime'], self.feEq/10, '--c', 
                label = r'$f_{e\ Eq}/10$', lw = 2)
        
        ### Set all sorts of plot parameters ###
        self.ax.set_ylim(top = highF, bottom = lowF)
            
        self.ax.legend(loc = legendLoc, bbox_to_anchor = None, fontsize = 10)
        if printTitle: 
            if instrument == 'WFR':
                self.ax.set_title('RBSP{} EMFISIS {} B field spectra '.format
                    (self.sc_id.upper(), instrument)
                    + '(diagonal elements) from {}'.format(
                    self.date.date().isoformat()))    
            else:
                self.ax.set_title('RBSP{} EMFISIS {} B field spectra '.format
                (self.sc_id.upper(), instrument) + ' from {}'.format(
                self.date.date().isoformat()))
        self.ax.set(ylabel = 'Hz')
        if printXlabel:
            self.ax.set_xlabel('UTC')

        self.ax.set_yscale('log')

        ### PLOT COLORBAR ###
        if plotCb in ['horizontal', 'h']:
            self.cb = plt.colorbar(cs, ticks = matplotlib.ticker.LogLocator(), 
                ax = self.ax, orientation='horizontal', cax = cax)
            # Move label to left to get out of the way of other plots.
            self.cb.ax.xaxis.set_label_coords(cbX[0], cbX[1], transform =
                self.ax.transAxes)

        elif plotCb in ['vertical', 'v']:
            self.cb = plt.colorbar(cs, ticks = matplotlib.ticker.LogLocator(), 
                ax = self.ax, orientation='vertical', cax = cax)

        self.cb.set_label(r'$nT^2/Hz$')
        self.cb.set_clim(vmin, vmax)

        # Plot the grid
        if grid:
            self.ax.grid(b=True, which='major', color='w', linestyle='-')
            self.ax.grid(b=True, which='minor', color='w', linestyle='--')
        if ax is None:
            gs.tight_layout(fig)
            plt.show()   
        return

if __name__ == '__main__':
    date = datetime(2017, 3, 31)
    #tBounds = [datetime(2017, 3, 31, 11, 10), datetime(2017, 3, 31, 11, 20)]
    tBounds = [datetime(2017, 3, 31, 11, 10), datetime(2017, 3, 31, 11, 25)]
    sc_id = 'A'

    pObj = EMFISISspectra(sc_id, date, tBounds = tBounds)
    pObj.loadWFRSpectra()
    pObj.loadMagEphem()
    
    #pObj.plotSpectra(plotCb = 'vertical', instrument = 'HFR')
    pObj.plotSpectra(plotCb = 'vertical', spectraMin = 10**-11, 
        spectraMax = 10**-8)
        
    #cLevels = np.power(10, -np.arange(9, 7, -0.1))

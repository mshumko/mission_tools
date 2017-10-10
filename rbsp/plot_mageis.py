import numpy as np
from datetime import datetime, timedelta
from spacepy import pycdf
import spacepy.datamodel
import glob, os, sys
import dateutil.parser
import matplotlib.pylab as plt
import matplotlib.colors as colors
import matplotlib.ticker
from matplotlib.patches import Rectangle
from matplotlib.pyplot import cm 
import matplotlib.gridspec as gridspec
import scipy.stats

# Try to import the directorories module
try:
    sys.path.insert(0, '/home/mike/Dropbox/0_grad_work')
    import directories
except SystemError:
    print('Could not import the directories.py file, '
        'please supply data directories manualy!')

class magEISspectra:
    def __init__(self, sc_id, date, **kwargs):
        """
        NAME:    magEISspectra(sc_id, date, **kwargs)
        USE:     This function will plot mageis flux.
        INPUT:   REQUIRED:
                    sc_id either 'A' or 'B', date is a 
                    datetime object of a date to plot. 
                 OPTIONAL:
                    ### FINISH WRITING THIS LATER ###

        AUTHOR:  Mykhaylo Shumko
        RETURNS: None
        MOD:     2017-07-04
        """
        self.sc_id = sc_id
        self.date = date
        splitDate = date.date().isoformat().split('-')
        self.joinedDate = ''.join(splitDate)

        self.dataLevel = kwargs.get('dataLevel', 2)
        self.mageisDir = kwargs.get('mageisDir', directories.mageis_dir(sc_id))
        self.magEphemDir = kwargs.get('magEphemDir', 
            directories.rbsp_magephem_dir(sc_id))
        self.tBounds = kwargs.get('tBounds', None) 
        return

    def loadMagEIS(self, **kwargs):
        """
        NAME:    loadMagEIS(self)
        USE:     This function will find and load the correct magEIS data, and 
                 filter by times if specified.
        INPUT:   OPTIONAL:
                 relType = 'int' - Release type (internal or public)
                 instrument = '' - Instrument to load, ('', 'LOW', 'M35', 
                 'M75', 'HIGH')
        AUTHOR:  Mykhaylo Shumko
        RETURNS: self.magEIS dictionary
        MOD:     2017-07-04
        """
        # Optional variables.
        relType = kwargs.get('relType', 'int')
        instrument = kwargs.get('instrument', '')
        self.highrate = kwargs.get('highrate', False)
        
        # Change data keys and file string lookup, depending on if you are 
        # plotting the highrate or not.
        if self.highrate:
            self.alphaKey = 'HighRate_Alpha360'
            self.fluxKey = 'HighRate'
            highrate = 'hr'
        else:
            self.alphaKey = 'FEDU_Unbinned_Alpha360'
            self.fluxKey = 'FEDU_Unbinned_0to360'
            highrate = ''

        # Find the file.
        searchStr = os.path.join(self.mageisDir, 
            'rbsp{}_{}_ect-mageis{}-*{}*L{}_{}*.cdf'.format(self.sc_id.lower(), 
            relType, instrument, highrate, self.dataLevel, self.joinedDate))
        magEISname = sorted(glob.glob(searchStr))
        assert len(magEISname) == 1, ('Error, none or multiple magEIS files '
            'found! \n Search str: {} \n {}'.format(searchStr, magEISname))
        magEISname = magEISname[0]
        # Load and copy the data so that it can be modified.
        self.magEISdata = pycdf.CDF(magEISname).copy() 
        
        # Now find filtered data indicies if tBounds is specified.
        if self.tBounds is not None:
            magEISInd = np.where((self.magEISdata['Epoch'][:] >= self.tBounds[0]) &
                (self.magEISdata['Epoch'][:] <= self.tBounds[1]))[0]
            assert len(magEISInd) > 0, ('ERROR: no filtered spectra found in '
                'the time range specified! Check tBounds keyword')
        else:
            magEISInd = range(len(self.magEISdata['Epoch']))

        N = len(self.magEISdata['Epoch'])

        # Filter the data.
        for key in self.magEISdata.keys():
            # If key is a function of time, filter it.
            if len(self.magEISdata[key]) == N:
                self.magEISdata[key] = self.magEISdata[key][magEISInd]
                
                
        # Define magEIS detector constants
        if instrument.lower() == 'low':
            if self.sc_id.upper() == 'A':
                self.Emid = [34, 54, 78, 108, 143, 182, 223] # keV
                self.Elow = [29, 46, 68, 95, 126, 164, 206] # keV
                self.Ehigh = [41, 66, 92, 126, 164, 204, 247] # keV
                # Units of (keV cm^2 sr)
                self.G0dE = [4.13E-2, 5.73E-2, 6.056E-2, 6.88E-2, 7.35E-2, 
                    6.90E-2, 5.98E-2]
                self.Ebins = [29, 41, 66, 92, 126, 164, 204, 247]
                
            if self.sc_id.upper() == 'B':
                self.Emid = [32, 51, 74, 101, 132, 168, 208] # keV
                self.Elow = [27, 43, 63, 88, 117, 152, 193] # keV
                self.Ehigh = [39, 63, 88, 117, 150, 188] # keV
                # Units of (keV cm^2 sr)
                self.G0dE = [4.33E-2, 5.41E-2, 5.926E-2, 6.605E-2, 6.460E-2,
                    6.23E-2, 5.96E-2]
                self.Ebins = [27, 39, 63, 88, 117, 150, 188]

        return self.magEISdata

    def calcDailyElectronEnergies(self, method = 'mode'):
        """
        NAME:    calcDayEnergies(self, method = 'mean')
        USE:     This function will add a new key, the energy channels for
                 electrons
        INPUT:   REQUIRED:
                    None
                 OPTIONAL:
                    method - How to calculate the enrgies, 'mean', 'median', 'mode'
        AUTHOR:  Mykhaylo Shumko
        RETURNS: self.magEphem['eEnergy']
        MOD:     2017-07-04

        ### MAYBE MAKE THIS MORE GENERAL WITH ENERGY DELTAS/WIDTHS? ###
        """
        assert method in ['mean', 'median', 'mode'], ('ERROR, statistical method'
        ' to calculate the electron energies is not "mean", "median", or "mode".')

        # Find error values and replace with nan
        validInd = np.where(self.magEISdata['FEDU_Energy'] == -1E31)
        self.magEISdata['FEDU_Energy'][validInd] = np.nan
    
        if method is 'mean':
            self.magEISdata['eEnergy'] = np.nanmean(self.magEISdata['FEDU_Energy'], 
                axis = 0)
        elif method is 'median':
            raise ValueError('Still not implemented!')
            self.magEISdata['eEnergy'] = np.nanmedian(self.magEISdata['FEDU_Energy'])
        elif method is 'mode':
            self.magEISdata['eEnergy'] = scipy.stats.mode(
                self.magEISdata['FEDU_Energy'], axis = 0, nan_policy = 'omit')[0][0]
        return self.magEISdata['eEnergy']    
        
    def resolveSpinTimes(self, spin_thresh = 0.1, est_spin = 10.9, n_sectors = None, flattenTime = False):
        """
        This funciton creates an array of times by spin sector.
        """
        # If number of pitch angle sectors are not specified, find it.
        if n_sectors is None:
            if self.highrate:
                n_sectors = 1000
            else:
                n_sectors = 64
            
        # If self.times array is not created, make it
        if not hasattr(self, "times"):
            self.times = np.repeat(self.magEISdata['Epoch'][:, np.newaxis],
                self.magEISdata[self.alphaKey].shape[1] ,axis = 1)
                 
        # Find valid Indicies. Minus one is to calculate dTspin.
        for i in range(self.magEISdata[self.alphaKey].shape[0]-1):
            # Logic that checks for data gaps. If dTSpin is outside 
            # est_spin +/-(1+spin_thresh), assume est_spin spin period
            dTspin = (self.times[i+1, 0] - self.times[i, 0]).total_seconds()
            if ((dTspin > est_spin*(1+spin_thresh)) & 
                (dTspin < est_spin*(1-spin_thresh)) ):
                dTspin = est_spin
                
            self.times[i, :n_sectors] = np.array([self.times[i, iAlpha] + 
                timedelta(seconds = iAlpha*dTspin/n_sectors) for iAlpha in 
                range(n_sectors)])
        # This just interpolates the last spin assuming the previous spin
        self.times[-1, :n_sectors] = np.array([self.times[-1, iAlpha] + 
                timedelta(seconds = iAlpha*dTspin/n_sectors) for iAlpha in 
                range(n_sectors)]) 
                 
        if flattenTime: # Flatten times for time series plotting.
            self.times = self.times.flatten()
        return self.times
        
        
    def plotHighRateTimeSeries(self, **kwargs):
        """
        NAME:    plotHighRateTimeSeries(self, **kwargs)
        USE:     This function plots the hight rate magEIS flux as a timeseries
        INPUT:   REQUIRED:
                    None
                 OPTIONAL:
                    'ax'= None: subplot object to plot on. If not specified, 
                        will create own figure.
                    'E_ch'= 0: Energy channel to plot
                    'downsampleAlpha' = 1: Plot the downsampleAlpha'th pitch 
                        angle for quicker plotting on slower machines.
                    'deflatTime' = True: Refromat the 'epoch' array so that each
                        pitch angle is attributed to a unique time stamp. Saved 
                        self.times
        AUTHOR:  Mykhaylo Shumko
        RETURNS: 
        MOD:     2017-07-04
        """
        est_spin = kwargs.get('est_spin', 10.9) # Seconds
        n_sectors = kwargs.get('n_sectors', None)
        spin_thresh = kwargs.get('spin_thresh', 0.1) # Fraction of the estimated spin period.
        ax = kwargs.get('ax', None)
        E_ch = kwargs.get('E_ch', None)
        deflatTime = kwargs.get('deflatTime', True) # alphaSpinTimes
        smooth = kwargs.get('smooth', 1)
        
        
        if n_sectors is None:
            if self.sc_id.upper() == 'A':
                n_sectors = 1000
            else:
                n_sectors = 64
        
        if ax is None:
            fig = plt.figure(figsize=(15, 10), dpi=80, facecolor = 'white')
            gs = gridspec.GridSpec(1,1)
            self.bx = fig.add_subplot(gs[0, 0], facecolor='white')
        else:
            self.bx = ax
        
        self.times = np.repeat(self.magEISdata['Epoch'][:, np.newaxis],
            self.magEISdata[self.alphaKey].shape[1] ,axis = 1)
            
        # Deflatten the 1d time array into 2d where each pitch angle is given 
        # a unique time stamp to where the spacecraft was pointing at the time.
        if deflatTime:
            self.resolveSpinTimes(spin_thresh, est_spin, n_sectors)
        
        if E_ch is None:
            E_ch = np.arange(7)
        elif not isinstance(E_ch, (list, np.ndarray)):
            E_ch = [E_ch]
            
        for ee in E_ch:
            flux = self.magEISdata[self.fluxKey][:, :n_sectors, ee].flatten()
            
            # Now smooth the flux
            flux = np.convolve(flux, np.ones(smooth)/smooth, mode='same')
            
            validF = np.where(flux != -1E31)[0]
            flatT = self.times[:, :n_sectors].flatten()
            
            self.bx.plot(flatT[validF], flux[validF], 
                label='{}-{} keV'.format(self.Elow[ee], self.Ehigh[ee]))
        self.bx.set(yscale='log')
        if ax is None:
            plt.show()    
        return self.bx

    def plotHighRateSpectra(self, **kwargs):
        """
        NAME:    plotHighRateSpectra()
        USE:     This function plots the hight rate magEIS flux.
        INPUT:   REQUIRED:
                    None
                 OPTIONAL:
                    'ax'= None: subplot object to plot on. If not specified, 
                        will create own figure.
                    'E_ch'= 0: Energy channel to plot
                    'scatterS' = 10: Size of scatter symbols to plot.
                    'vmin' = None: Min flux value range for colorbar
                    'vmax' = None: Max flux value range for colorbar
                    'pltXlabel' = True: Wether or not to label the xaxis as UTC.
                    'pltTitle' = True: If True, will print RBSP-{} magEIS 
                        HigHRate Data.
                    'downsampleAlpha' = 1: Plot the downsampleAlpha'th pitch 
                        angle for quicker plotting on slower machines.
                    'deflatTime' = True: Refromat the 'epoch' array so that each
                        pitch angle is attributed to a unique time stamp. Saved 
                        self.times
        AUTHOR:  Mykhaylo Shumko
        RETURNS: 
        MOD:     2017-07-04
        """
        ax = kwargs.get('ax', None)
        E_ch = kwargs.get('E_ch', 0)
        scatterS = kwargs.get('scatterS', 10)
        vmin = kwargs.get('vmin', None)
        vmax = kwargs.get('vmax', None)
        pltXlabel = kwargs.get('pltXlabel', True)
        pltTitle = kwargs.get('pltTitle', True)
        downsampleAlpha = kwargs.get('downsampleAlpha', 1)
        deflatTime = kwargs.get('deflatTime', True) # alphaSpinTimes
        spin_thresh = kwargs.get('spin_thresh', 0.1) # Fraction of the estimated spin period.
        est_spin = kwargs.get('est_spin', 10.9) # Seconds
        n_sectors = kwargs.get('n_sectors', None)
        cax = kwargs.get('cax', None)
        aspect = kwargs.get('cAspect', None)
        plotCb = kwargs.get('plotCb', True)
        
        if n_sectors is None:
            if self.sc_id.upper() == 'A':
                n_sectors = 1000
            else:
                n_sectors = 64
        
        if ax is None:
            fig = plt.figure(figsize=(15, 10), dpi=80, facecolor = 'white')
            gs = gridspec.GridSpec(1,1)
            self.ax = fig.add_subplot(gs[0, 0], facecolor='black')
        else:
            self.ax = ax
        
        self.times = np.repeat(self.magEISdata['Epoch'][:, np.newaxis],
            self.magEISdata[self.alphaKey].shape[1] ,axis = 1)
            
        # Deflatten the 1d time array into 2d where each pitch angle is given 
        # a unique time stamp to where the spacecraft was pointing at the time.
        if deflatTime:
            self.resolveSpinTimes(spin_thresh, est_spin, n_sectors)
        
        flux = self.magEISdata[self.fluxKey][:, :, E_ch]
        
        self.sc = self.ax.scatter(self.times[:, ::downsampleAlpha], 
            self.magEISdata[self.alphaKey][:, ::downsampleAlpha], 
            c = flux[:, ::downsampleAlpha], norm=matplotlib.colors.LogNorm(),
            cmap = plt.get_cmap('rainbow'), 
            vmin = vmin, vmax = vmax, s = scatterS)
        if plotCb:
            if aspect is None:
                self.cb = plt.colorbar(self.sc, ax = self.ax, cax = cax,
                orientation='vertical', label = 'counts/sec')
            else:
                self.cb = plt.colorbar(self.sc, ax = self.ax, cax = cax,
                orientation='vertical', aspect = aspect, 
                    label = 'counts/sec')

        self.ax.set_xlim(self.times[0,0], self.times[-1, -1])
        self.ax.set_ylim(0, 360)
        
        self.ax.set_ylabel('Pitch angle (deg)')
        if pltXlabel: self.ax.set_xlabel('UTC')
        if pltTitle: self.ax.set_title('RBSP-{} magEIS HighRate'
            ' Data'.format(self.sc_id.upper()))
        
        # Set pitch angle ticks
        alpha_tick_locator = matplotlib.ticker.FixedLocator(range(0, 360, 10))
        self.ax.yaxis.set_minor_locator(alpha_tick_locator)
        self.ax.tick_params(which = 'minor')
        
        if ax is None:
            plt.show()
        return self.ax     

    def plotSpinAvgFlux(self, ax = None, **kwargs):
        """
        NAME:    plotSpinAvgFlux(self, ax = None, **kwargs)
        USE:     Plots the spin averaged electron flux from magEIS in line format.
        INPUT:   REQUIRED:
                    None
                 OPTIONAL:
                    ax = None: A subplot object. If none specied, it will 
                        create one.
                    yscale = 'log': y scale of the plot.
                    ymin = 1: lower flux bound.
                    pltLegendLoc = 1: If False, will not plot legend, otherwise 
                        will plot in the location given by pltLegendLoc as 
                        documented in matplotlib.pylab.
                    pltXlabel = True: Plot or not the x label
                    pltTitle = True: Plot or not the title.
        AUTHOR:  Mykhaylo Shumko
        RETURNS: Subplot object, self.ax
        MOD:     2017-07-06
        """
        yscale = kwargs.get('yscale', 'log')
        ymin = kwargs.get('ymin', 1)
        pltLegendLoc = kwargs.get('pltLegendLoc', 1)
        pltXlabel = kwargs.get('pltXlabel', True)
        pltTitle = kwargs.get('pltTitle', True)

        if ax is None:
            fig = plt.figure(figsize=(15, 10), dpi=80, facecolor = 'white')
            gs = gridspec.GridSpec(1,1)
            self.ax = fig.add_subplot(gs[0, 0], facecolor='w')
        else:
            self.ax = ax
        c = ['r', 'b', 'g', 'c', 'm', 'k']

        # If the energy channel is errorous, skip it.
        nE = np.where(~np.isnan(self.magEISdata['eEnergy']))[0]
        c = iter(cm.rainbow(np.linspace(0, 1, len(nE))))
        # Plot the sping-averaged electron flux (FESA)
        for E in nE:
            # Find the non-error flux values
            validInd = np.where(self.magEISdata['FESA'][:, E] != -1E31)
            self.ax.plot(self.magEISdata['Epoch'][validInd[0]], 
                self.magEISdata['FESA'][validInd[0], E], c = next(c), 
                label = '{} keV'.format(int(self.magEISdata['eEnergy'][E])))
        
        if pltLegendLoc:
            self.ax.legend(loc = pltLegendLoc)
        self.ax.set_yscale(yscale)
        self.ax.set_ylim(bottom = ymin)
        self.ax.set_ylabel(r'Electron flux $(cm^2 \ sr \ s \ keV)^{-1}$')
        if pltTitle:
            self.ax.set_title('RBSP{} magEIS from {} '.format
                (self.sc_id.upper(), self.date.date().isoformat()))
        if pltXlabel:
            self.ax.set_xlabel('UTC')
        if ax is None:
            plt.show()
        return self.ax
        
    def plotUnidirectionalFlux(self, alpha, energy = None, ax = None, **kwargs):
        """
        NAME:    plotUnidirectionalFlux(self, ax = None)
        USE:     Plots the unidirectional electron flux from magEIS as a lineplot.
        INPUT:   REQUIRED:
                    None
                 OPTIONAL:
                    ax = None: A subplot object. If none specied, it will 
                        create one.
                    yscale = 'log': y scale of the plot.
                    ymin = 1: lower flux bound.
                    pltLegendLoc = 1: If False, will not plot legend, otherwise 
                        will plot in the location given by pltLegendLoc as 
                        documented in matplotlib.pylab.
                    pltXlabel = True: Plot or not the x label
                    pltTitle = True: Plot or not the title.
        AUTHOR:  Mykhaylo Shumko
        RETURNS: Subplot object, self.ax
        MOD:     2017-07-06
        """
        yscale = kwargs.get('yscale', 'log')
        ymin = kwargs.get('ymin', 1)
        pltLegendLoc = kwargs.get('pltLegendLoc', 1)
        pltXlabel = kwargs.get('pltXlabel', True)
        pltTitle = kwargs.get('pltTitle', True)

        if ax is None:
            fig = plt.figure(figsize=(15, 10), dpi=80, facecolor = 'white')
            gs = gridspec.GridSpec(1,1)
            self.ax = fig.add_subplot(gs[0, 0], facecolor='w')
        else:
            self.ax = ax
            
        # Find indicies closest to user alpha input
        self._findAlphaIdx(alpha)
        
        if energy is None:
            self.energyIdx = np.arange(len(self.magEISdata['eEnergy']))
        else:
            self._findEnergyIdx(energy)
            
        # Replace the error values with nans
        self.magEISdata['FEDU_0to180'][self.magEISdata['FEDU_0to180'] == -1e31] = np.nan
        
        # Plot all of the energies and pitch angles.
        for e in self.energyIdx:
            for a in self.alphaIdx:
                self.ax.plot(self.magEISdata['Epoch'][:], 
                    self.magEISdata['FEDU_0to180'][:, a, e],
                    label = r'E = {} keV, $\alpha$ = {}'.format(
                    int(self.magEISdata['eEnergy'][e]), 
                    int(self.magEISdata['FEDU_Alpha'][a])) )
        if pltLegendLoc:
            self.ax.legend(loc = pltLegendLoc)    
        if pltXlabel:
            self.ax.set_xlabel('UTC')
        self.ax.set_ylabel(r'Electron flux $(cm^2 \ sr \ s \ keV)^{-1}$')  
        self.ax.set_yscale(yscale)   
        
        if ax is None:
            plt.show()
        
    def _findAlphaIdx(self, alpha):
        """
        NAME:    _findAlphaIdx(self, alpha)
        USE:     Finds the indicies of pitch angles closest to the pitch angles
                 specified by the alpha array. 
        INPUT:   REQUIRED:
                    alpha: an array of input 
        AUTHOR:  Mykhaylo Shumko
        RETURNS: alphaIdx: Pitch angle index array
        MOD:     2017-07-06
        """
        # Copy the FEDU_Alpha array and insert nan's for error values of -1E31.
        FEDU_Alpha = np.array(self.magEISdata['FEDU_Alpha'])
        FEDU_Alpha[FEDU_Alpha == -1E31] = np.nan
        
        if isinstance(alpha, int): # If alpha is an int
            # find indicies of the FEDU_Alpha closest to user alpha input.
            idx = np.nanargmin(np.abs(FEDU_Alpha - alpha))
            self.alphaIdx = np.array([idx])
        else: # If alpha is an array.
            self.alphaIdx = -1*np.ones(len(alpha), dtype = np.int)
            # Iterate, and do the same as the other case in the if statement.
            for i in range(len(self.alphaIdx)):
                self.alphaIdx[i] = np.nanargmin(np.abs(FEDU_Alpha - alpha[i]))
        return self.alphaIdx
        
    def _findEnergyIdx(self, energy):
        """
        NAME:    _findAlphaIdx(self, alpha)
        USE:     Finds the indicies of energies closest to the energies
                 specified by the 'eEnergy' array
        INPUT:   REQUIRED:
                    alpha: an array of input 
        AUTHOR:  Mykhaylo Shumko
        RETURNS: alphaIdx: Energy index array
        MOD:     2017-07-06
        """
        # Copy the FEDU_Alpha array and insert nan's for error values of -1E31.
        FEDU_energy = np.array(self.magEISdata['eEnergy'])

        if isinstance(energy, int): # If alpha is an int
            # find indicies of the FEDU_Alpha closest to user alpha input.
            idx = np.nanargmin(np.abs(FEDU_energy - energy))
            self.energyIdx = np.array([idx])
        else: # If alpha is an array.
            self.energyIdx = -1*np.ones(len(energy), dtype = np.int)
            # Iterate, and do the same as the other case in the if statement.
            for i in range(len(self.energyIdx)):
                self.energyIdx[i] = np.nanargmin(np.abs(FEDU_energy - energy[i]))
        return self.energyIdx
    
    def loadMagEphem(self, Bmodel = 'TS04D'):
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
        RETURNS: self.magEphem dictionary
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
        self.magEphem['Loss_Cone_Alpha_n'] = self.magEphem['Loss_Cone_Alpha_n'][magInd]
        self.magEphem['Loss_Cone_Alpha_s'] = self.magEphem['Loss_Cone_Alpha_s'][magInd]
        self.magEphem['BoverBeq'] = self.magEphem['BoverBeq'][magInd]
        return

if __name__ == '__main__':
    rb_id = 'A'
    #alpha = np.array([5, 90, 180])
    date = datetime(2017, 3, 31)
    #alpha = [10, 90, 170]
    #energy = [35, 50, 100]
    fluxObj = magEISspectra(rb_id, date, dataLevel = 3)
#    fluxObj.tBounds = [datetime(2017, 3, 31, 11, 15), datetime(2017, 3, 31, 11, 25)]
    fluxObj.tBounds = [datetime(2017, 3, 31, 11, 15), datetime(2017, 3, 31, 11, 20)]
    fluxObj.loadMagEIS(instrument = 'LOW', highrate = True)
    fluxObj.plotHighRateTimeSeries(smooth = 10)
    #fluxObj.plotHighRateSpectra(E_ch = 1, scatterS = 50)

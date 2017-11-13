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

class PlotHighrate:
    def __init__(self, sc_id, date, **kwargs):
        """
        NAME:    PlotHighrate(self, sc_id, date, **kwargs)
        USE:     This class plots the MagEIS highrate data. 
        INPUT:   REQUIRED:
                    sc_id: Either 'A' or 'B' (case doesnt matter).abs
                    date: A datetime object for title purposes.
                 OPTIONAL:
                    est_spin = 10.9: Estimated spin period
                    n_sectors = None: Number of highrate pitch angle 
                        sectors (1000 for RBSP-A, 64 for RBSP-B).
                    spin_thresh = 0.1: The allowed threshold to
                        determine the true spin rate of RBSP.
                    alphaKey = HighRate_Alpha360: Which pitch angle
                        data product to plot in the pitch angle
                        scatter plot function.
                        HighRate_Alpha is pitch angles 0-180 deg,
                        and HighRate_Alpha360 is from 0-360 deg.
        AUTHOR:  Mykhaylo Shumko
        RETURNS: None
        MOD:     2017-11-10
        """
        self.sc_id = sc_id
        # Modify these parameters if working with unusual data.
        self.est_spin = kwargs.get('est_spin', 10.9) # Spin period in seconds
        self.n_sectors = kwargs.get('n_sectors', None) # Number of highrate sectors
        if self.n_sectors is None:
            if self.sc_id.lower() == 'a':
                self.n_sectors = 1000
            else:
                self.n_sectors = 64
        # Threshold for self.est_spin to resolve the spin time.
        self.spin_thresh = kwargs.get('spin_thresh', 0.1)
        # Set data keys for pitch angle
        self.alphaKey = kwargs.get('alphaKey', None)
        if self.alphaKey is None:
            if self.sc_id.lower() == 'a':
                self.alphaKey = 'HighRate_Alpha'
            else:
                self.alphaKey = 'FEDU_Unbinned_Alpha'

        self.fluxKey = kwargs.get('fluxKey', None)
        if self.fluxKey is None:
            if self.sc_id.lower() == 'a':
                self.fluxKey = 'HighRate'
            else:
                self.fluxKey = 'FEDU_Unbinned_0to360'
        return None

    def _resolveSpinTimes(self, deflatTime):
        """
        NAME:    _resolveSpinTimes(self, deflatTime)
        USE:     Resolves the data time stamps (one per spin)
                 to (one per pitch angle). Effectively this 
                 fills every spin period with self.n_sectors
                 time stamps.
        INPUT:   REQUIRED:
                    deflatTime: If True, will fill in the time
                    stamps, if false, will expand the time array
                    into 2D, so each spin can be easily plotted.
                    This seems useless for now...
                 OPTIONAL:
                    None
        AUTHOR:  Mykhaylo Shumko
        RETURNS: self.times (a 2D time array for plotting)
        MOD:     2017-11-10
        """ 
        # If self.times array is not created, make it
        # if not hasattr(self, "times"):
        self.times = np.repeat(self.magEISdata['Epoch'][:, np.newaxis],
            self.magEISdata[self.alphaKey].shape[1] ,axis=1)
        if deflatTime:         
            # Find valid Indicies. Minus one is to calculate dTspin.
            for i in range(self.magEISdata[self.alphaKey].shape[0]-1):
                # Logic that checks for data gaps. If dTSpin is outside 
                # est_spin +/-(1+spin_thresh), assume est_spin spin period
                dTspin = (self.times[i+1, 0] - self.times[i, 0]).total_seconds()
                if ((dTspin > self.est_spin*(1+self.spin_thresh)) & 
                    (dTspin < self.est_spin*(1-self.spin_thresh)) ):
                    dTspin = self.est_spin
                #print(self.sc_id, self.n_sectors)
                self.times[i, :self.n_sectors] = np.array([self.times[i, iAlpha] + 
                    timedelta(seconds = iAlpha*dTspin/self.n_sectors) for iAlpha in 
                    range(self.n_sectors)])
            # This just interpolates the last spin assuming the previous spin
            self.times[-1, :self.n_sectors] = np.array([self.times[-1, iAlpha] + 
                    timedelta(seconds = iAlpha*dTspin/self.n_sectors) for iAlpha in 
                    range(self.n_sectors)]) 
                 
        # Flatten times for time series plotting.
        self.times = self.times[:, :self.n_sectors].flatten()
        return self.times

    def getFluxTimeseries(self, fluxConvert=True, smooth=1):
        """
        NAME:    getFluxTimeseries(self, fluxConvert=True, smooth=1)
        USE:     This function calculates the flux or count timeseries.
        INPUT:   REQUIRED:
                    None
                 OPTIONAL:
                    fluxConvert=True: Convert to units of flux using
                        Seth's G0dE factors.
                    smooth = 1: Data points for running average (concolution
                        kernel size. If 1, do not convolve) 
        AUTHOR:  Mykhaylo Shumko
        RETURNS: The self.flux (or count) timeseries array of (nTime, nEnergy)
        MOD:     2017-11-10
        """
        nT, nS, nE = self.magEISdata[self.fluxKey].shape
        self.flux = np.nan*np.ones((nT*self.n_sectors, nE), dtype=float)
        if self.instrument.lower() == 'low':
            nE = 7 # RBSP-B LOW detector has 7 working energy channels

        self._resolveSpinTimes(True)

        for ee in range(nE): 
            self.flux[:, ee] = self.magEISdata[self.fluxKey][:, :self.n_sectors, ee].flatten()
            if fluxConvert:
                self.flux[:, ee] /= self.G0dE[ee]
            if smooth > 1:
                self.flux[:, ee] = np.convolve(self.flux[:, ee], np.ones(smooth)/smooth, 
                    mode='same')
        return self.times, self.flux

    def plotTimeseries(self, **kwargs):
        """
        NAME:    plotTimeseries(self, **kwargs)
        USE:     Plots the MagEIS highrate timeseries.
        INPUT:   REQUIRED:
                    None
                 OPTIONAL:
                    ax = None: Subplot object to plot on. Will create a bx class attribute
                    E_ch = None: An array of channel numbers, or an integer channel to plot.
                    deflatTime = True: Calculate the time stamps for each pitch angle
                        sector. If False, this is pretty useless.
                    smooth = 1: Running average size in data points. If 1, no running average
                        is calculated.
                    chLegend = True: Plot channel legend
                    pltFlux = True: Plot flux or counts
                    pltLabels = True: Add x and y labels and add basic title to plot.
        AUTHOR:  Mykhaylo Shumko
        RETURNS: self.bx (subplot object)
        MOD:     2017-11-10
        """
        # Plotting parameters
        ax = kwargs.get('ax', None)
        E_ch = kwargs.get('E_ch', None)
        deflatTime = kwargs.get('deflatTime', True)
        smooth = kwargs.get('smooth', 1)
        chLegend = kwargs.get('chLegend', True)     
        pltFlux = kwargs.get('pltFlux', True) 
        pltLabels = kwargs.get('pltLabels', True)

        if ax is None:
            fig = plt.figure(figsize=(15, 10), dpi=80, facecolor = 'white')
            gs = gridspec.GridSpec(1,1)
            self.bx = fig.add_subplot(gs[0, 0], facecolor='white')
        else:
            self.bx = ax
    
        # Deflatten the 1d time array into 2d where each pitch angle is given 
        # a unique time stamp to where the spacecraft was pointing at the time.
        self._resolveSpinTimes(deflatTime)
        # Get flux (or counts)
        self.getFluxTimeseries(fluxConvert=pltFlux, smooth=smooth) 
         # Get energy channels to plot.
        if E_ch is None:
            E_ch = np.arange(7)
        elif not isinstance(E_ch, (list, np.ndarray)):
            E_ch = [E_ch]
            
        for ee in E_ch:
            validF = np.where(self.flux[:, ee] > 100)[0]
            self.bx.plot(self.times[validF], self.flux[validF, ee], 
                label='{}-{} keV'.format(self.Elow[ee],
                self.Ehigh[ee]))
        # Label and beautify the plot.
        self.bx.set(yscale='log')
        if chLegend: # Add legend
            self.fluxLegend = self.bx.legend()
        if pltLabels: # Add plot labels
            if pltFlux:
                fluxLabel = r'Electron flux $(cm^2 \ sr \ s \ keV)^{-1}$'
            else:
                fluxLabel = r'Electron counts/s'
            self.bx.set(xlabel='UTC', ylabel=fluxLabel, title='MagEIS-{} Highrate from {}'.format(
                self.sc_id.upper(), self.date.date()))
        if ax is None: # If ax was not passed, show plot
            plt.show()    
        return self.bx

    def plotAlpha(self, **kwargs):
        """
        NAME:    plotAlpha(self, **kwargs)
        USE:     This function plots the MagEIS highrate flux or counts.
        INPUT:   REQUIRED:
                    None
                 OPTIONAL:
                    'ax'= None: subplot object to plot on. If not specified, 
                        will create own figure.
                    'E_ch'= 0: Energy channel to plot
                    'scatterS' = 10: Size of scatter symbols to plot.
                    'cmin' = None: Min flux value range for colorbar
                    'cmax' = None: Max flux value range for colorbar
                    smooth = 1: Size in data points to do a running average over.
                    pltLabels = True: Wether or not to label ans title the axis.
                    downsampleAlpha = 1: Plot the downsampleAlpha'th pitch 
                        angle for quicker plotting on slower machines.
                    deflatTime = True: Refromat the 'epoch' array so that each
                        pitch angle is attributed to a unique time stamp. Saved in 
                        self.times
                    plotCb = True: Switch to plot colorbar.
                    cax = None: Subplot to plot the colorbar on.
                    aspect = None: Colorbar aspect rario
        AUTHOR:  Mykhaylo Shumko
        RETURNS: self.ax
        MOD:     2017-11-10
        """
        ax = kwargs.get('ax', None)
        E_ch = kwargs.get('E_ch', 1)
        scatterS = kwargs.get('scatterS', 10)
        cmin = kwargs.get('cmin', None)
        cmax = kwargs.get('cmax', None)
        smooth = kwargs.get('smooth', 1)
        pltLabels = kwargs.get('pltLabels', True)
        pltFlux = kwargs.get('pltFlux', True)
        downsampleAlpha = kwargs.get('downsampleAlpha', 1)
        deflatTime = kwargs.get('deflatTime', True) # alphaSpinTimes
        cax = kwargs.get('cax', None)
        aspect = kwargs.get('cAspect', None)
        plotCb = kwargs.get('plotCb', True)
        
        if ax is None:
            fig = plt.figure(figsize=(15, 10), dpi=80, facecolor = 'white')
            gs = gridspec.GridSpec(1,1)
            self.ax = fig.add_subplot(gs[0, 0], facecolor='black')
        else:
            self.ax = ax

        # Deflatten the 1d time array into 2d where each pitch angle is given 
        # a unique time stamp to where the spacecraft was pointing at the time.
        self._resolveSpinTimes(deflatTime)
        # Get flux (or counts)
        self.getFluxTimeseries(fluxConvert=pltFlux, smooth=smooth) 
        # Get flattened pitch angles
        alphas = self.magEISdata[self.alphaKey][:, :self.n_sectors].flatten()
        # Do the scatter plot.       
        self.sc = self.ax.scatter(self.times[::downsampleAlpha], 
            alphas[::downsampleAlpha], 
            c=self.flux[::downsampleAlpha, E_ch], norm=matplotlib.colors.LogNorm(),
            cmap=plt.get_cmap('rainbow'), 
            vmin=cmin, vmax=cmax, s=scatterS)
        if plotCb:
            if aspect is None:
                try:
                    self.cb = plt.colorbar(self.sc, ax=self.ax, cax=cax,
                        orientation='vertical', label=r'$(keV \ cm^2 \ sr \ s)^{-1}$')
                except ValueError:
                    print('****No valid data found.****')
                    raise
            else:
                try:
                    self.cb = plt.colorbar(self.sc, ax=self.ax, cax=cax,
                        orientation='vertical', aspect=aspect, 
                            label=r'$(keV \ cm^2 \ sr \ s)^{-1}$')
                except ValueError:
                    print('****No valid data found.****')
                    raise

        if ax is None:
            self.ax.set_xlim(self.times[0], 
                self.times[-1])
        if self.alphaKey == 'HighRate_Alpha360':
            self.ax.set_ylim(0, 360)
        else:
            self.ax.set_ylim(0, 180)
        if pltLabels:
            self.ax.set(ylabel='Pitch angle (deg)', xlabel='UTC', 
                title='RBSP-{} magEIS HighRate from {}'.format(
                self.sc_id.upper(), self.date.date()))
        
        # Set pitch angle ticks
        alpha_tick_locator = matplotlib.ticker.FixedLocator(range(0, 360, 10))
        self.ax.yaxis.set_minor_locator(alpha_tick_locator)
        self.ax.tick_params(which='minor')
        
        if ax is None:
            plt.show()
        return self.ax   

class PlotRel03:
    def __init__(self, sc_id, date, **kwargs):
        """
        NAME:    PlotRel03(self, sc_id, date, **kwargs)
        USE:     This class plots the MagEIS rel03 data. 
        INPUT:   REQUIRED:
                    sc_id: Either 'A' or 'B' (case doesnt matter).abs
                    date: A datetime object for title purposes.
                 OPTIONAL:
                    fluxKey = 'FEDU': Electron flux key to plot.
        AUTHOR:  Mykhaylo Shumko
        RETURNS: None
        MOD:     2017-11-11
        """
        # Get basic key parameters
        self.fluxKey = kwargs.get('fluxKey', 'FEDU')
        # Calculate the energy channel values for the day.
        self._calcDailyElectronEnergies() 
        return None

    def plotEpinAvgSpectra(self, ax=None, **kwargs):
        """
        NAME:    plotEpinAvgSpectra(self, ax=None, **kwargs)
        USE:     Plots the spin averaged electron flux spectra from magEIS.
        INPUT:   REQUIRED:
                    None
                 OPTIONAL:
                    ax = None: A subplot object. If none specied, it will 
                        create one.
                    yscale = 'log': y scale of the plot.
                    pltLabel = True: Plot labels and titles
                    pltColorbar = True: Plot flux colorbar
                    maxValidE = 21: Maximum valid energy channel
                    cax = None: Colorbar axis
        AUTHOR:  Mykhaylo Shumko
        RETURNS: Subplot object, self.ax
        MOD:     2017-11-11
        """
        yscale = kwargs.get('yscale', 'log')
        pltLabel = kwargs.get('pltXlabel', True)
        pltColorbar = kwargs.get('pltColorbar', True)
        maxValidE = kwargs.get('maxValidE', 21) # Max valid energy channel
        cax = kwargs.get('cax', None)

        if ax is None:
            fig = plt.figure(figsize=(15, 10), dpi=80, facecolor = 'white')
            gs = gridspec.GridSpec(1,1)
            self.ax = fig.add_subplot(gs[0, 0], facecolor='w')
        else:
            self.ax = ax
    
        # Calculate the mean flux as a function of time and energy
        self.magEISdata[self.fluxKey][self.magEISdata[self.fluxKey] == -1e31] = np.nan
        if 'FESA' not in self.magEISdata.keys(): # Create spin-avraged flux if there is none.
            self.magEISdata['FESA'] = np.nanmean(self.magEISdata['FEDU'], axis=1)

        # Plot spectra
        wTT, wEE = np.meshgrid(range(len(self.magEISdata['Epoch'])),
            self.magEISdata['eEnergy'][:maxValidE])
        # Fill wTT with datetimes instead of indicies
        wTT = np.broadcast_to(self.magEISdata['Epoch'], wEE.shape)
        cs = self.ax.pcolormesh(wTT, wEE, self.magEISdata['FESA'].T[:maxValidE, :], 
            cmap = plt.get_cmap('gnuplot2'), norm=colors.LogNorm(),
            vmin=np.nanmin(self.magEISdata['FESA'].T[:maxValidE, :]), 
            vmax=np.nanmax(self.magEISdata['FESA'].T[:maxValidE, :]))
        
        # Beautify and label plot
        self.ax.set_yscale(yscale)
        if pltLabel:
            self.ax.set(xlabel='UTC', ylabel='Energy (keV)', 
                title='MagEIS-{} energy spectra {}'.format(
                self.sc_id.upper(), self.date.date()))
        if pltColorbar:
            plt.colorbar(cs, label=r'Electron flux $(cm^2 \ sr \ s \ keV)^{-1}$', 
                ax=self.ax, cax=cax, ticks=matplotlib.ticker.LogLocator())
        if ax is None:
            plt.show()
        return self.ax

    def plotUnidirectionalFlux(self, alpha, energy=None, ax=None, **kwargs):
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
        MOD:     2017-11-10
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
        self.magEISdata[self.fluxKey][self.magEISdata[self.fluxKey] == -1e31] = np.nan
        
        # Plot all of the energies and pitch angles.
        for e in self.energyIdx:
            for a in self.alphaIdx:
                self.ax.plot(self.magEISdata['Epoch'][:], 
                    self.magEISdata[self.fluxKey][:, a, e],
                    label = r'E = {} keV, $\alpha$ = {}'.format(
                    int(self.magEISdata['eEnergy'][e]), 
                    int(self.magEISdata['FEDU_Alpha'][a])) )
        if pltLegendLoc:
            self.fluxLegend = self.ax.legend(loc = pltLegendLoc)    
        if pltXlabel:
            self.ax.set_xlabel('UTC')
        self.ax.set_ylabel(r'Electron flux $(cm^2 \ sr \ s \ keV)^{-1}$')  
        self.ax.set_yscale(yscale)   
        
        if ax is None:
            plt.show()
        return self.ax

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

    def _calcDailyElectronEnergies(self, method='mode'):
        """
        NAME:    _calcDayEnergies(self, method='mode')
        USE:     This function will add a new key, the energy channels for
                 electrons
        INPUT:   REQUIRED:
                    None
                 OPTIONAL:
                    method = 'mode': How to calculate the energies, 
                        Can use mean, median, or mode.
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

class PlotMageis(PlotHighrate, PlotRel03):
    def __init__(self, sc_id, date, dtype, **kwargs):
        """

        """
        # Initialize basic variables and paths
        self.date = date
        self.sc_id = sc_id
        self.dtype = dtype
        self.mageisDir = kwargs.get('mageisDir', directories.mageis_dir(sc_id))
        self.magEphemDir = kwargs.get('magEphemDir', 
            directories.rbsp_magephem_dir(sc_id))
        self.tRange = kwargs.get('tRange', None) 
        self.dataLevel = kwargs.get('dataLevel', 3)

        assert dtype in ['highrate', 'rel03'], 'dtype must be highrate or rel03'
        self.instrument = kwargs.get('instrument', None)
        if (dtype == 'highrate') and (self.instrument is None):
            raise TypeError('Need to specify instrument for highrate data. '
                '(LOW, M35, M75, HIGH)')
        self._loadmageis() # Load MagEIS data

        if dtype == 'highrate': # Initialize parent highrate class
            PlotHighrate.__init__(self, sc_id, date, **kwargs)
            self._getDetectorParams() # Get instrument parameters
        else: # Initialize parent rel03 class
            PlotRel03.__init__(self, sc_id, date, **kwargs)
        return

    def loadMagEphem(self, Bmodel='TS04D'):
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
        splitDate = self.date.date().isoformat().split('-')
        joinedDate = ''.join(splitDate)
        magEphemName = sorted(glob.glob(os.path.join(self.magEphemDir,
            ('rbsp{}*MagEphem*{}*{}*.txt'.format(self.sc_id.lower(), 
            Bmodel, joinedDate)))))
        assert len(magEphemName) == 1, ('Error, none or multiple'
        ' magnetic ephemeris files found!')
        magEphemName = magEphemName[0]
        self.magEphem = spacepy.datamodel.readJSONheadedASCII(magEphemName)
        
        # Convert times
        self.magEphem['DateTime'] = np.array(list(map(
            lambda x: dateutil.parser.parse(x).replace(tzinfo=None),
            self.magEphem['DateTime'])))

        # Filter data by times
        if self.tRange is not None:
            magInd = np.where((self.magEphem['DateTime'] >= self.tRange[0]) &
                (self.magEphem['DateTime'] <= self.tRange[1]))[0]
            assert len(magInd) > 0, ('ERROR: no filtered magEphem found'
                ' in the time range specified! Check tBounds keyword')
        else:
            magInd = range(len(self.magEphem['DateTime']))

        self.magEphem['DateTime'] = self.magEphem['DateTime'][magInd]
        self.magEphem['Bmin'] = self.magEphem['Bmin_gsm'][magInd, 3]
        self.magEphem['Lstar'] = self.magEphem['Lstar'][magInd, :]
        self.magEphem['EDMAG_MLT'] = self.magEphem['EDMAG_MLT'][magInd]
        self.magEphem['Loss_Cone_Alpha_n'] = self.magEphem['Loss_Cone'
            '_Alpha_n'][magInd]
        self.magEphem['Loss_Cone_Alpha_s'] = self.magEphem['Loss_Cone'
            '_Alpha_s'][magInd]
        self.magEphem['BoverBeq'] = self.magEphem['BoverBeq'][magInd]
        self.magEphem['EDMAG_MLAT'] = self.magEphem['EDMAG_MLAT'][magInd]
        return

    def _getDetectorParams(self):
        # Define magEIS detector constants
        
        if self.instrument.lower() == 'low':
            if self.sc_id.lower() == 'a':
                self.Emid = [34, 54, 78, 108, 143, 182, 223] # keV
                self.Elow = [29, 46, 68, 95, 126, 164, 206] # keV
                self.Ehigh = [41, 66, 92, 126, 164, 204, 247] # keV
                # Units of (keV cm^2 sr)
                self.G0dE = [4.13E-2, 5.73E-2, 6.056E-2, 6.88E-2, 7.35E-2, 
                    6.90E-2, 5.98E-2]
                self.Ebins = [29, 41, 66, 92, 126, 164, 204, 247]
                
            if self.sc_id.lower() == 'b':
                self.Emid = [32, 51, 74, 101, 132, 168, 208] # keV
                self.Elow = [27, 43, 63, 88, 117, 152, 193] # keV
                self.Ehigh = [39, 63, 88, 117, 150, 188, None] # keV 
                # Units of (keV cm^2 sr)
                self.G0dE = [4.33E-2, 5.41E-2, 5.926E-2, 6.605E-2, 6.460E-2,
                    6.23E-2, 5.96E-2]
                self.Ebins = [27, 39, 63, 88, 117, 150, 188]
        return


    def _loadmageis(self):
        # Find the file.
        splitDate = self.date.date().isoformat().split('-')
        joinedDate = ''.join(splitDate)
        if self.dtype == 'highrate':
            if self.sc_id.lower() == 'a':
                searchStr = os.path.join(self.mageisDir, 
                    'rbsp{}_int_ect-mageis{}-hr-L{}_{}*.cdf'.format(self.sc_id.lower(), 
                    self.instrument.upper(), self.dataLevel, joinedDate))
            else:
                searchStr = os.path.join(self.mageisDir, 
                    'rbsp{}_int_ect-mageis{}-L{}_{}*.cdf'.format(self.sc_id.lower(), 
                    self.instrument.upper(), self.dataLevel, joinedDate))
        if self.dtype == 'rel03':
            searchStr = os.path.join(self.mageisDir, 
                'rbsp{}_rel03_ect-mageis-L{}_{}*.cdf'.format(self.sc_id.lower(), 
                self.dataLevel, joinedDate))

        magEISname = sorted(glob.glob(searchStr))
        assert len(magEISname) == 1, ('Error, none or multiple magEIS files '
            'found! \n Search str: {} \n {}'.format(searchStr, magEISname))
        magEISname = magEISname[0]
        # Load and copy the data so that it can be modified.
        self.magEISdata = pycdf.CDF(magEISname).copy() 

        if self.tRange is not None:
            self._filtermageis()

    def _filtermageis(self):
        # Now find filtered data indicies if tBounds is specified.
        magEISInd = np.where((self.magEISdata['Epoch'][:] >= self.tRange[0]) &
            (self.magEISdata['Epoch'][:] <= self.tRange[1]))[0]
        assert len(magEISInd) > 0, ('ERROR: no filtered spectra found in '
            'the time range specified! Check tBounds keyword')
        N = len(self.magEISdata['Epoch'])
        # Filter the data.
        for key in self.magEISdata.keys():
            # If key is a function of time, filter it.
            if len(self.magEISdata[key]) == N:
                self.magEISdata[key] = self.magEISdata[key][magEISInd]

if __name__ == '__main__':
    rb_id = 'A'
    date = datetime(2017, 3, 31)
    tRange = [datetime(2017, 3, 31, 11, 15), datetime(2017, 3, 31, 11, 25)]

    # Test rel03 data plotting
    dtype = 'rel03'
    rel03Obj = PlotMageis(rb_id, date, dtype, tRange=tRange)
    rel03Obj.plotUnidirectionalFlux(90)
    rel03Obj.plotEpinAvgSpectra()

    # Test highrate data plotting
    instrument = 'LOW'
    dtype = 'highrate'
    hrObj = PlotMageis(rb_id, date, dtype, tRange=tRange, instrument=instrument,
        alphaKey='HighRate_Alpha')
    hrObj.plotAlpha()
    hrObj.plotTimeseries()

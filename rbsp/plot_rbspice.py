# This class contains the scripts to plot RBSPICE data.
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
from datetime import datetime

from spacepy import pycdf

# Try to import the directorories module
try:
    sys.path.insert(0, '/home/mike/Dropbox/0_grad_work')
    import directories
except SystemError:
    print('Could not import the directories.py file, '
        'please supply data directories manualy!')

class plot_rbspice:
    def __init__(self, sc_id, date, **kwargs):
        """
        NAME:    plot_rbspice(sc_id, date, **kwargs)
        USE:     This class will load and plot RBSPICE data.
        INPUT:   REQUIRED:
                    sc_id either 'A' or 'B', date is a 
                    datetime object of a date to plot. 
                 OPTIONAL:
                    ### FINISH WRITING THIS LATER ###

        AUTHOR:  Mykhaylo Shumko
        RETURNS: None
        MOD:     2017-10-17
        """
        self.sc_id = sc_id
        self.date = date
        splitDate = date.date().isoformat().split('-')
        self.joinedDate = ''.join(splitDate)

        self.dataLevel = kwargs.get('dataLevel', 3)
        self.rbspiceDir = kwargs.get('rbspiceDir',
            directories.rbspice_dir(sc_id))
        self.magEphemDir = kwargs.get('magEphemDir', 
            directories.rbsp_magephem_dir(sc_id))
        self.tBounds = kwargs.get('tBounds', None) 
        return
    
    def loadData(self, **kwargs):
        """
        NAME:    loadData(self, **kwargs)
        USE:     This function will load and filter RBSPICE
                 data by time.
        INPUT:   REQUIRED:
                    None
                 OPTIONAL:
                    ### FINISH WRITING THIS LATER ###
        AUTHOR:  Mykhaylo Shumko
        RETURNS: None
        MOD:     2017-10-17
        """
        # Optional variables.
        dataType = kwargs.get('dataType', 'ESRHELT').upper()
        
        # Find and load the data
        searchStr = os.path.join(self.rbspiceDir, 
            'rbsp-{}-rbspice_lev-{}_{}_{}_*.cdf'.format(self.sc_id.lower(),
             self.dataLevel, dataType, self.joinedDate))
        # rbsp-a-rbspice_lev-3_ESRHELT_20170331_v1.1.9-01.cdf
        rbspiceName = sorted(glob.glob(searchStr))
        assert len(rbspiceName) == 1, ('Error, none or multiple RBSPICE files '
            'found! \n Search str: {} \n {}'.format(searchStr, rbspiceName))
        rbspiceName = rbspiceName[0]
        self.rbspicedata = pycdf.CDF(rbspiceName).copy() 
        
        # Now find filtered data indicies if tBounds is specified.
        if self.tBounds is not None:
            idt = np.where((self.rbspicedata['Epoch'][:] >= self.tBounds[0]) &
                (self.rbspicedata['Epoch'][:] <= self.tBounds[1]))[0]
            assert len(idt) > 0, ('ERROR: no filtered data found in '
                'the time range specified! Check tBounds keyword')
        else:
            idt = range(len(self.rbspicedata['Epoch']))
            
        # Now filter the data that will be used for plotting
        self.rbspicedata['Epoch'] = self.rbspicedata['Epoch'][idt]
        self.rbspicedata['Alpha'] = self.rbspicedata['Alpha'][idt, :]
        self.rbspicedata['FEDU_Alpha'] = self.rbspicedata['FEDU_Alpha'][idt, :] 
        self.rbspicedata['FEDU_AlphaRange'] = (
            self.rbspicedata['FEDU_AlphaRange'][idt, :, :] )
        # FEDU units: Counts/(MeV-cm^2-s-sr)
        self.rbspicedata['FEDU'] = self.rbspicedata['FEDU'][idt, :, :]  
        return
        
    def plotSpectra(self, **kwargs):
        """
        NAME:    plotSpectra(self, **kwargs)
        USE:     This function will plot the RBSPICE spectra for a particular 
                 telescope (energy and time on y and x axes, and color
                 is the FEDU amplitude).
        INPUT:   REQUIRED:
                    None
                 OPTIONAL:
                    tel = 0: telescope to plot
                    ax = None: Subplot object to plot on
                    cmin/cmax = None: colorbar min and max, 
                        will find automatically if none
                    cax = None: Colorbax axis
                    logE = True: Log the energy axis     
        AUTHOR:  Mykhaylo Shumko
        RETURNS: None
        MOD:     2017-10-17
        """
        tel = kwargs.get('tel', 0) 
        ax = kwargs.get('ax', None)
        cmin = kwargs.get('cmin', None)
        cmax = kwargs.get('cmax', None)
        cax = kwargs.get('cax', None)
        logE = kwargs.get('logE', True)
        
        if ax is None:
            fig = plt.figure(figsize=(15, 10), dpi=80, facecolor = 'white')
            gs = gridspec.GridSpec(1,1)
            self.ax = fig.add_subplot(gs[0, 0], facecolor='k')
        else:
            self.ax = ax
            
        Ebins = 1000*np.array(self.rbspicedata['FEDU_EnergyRange'])[tel, 0, :]
        Ebins = np.append(Ebins, 1000*np.array(
            self.rbspicedata['FEDU_EnergyRange'])[tel, 1, -1])
        # Create a meshgrid, and populate it with energy and time
        sTT, sEE = np.meshgrid(range(len(self.rbspicedata['Epoch'])), Ebins)
        # Now replace the time indicies by time stamps
        sTT = np.broadcast_to(self.rbspicedata['Epoch'], sTT.shape)
        flux = self.rbspicedata['FEDU'][:, tel, :] 
        # Calculate the color bounds if not given 
        validFlux = np.where(flux > 0)
        if cmin is None:
            # Find the minimum flux greater than 0.
            cmin = np.min(flux[validFlux])
        if cmax is None:
            # Find the maximum flux.
            cmax = np.max(flux[validFlux])
        
        # Now plot the spectra
        cs = self.ax.pcolormesh(sTT, sEE, 
            np.transpose(flux)/1000, 
            cmap = plt.get_cmap('gnuplot2'), norm=colors.LogNorm(), 
            vmin=cmin, vmax=cmax)
        cb = plt.colorbar(cs, ax=self.ax, cax=cax, 
            label=r'Flux $(keV \ cm^2 \ s \ sr)^-1$') 
        self.ax.set(xlabel='UTC', ylabel = 'Energy (keV) \n tel{}'.format(tel))
        if logE:
            self.ax.set_yscale('log')  
        if ax is None:
            gs.tight_layout(fig)
            plt.show()
        return self.ax
        
    def plotTelecopeAlphaScatter(self, Ech, **kwargs):
        """
        NAME:    plotTelecopeAlphaScatter(self, **kwargs)
        USE:     This function will plot the RBSPICE spectra for a particular 
                 telescope (energy and time on y and x axes, and color
                 is the FEDU amplitude).
        INPUT:   REQUIRED:
                    Ech: Energy channel to plot
                 OPTIONAL:
                    ax = None: Subplot object to plot on
                    tel = 0: telescope to plot
        AUTHOR:  Mykhaylo Shumko
        RETURNS: None
        MOD:     2017-10-17
        """
        ax = kwargs.get('ax', None)
        cmin = kwargs.get('cmin', None)
        cmax = kwargs.get('cmax', None)
        cax = kwargs.get('cax', None)
        logE = kwargs.get('logE', True)
        Elabel = kwargs.get('Elabel', True)
                
        flux = self.rbspicedata['FEDU'][:, :, Ech]/1000
        
        # If Ech is an array, average over the energy channels.
        if hasattr(Ech, '__len__'):
            flux = np.sum(flux, axis=2)/len(Ech)
        # Calculate the color bounds if not given 
        validFlux = np.where(flux > 0)
        if cmin is None:
            # Find the minimum flux greater than 0.
            cmin = np.min(flux[validFlux])
        if cmax is None:
            # Find the maximum flux.
            cmax = np.max(flux[validFlux])
        
        if ax is None:
            fig = plt.figure(figsize=(15, 10), dpi=80, facecolor = 'white')
            gs = gridspec.GridSpec(1,1)
            self.ax = fig.add_subplot(gs[0, 0], facecolor='k')
        else:
            self.ax = ax
            
        for tel in range(len(self.rbspicedata['Telescope'])):
            sc = self.ax.scatter(self.rbspicedata['Epoch'], 
                self.rbspicedata['Alpha'][:, tel], 
                c=flux[:, tel], 
                norm=colors.LogNorm(), vmin=cmin, vmax=cmax)
                
        cb = plt.colorbar(sc, ax=self.ax, cax=cax, 
            label=r'Flux $(keV \ cm^2 \ s \ sr)^-1$')
                
        self.ax.set_xlim(self.rbspicedata['Epoch'][0], 
            self.rbspicedata['Epoch'][-1])
        self.ax.set(ylabel=r'$\alpha_{sc}$', xlabel='UTC')
        
        # Write the energy channel in the plot
        if Elabel:
            if hasattr(Ech, '__len__'):
                labelText = '{}-{} keV'.format(
                    round(1000*self.rbspicedata['FEDU_EnergyRange'][0, 0, 
                    Ech[0]]), 
                    round(1000*self.rbspicedata['FEDU_EnergyRange'][0, 1, 
                    Ech[-1]]))
            else:
                labelText = '{}-{} keV'.format(
                    round(1000*self.rbspicedata['FEDU_EnergyRange'][0, 0, Ech]), 
                    round(1000*self.rbspicedata['FEDU_EnergyRange'][0, 1, Ech]))
            self.ax.text(1, 1, labelText, size=20, transform=self.ax.transAxes,
                ha="right", va="top",
                bbox=dict(boxstyle="square", facecolor='white'))

        if ax is None:
            gs.tight_layout(fig)
            plt.show()
        return self.ax
        
        
    def _makeAlphaTimeVertecies_(self, times, alphas):
        vertex = np.array(
                    [[pPerp[i, j], pParallel[i, j]], 
                    [pPerp[i+1, j], pParallel[i+1, j]], 
                    [pPerp[i+1, j+1], pParallel[i+1, j+1]], 
                    [pPerp[i, j+1], pParallel[i, j+1]]] )
        
if __name__ == '__main__':
    rb_id = 'A'
    tBounds = [datetime(2017, 3, 31, 11, 17, 8), 
        datetime(2017, 3, 31, 11, 17, 20)]
    pltObj = plot_rbspice(rb_id, tBounds[0], tBounds=tBounds)
    pltObj.loadData()
    
    fig = plt.figure(figsize=(15, 10), dpi=80, facecolor = 'white')
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0, 0], facecolor='k')
    
#    pltObj.plotSpectra(tel=5, ax=ax, logE=False, cmin=1E4, cmax=1E6)
#    ax.set_ylim(25, 100)
#    plt.show()

    pltObj.plotTelecopeAlphaScatter(np.arange(0, 20), ax=ax)
    plt.show()

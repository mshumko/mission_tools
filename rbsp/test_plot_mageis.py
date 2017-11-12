# This script does numberous unittests on the plot_mageis.py library to determine
# if it works.
import unittest
import matplotlib.pyplot as plt
from datetime import datetime

import plot_mageis

RB_ID = 'A'
DATE = datetime(2017, 3, 31)
T_BOUNDS = [datetime(2017, 3, 31, 11, 15), datetime(2017, 3, 31, 11, 20)]

class TestMagEISLib(unittest.TestCase):

    def test_highrate(self):
        fig, ax = plt.subplots()
        fluxObj = plot_mageis.magEISspectra(RB_ID, DATE, dataLevel=3)
        fluxObj.tBounds = T_BOUNDS
        fluxObj.loadMagEIS(instrument='LOW', highrate=True)
        fluxObj.plotHighRateTimeSeries(ax=ax, smooth=10)

    # def test_int(self):
    #     fig, ax = plt.subplots()
    #     fluxObj = plot_mageis.magEISspectra(RB_ID, DATE, dataLevel = 3)
    #     fluxObj.tBounds = T_BOUNDS
    #     fluxObj.loadMagEIS(highrate=False, relType='int')
    #     fluxObj.plotHighRateTimeSeries(ax=ax, smooth = 10, n_sectors=11)
        
    
    def test_rel03(self):
        fig, ax = plt.subplots()
        fluxObj = plot_mageis.magEISspectra(RB_ID, DATE, dataLevel = 3)
        fluxObj.tBounds = T_BOUNDS
        fluxObj.loadMagEIS(highrate=False, relType='rel03')
        fluxObj.plotHighRateTimeSeries(ax=ax, smooth = 10, n_sectors=11)
        #fluxObj.plotHighRateSpectra(E_ch = 1, scatterS = 50)
        
        

if __name__ == '__main__':
    unittest.main()

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
        for rb_id in ['A', 'b']:
            fig, ax = plt.subplots()
            hrObj = plot_mageis.PlotMageis(rb_id, DATE, 'highrate', tRange=T_BOUNDS, 
                instrument='low')
            hrObj.plotAlpha(ax=ax)
            hrObj.plotTimeseries(ax=ax)

    #     rb_id = 'A'
    # date = datetime(2017, 3, 31)
    # tRange = [datetime(2017, 3, 31, 11, 15), datetime(2017, 3, 31, 11, 25)]
    
    def test_rel03(self):
        for rb_id in ['A', 'b']:
            fig, ax = plt.subplots()
            rel03Obj = plot_mageis.PlotMageis(rb_id, DATE, 'rel03', tRange=T_BOUNDS)
            rel03Obj.plotUnidirectionalFlux(90, ax=ax)
            rel03Obj.plotEpinAvgSpectra(ax=ax)
        # fluxObj = plot_mageis.magEISspectra(RB_ID, DATE, dataLevel = 3)
        # fluxObj.tBounds = T_BOUNDS
        # fluxObj.loadMagEIS(highrate=False, relType='rel03')
        # fluxObj.plotHighRateTimeSeries(ax=ax, smooth = 10, n_sectors=11)
        # #fluxObj.plotHighRateSpectra(E_ch = 1, scatterS = 50)
        
        

if __name__ == '__main__':
    unittest.main()

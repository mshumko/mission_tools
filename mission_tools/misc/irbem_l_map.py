import IRBEM
import numpy as np
import matplotlib.pyplot as plt
import itertools
import datetime

### MODEL PARAMETERS ###
kp = 20
kext = 'T89'
time = datetime.datetime.now()
lons = np.linspace(-180, 180, 200)
lats = np.linspace(-90, 90, 200)
alt = 500

lon, lat = np.meshgrid(lons, lats)
lon, lat = lon.T, lat.T
L = np.zeros((len(lons), len(lats)), dtype=float)

### RUN IRBEM ###
magModel = IRBEM.MagFields(kext=kext)

for i,j in itertools.product(range(len(lons)), range(len(lats))):
    #print(i,j)
    X = {'x1':alt, 'x2':lat[i, j], 
         'x3':lon[i, j], 'dateTime':time}
    maginput = {'Kp':kp}
    magModel.make_lstar(X, maginput)
    L[i, j] = magModel.make_lstar_output['Lm'][0]
    
L = np.abs(L)
ix, iy = np.where(L == 1E31)
for i, j in zip(ix, iy):
    L[i, j] = np.nan
    
### SAVE DATA ###
np.save('irbem_l_lons', lon)
np.save('irbem_l_lats', lat)
np.save('irbem_l_l', L)

### PLOT DATA ###
levels = np.arange(2, 10, 2)
CS = plt.contour(lon, lat, L, levels=levels, colors='k')
plt.clabel(CS, inline=1, fontsize=10)
plt.show()

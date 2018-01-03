# This is the test code for the oblate epheroid calculations
import math
import copy
import numpy as np
import matplotlib.pyplot as plt

# WGS-72
f = 1/298.26 # Earth's flattening factor.
ee = 2*f - f**2

Re = 6378.135 # km

# The spherical results
X = (1917.5314168385216, 1549.8622520944375, -6367.881775563029)
XGEO = {'Alt': 450.4020353083497, 'Lat': -68.834230973252787, 'Lon': -99.344940100715803}

def _calcLat(X, thresh=0.001, maxIter=10):
    """
    This function iterates over the method described in CelesTrak's 
    Orbital Coordinate Systems, Part III algorithm to determine the geodetic 
    latitude and altitude for an oblate spheroid geometry.
    """
    # Calculate the latitude in spherical geometry first
    magX = np.sqrt(X[0]**2 + X[1]**2 + X[2]**2)        
    phiPrime = np.rad2deg(np.arcsin(X[2]/magX))
    phiPrime2 = copy.copy(phiPrime)

    phi = copy.copy(phiPrime)

    i=0
    while np.abs(phiPrime) > thresh and i < maxIter:
        i += 1
        phi = phiPrime
        C = np.sqrt(1 - ee*np.sin(np.deg2rad(phi))**2 )**(-1) 
        R = Re*np.cos(np.deg2rad(phi))
        z = Re*np.sin(np.deg2rad(phi))

        phiPrime = np.arctan2(z + Re*C*ee*np.sin(np.deg2rad(phi)), R)    
        phiPrime2 += phiPrime
        print(phiPrime, phiPrime2) 

        alt = R/np.cos(np.deg2rad(phiPrime2)) - Re*(C-1)
        print(R/np.cos(np.deg2rad(phiPrime)), alt, C)

    return

_calcLat(X)

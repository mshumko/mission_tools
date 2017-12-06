# This script makes a file that calculates the cross-spacecraft separation and 
# in-track lag.

import numpy as np
import csv

Re = 6371 # km.

def crowFliesDist(lat, lon, alt):
    """ 
    This function uses the haverside foruma to calculate the great-circle
    distance between two points.
    
    lat/lon have to be dimentions of 2 or 2*nT. alt dimentions should be 1 or nT
    """
    # Convert degrees to radians
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)
    dlat = np.subtract(lat[0], lat[1])
    dlon = np.subtract(lon[0], lon[1])
    # Calculate distance
    a = np.sin(dlat/2)**2 + np.cos(lat[0])*np.cos(lat[1])*np.sin(dlon/2)**2
    c = 2*np.arctan2(np.sqrt(a), np.sqrt(1-a))
    return (Re+alt)*c
    
# Testing
if __name__ == '__main__':
    lat = [50, 58]
    lon = [15, 3]
    print(crowFliesDist(lat, lon, 0))

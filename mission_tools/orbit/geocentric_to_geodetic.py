#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 08:27:51 2018
@author: dianamadera
"""
import math
#Geocentric to geodetic LLA
def geocentric_to_geodetic(geoc_lat,geoc_lon,geoc_h):
    #Tolerances:
    tol = 10^(-10) # in radians
    kmax = 100
    k = 1
    #Earth parameters:
    R_earth = 6378137 #in m
    eE = 0.0818 #considering the shape of Earth to be an oblate ellipsoid with eccentricity eE (from side):
    geod_lat_guess = math.radians(geoc_lat) #using geocentric lat as first guess for geodetic lat (in radians)
    geod_lon = geoc_lon #the geodetic longitude is same as geocentric longitude    
    h = geoc_h # h is height above sphere (geocentric height) in m.
    #Projections of vector r_ECEF  (ECEF = Earth Centered Earth Fixed coordinate system)
    #rX = (R_earth+h)*math.cos(math.radians(geoc_lat))*math.cos(math.radians(geoc_lon)) #Projection onto X
    #rY = (R_earth+h)*math.cos(math.radians(geoc_lat))*math.sin(math.radians(geoc_lon)) #Projection onto Y
    #rXY = math.sqrt(math.pow(rX,2) + math.pow(rY,2)) #Projection onto XY, which can be simplified to:
    rXY = (R_earth+h)*math.cos(math.radians(geoc_lat)) #Projection onto XY
    rZ = (R_earth+h)*math.sin(math.radians(geoc_lat)) #Projection onto Z
    geod_lat_old = geod_lat_guess
    N = R_earth/math.sqrt(1-(math.pow(eE*math.sin(geod_lat_old),2)))
    geod_lat = math.atan2(rZ + N*math.pow(eE,2)*math.sin(geod_lat_old),rXY)
    while abs(geod_lat - geod_lat_old) > tol and k < kmax:
        geod_lat_old = geod_lat
        N = R_earth/math.sqrt(1-(math.pow(eE*math.sin(geod_lat_old),2)))
        geod_lat = math.atan2(rZ + N*math.pow(eE,2)*math.sin(geod_lat_old),rXY)
        k = k+1   
    geod_h = rXY/math.cos(geod_lat) - N
    geod_lat = math.degrees(geod_lat)
    return geod_lat, geod_lon, geod_h

if __name__ == '__main__':
    #To call it:
    #[geod_lat, geod_lon, geod_h] = geocentric_to_geodetic(geoc_lat,geoc_lon,geoc_h)   
    #geoc_lat, geoc_lon have to be entered in degrees and geoc_h in meters
    #returns geod_lat, geod_lon in degrees and geod_h in meters 
    import timeit
    print(timeit.timeit(stmt='geocentric_to_geodetic(60, 0, 500000)', 
        number=50000, setup="from __main__ import geocentric_to_geodetic"))    
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 19:38:09 2017

@author: mike
"""
import dateutil.parser
import datetime

# ECI to lat, lon, alt Conversion
import numpy as np
Re = 6378.135 # km

def eci2gdz(X, dateTime, **kwargs):
    """
    Converts SGP4 ECI-TEME coordinates to GDZ cooridnates (lat, long, alt)
    of the spacecraft subpoint.
    
    info found here: https://celestrak.com/columns/v02n01/ and 
    http://celestrak.com/columns/v02n02/
    https://celestrak.com/publications/AIAA/2006-6753/faq.asp
    
    X  = (x, y, z) is the input vector in TEME coordinates.
    
    To vectorize the operation, feed in numpy arrays with data. Lists will not
    work!
    """
    sphere_shape = kwargs.get('sphere_shape', True)
    
    out = {}
    # If you give it an array of dateTimes
    assert type(dateTime) is not list, 'dateTime input cannot be a list! Use numpy arrays'
    if type(dateTime) is np.ndarray:
        gst = [Greenwich_sidereal_time(i) for i in dateTime]
    else:
        gst = Greenwich_sidereal_time(dateTime)
    if sphere_shape:
        magX = np.sqrt(X[0]**2 + X[1]**2 + X[2]**2)        
        out['Alt'] = magX - Re
        out['Lat'] = np.rad2deg(np.arcsin(X[2]/magX))
        
        # There is got to be a clever way to do a modulo, but this is the
        # quick fix.
        deg = np.rad2deg(np.arctan2(X[1],X[0]) - gst)
        while deg < -180:
            deg += 360
        while deg > 180:
            deg -= 360
        out['Lon'] = deg
                             
        #out['Lon'] = (deg % 180)
    else:
        # Code to calculate the subpoint for the oblate spheroid shape
        pass
    
    return out
    
def gdz2eci(lat, long, alt, dateTime, **kwargs):
    """
    Converts from GDZ coordinates to ECI
    
    To vectorize the operation, feed in numpy arrays with data. Lists will not
    work!
    """
    sphere_shape = kwargs.get('sphere_shape', True)
    
    # If you give it an array of dateTimes
    assert type(dateTime) is not list, 'dateTime input cannot be a list! Use numpy arrays'   
    if type(dateTime) is np.ndarray:
        gst = [Greenwich_sidereal_time(i) for i in dateTime]
    else:
        gst = Greenwich_sidereal_time(dateTime)
    if sphere_shape:
        theta = np.rad2deg(gst) + long
        x = (Re + alt)*np.cos(np.deg2rad(lat))*np.cos(np.deg2rad(theta))
        y = (Re + alt)*np.cos(np.deg2rad(lat))*np.sin(np.deg2rad(theta))
        z = (Re + alt)*np.sin(np.deg2rad(lat))
    else:
        # Code to calculate the subpoint for the oblate spheroid shape
        pass
    
    return np.array([x, y, z])
    
R = lambda z, alt: (Re + alt)*np.sqrt(1 - (z/Re)**2)
    
def Julian_date_of_year(year):
    """
    Code taken from Page 61 of Astronomical Algorithms by Jean Meeus or the
    implementation from https://celestrak.com/columns/v02n02/
    
    Calculates the juian date of January 0.0 the input year.
    
    This function works!
    """
    year -= 1
    A = year//100
    B = 2 - A + A//4
    return int(365.25 * year) + int(30.6001 * 14) + 1720994.5 + B
    
def Julian_date(date):
    """
    This function calculates the julian date using a UTC date time.
    date is an iso formated UTC date time.
    
    If date is a ISO formatted string, it will be parsed by 
    dateutil.parser.parse().
    
    This function works!
    """
    # Parse datetime if necessary
    if type(date) is not datetime.datetime:
        date = dateutil.parser.parse(date)
    dateTup = date.timetuple()
    doy = dateTup.tm_yday # day of year
    
    # fraction of the day
    dFrac = (date - datetime.datetime.combine(date.date(), \
    datetime.datetime.min.time())).total_seconds()/86400.0
    
    return Julian_date_of_year(dateTup.tm_year) + doy + dFrac

def Greenwich_sidereal_time(dateTime):
    #We = 7.29211510E-5 #radians/second

    # Parse datetime if necessary
    if type(dateTime) is not datetime.datetime:
        dateTime = dateutil.parser.parse(dateTime)
        
        # Find number of seconds in today.
    #dt = (dateTime - datetime.datetime.combine(dateTime.date(), \
    #datetime.datetime.min.time())).total_seconds()
    #print(dt)
    #dateTup = dateTime.timetuple()
    
    jd = Julian_date(dateTime)
    UT = (jd + 0.5) % 1 # Return the fractional part
    jd -= UT
    Tu = (jd - 2451545.0)/36525
    thetaG0 = 24110.54841 + 8640184.812866*Tu + 0.093104*Tu**2 - 6.2E-6*Tu**3
    # We*dt
    return 2*np.pi*((thetaG0 + 86400.0*1.00273790934*UT) % 86400)/86400.0
    
if __name__ == '__main__':
    print('\n Running ECI to GDZ coordinate transformation test harness. \n')
    dateTime = np.array(['1995-10-1T09:00:00', '1995-10-1T09:00:00'])
    # Now run an example! 
    lon = np.array([-75, -75])
    lat = np.array([40, 40])
    alt = np.array([0, 0])
    
    X = gdz2eci(lat, lon, alt, dateTime)
    GDZCoords = eci2gdz(X, dateTime)

    print('Convert from lat = ', lat, ' lon = ', lon, ' alt = ', alt, ' to ECI')
    print(X, '\n')
    print('Convert from ECI to GDZ')
    print(GDZCoords, '\n')
    
    print('Differences are (alt, lat long) = (', alt - GDZCoords['Alt'], \
    lat - GDZCoords['Lat'], lon - GDZCoords['Lon'], ')')
    
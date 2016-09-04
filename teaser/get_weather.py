#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 02 10:49:18 2016

@author: tsz
"""

from __future__ import division
import numpy as np

import sun

def get_weather(filename, beta, gamma, albedo=0.2, timeZone=1,
                altitude=0, location=(49.5, 8.5)):
    
    """
    Parse weather related data from TRY file

    Arguments
    ---------
    filename: string
              name of the TRY file to be parsed (not all of the data is used)
    beta: array-like
          angles between exterior wall areas and the ground
    gamma: array-like
           orientations of the exterior wall areas (N = 0, E = 90, S = 180, W = 270)
    
    Returns
    -------
    rad_sky: numpy ndarray
             radiation coming from sky
    rad_earth: numpy ndarray
               radiation coming from ground
    temp: numpy ndarray
          air temperature
    sun_rad: numpy ndarray
             radiation on tilted surface areas for each orientation
    """    
    
    weather_data = np.loadtxt(filename, skiprows=38, usecols=(8,13,14,16,17))
    temp = weather_data[:,0]        # temperature
    sun_dir = weather_data[:,1]     # direct sun radiation
    sun_diff = weather_data[:,2]    # diffuse sun radiation
    rad_sky = weather_data[:,3]     # irradiation from sky
    rad_earth = weather_data[:,4]   # irradiation from land surface
    
#    # Parse all lines 
#    with open(filename) as f:
#        # Read all lines at once
#        all_lines = f.readlines()
#        
#        # get location data
#        if re.match("Lage", all_lines[2]) != None:
#            i = all_lines[2].replace("<- B.","").replace("<- L.","").replace(" ","").split("N")
#            j = filter(None, re.split("[° \']+",i[0].replace("Lage:","")))
#            latitude = float(j[0])+float(j[1])/60
#            i = i[1].split("O")
#            j = filter(None, re.split("[° \']+",i[0]))
#            longitude = float(j[0])+float(j[1])/60
#            altitude = float(i[1].replace("Meter\xfcber",""))

    # Sun radiation on surfaces
    sun_rad  = sun.getSolarGains(0, 3600, weather_data.shape[0],
                                 timeZone = timeZone,
                                 location = location,
                                 altitude = altitude,
                                 beta     = beta,
                                 gamma    = gamma,
                                 beam     = sun_dir,
                                 diffuse  = sun_diff,
                                 albedo   = albedo)
    
    return rad_sky, rad_earth, temp + 273.15, sun_rad

def get_betaGamma(orientationswallshorizontal, offset = 0):
    
    beta = orientationswallshorizontal        
    n = len(beta)
    
    gamma = np.array([0,90,180,270]) + offset
    
    if n == 4:
        pass
    elif n == 5:
        gamma.append(0)
    elif n == 6:
        # in the current Teaser data file: beta = [45,90,90,45,90,90]
        gamma = -np.array([0,0,90,0,180,270])
        
    return beta, gamma

if __name__ == "__main__":
    filename = "TRY2010_12_Jahr.dat"
    beta, gamma = get_betaGamma([45, 90, 90, 45, 90, 90])
    res = get_weather(filename, beta, gamma)
    
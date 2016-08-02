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
    
    weather_data = np.loadtxt(filename, skiprows=38, usecols=(8,13,14,16,17))
    temp = weather_data[:,0]
    sun_dir = weather_data[:,1]
    sun_diff = weather_data[:,2]
    rad_sky = weather_data[:,3]
    rad_earth = weather_data[:,4]

    # Sun radiation on surfaces
    sun_rad  = sun.getSolarGains(0, 3600, weather_data.shape[0],
                                 timeZone=timeZone,
                                 location=location,
                                 altitude=altitude,
                                 beta=beta,
                                 gamma=gamma,
                                 beam = sun_dir,
                                 diffuse = sun_diff,
                                 albedo = albedo)
    
    return rad_sky, rad_earth, temp, sun_rad

if __name__ == "__main__":
    filename = "TRY2010_12_Jahr.dat"
    beta = [45, 90, 90, 45, 90, 90]
    gamma = -np.array([0, 0, 90, 0, 180, 270])
    res = get_weather(filename, beta, gamma)
    
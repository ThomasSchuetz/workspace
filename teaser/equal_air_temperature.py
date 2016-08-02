#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 02 12:52:24 2016

@author: tsz
"""
from __future__ import division

import numpy as np

def eqAirTempVDI(weatherData, houseData, solarRad_in, sunblindsig, withLongwave=True):      
#    Parameters partialEqAirTemp 
    aowo      = houseData["aowo"] # Coefficient of absorption of the outer walls
    eowo      = houseData["epso"] # Coefficient of emission   of the outer walls
    n         = len(houseData["orientationswallshorizontal"]) # Number of orientations (without ground)
    T_ground  = houseData["temperatureground"] # Temperature of the ground in contact with ground slab
    
    wf_wall   = houseData["weightfactorswall"]    # Weight factors of the walls
    wf_win    = houseData["weightfactorswindow"]  # Weight factors of the windows
    wf_ground = houseData["weightfactorground"]   # Weight factor of the ground (0 if not considered)
    
    T_air   = weatherData["temperature"] # outdoor air temperature
    E_sky   = weatherData["solar_irrad_sky"] # Irradiation from sky
    E_earth = weatherData["solar_irrad_earth"] # Irradiation from land surface
    
#    Parameters EqAirTempVDI
    alphaowo = houseData["alphaowo"] # Outer wall's coefficient of heat transfer (outer side)
    orientationswallshorizontal = np.array(houseData["orientationswallshorizontal"]) # orientations of the walls against the vertical (wall,roof)    
    
#    Initialization
    timesteps = np.shape(solarRad_in)[1]
    T_eqLW       = np.zeros((n,timesteps))
    T_eqWin      = np.zeros((n,timesteps))
    T_eqWall     = np.zeros((n,timesteps))    
    
#    Calculation of single-index values
    T_earth = ((-E_earth / (0.93 * 5.67)) ** 0.25) * 100 -273.15
    T_sky   = ((E_sky / (0.93 * 5.67)) ** 0.25) * 100    -273.15
    
    alpharad = (E_sky + (E_earth / 0.93)) / (T_sky - T_earth)
    alpharad[np.abs(E_sky + E_earth) < 0.1] = 5
    
    phiprivate = 0.5 * (np.ones(n) + np.cos(orientationswallshorizontal * (np.pi / 180)))
    
    T_eqSW = solarRad_in * aowo / alphaowo
    
#    Calculation of multi-index values
    for i in range(n):
       T_eqLW[i,:] = ((T_earth - T_air) * (1 - phiprivate[i]) + (T_sky - T_air) * phiprivate[i]) * eowo * alpharad / alphaowo
       
       if withLongwave:
           T_eqWin[i,:]  = T_air + T_eqLW[i,:] * np.abs(sunblindsig[i,:] - 1)
           T_eqWall[i,:] = T_air + T_eqLW[i,:] + T_eqSW[i,:]
       else:
           T_eqWin[i,:]  = T_air
           T_eqWall[i,:] = T_air + T_eqSW[i,:]

#    Equal air temperature
    equalAirTemp = np.dot(wf_wall, T_eqWall) + np.dot(wf_win, T_eqWin) + (T_ground - 273.15) * wf_ground

    return equalAirTemp, alpharad


if __name__ == "__main__":
    import pickle as pkl

    with open("inputs_eqAir.pkl", "rb") as fin:
        raw_inputs = pkl.load(fin)
        houseData = pkl.load(fin)
        solarRad_in = pkl.load(fin)
        sunblindsig = pkl.load(fin)

    equalAirTemp, alphaRad = eqAirTempVDI(raw_inputs, houseData, solarRad_in, sunblindsig)
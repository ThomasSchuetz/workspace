#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 02 11:59:46 2016

@author: tsz
"""
from __future__ import division
import numpy as np

def sunblinds(params, solarBeforeSunblinds):
    """ 
    Calculate solar input for reduced order model and for eqAirTemp through sunblinds
    
    Arguments
    ---------
    params: dictionary
            contains misc. input parameters
            - gWin: array-like
                    Total energy transmittances if sunblind is CLOSED
            - Imax: int/float
                    Intensity at which the sunblind closes
            - AWin: array-like
                    Window area for each orientation
    solarBeforeSunblinds: numpy ndarray
                          radiation before consideration of sunblinds
    
    Returns
    -------
    sunblindsig: array-like 
                 state of shading, reduction factor due to shading (0 if sunblinds are completely opened)
    """
    n           = len(params["gWin"]) # Number of orientations (without ground)
    I_max       = params["Imax"] # Intensity at which the sunblind closes
    gsunblind   = params["gWin"] # Total energy transmittances if sunblind is closed
    Aw          = params["AWin"] # Window area for each orientation
    
    T = solarBeforeSunblinds.shape[1]
    sunblindsig = np.zeros((n,T))
    afterSunblind = solarBeforeSunblinds

    for i in range(n):
        ind_greater_Imax = solarBeforeSunblinds[i,:] > I_max
        sunblindsig[i,ind_greater_Imax] = 1 - gsunblind[i]
        afterSunblind[i,ind_greater_Imax] = solarBeforeSunblinds[i,ind_greater_Imax] * gsunblind[i]

    # weight solar irradiation by window area
    sum_Aw = max(0.0001, sum(Aw))
    weighted_sum = np.dot(Aw, afterSunblind) / sum_Aw

    return weighted_sum, sunblindsig

if __name__ == "__main__":
    houseData = {"Imax": 150, 
                 "AWin": [1.56, 4.0, 7.64, 3.125, 12.6, 6.9],
                 "gWin": [1.0, 0.8, 0.4, 0.3, 1.0]}
    
    np.random.seed(0)
    solarBeforeSunblinds = np.random.randint(low=0,
                                             high=1000,
                                             size=(len(houseData["AWin"]), 8760))
    
    weighted_sum, sunblindsig = sunblinds(houseData, solarBeforeSunblinds)
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 02 12:52:24 2016

@author: tsz
"""
from __future__ import division

import numpy as np

def equal_air_temp(HSol, TBlaSky, TDryBul, sunblind, params):
    """
    Computes the equivalent air temperature on exterior walls without considering windows seperately
    (several approaches possible here)
    
    Arguments
    ---------
    HSol: numpy ndarray 
           solar radiation per unit area
    TBlaSky: numpy ndarray
             black-body sky temperature
    TDryBul: numpy ndarray
             dry bulb temperature
    sunblind: numpy ndarray
              opening factor of sunblinds for each direction (0 = open to 1 = closed, sunblinds.sunblindsig)
    params: dictionary
            misc. constant input parameters
            - eExt: float
                    coefficient of emission of exterior walls (outdoor)
            - aExt: float
                    coefficient of absorption of exterior walls (outdoor)
            - alpha_rad_wall: int/float
                              heat transfer coefficient
            - alpha_wall_out: int/float
                              heat transfer coefficient
            - wfWall: float
                      weight factors of the walls
            - wfWin: float
                     weight factors of the windows
            - wfGro: float
                     weight factor of the ground (0 if not considered)
            - T_Gro: int/float
                     constant ground temperature
            - withLongwave: boolean
                            True if longwave radiation is considered

    Returns
    -------
    TEqAir: numpy ndarray
            equivalent air temperature on exterior walls for convective heat entry
    TEqAirWin: numpy ndarray
               equivalent air temperature on exterior windows for convective heat entry
    """
    # Read parameters to improve readability in the equations
    eExt = params["eExt"] # coefficient of emission of exterior walls (outdoor)
    aExt = params["aExt"] # coefficient of absorption of exterior walls (outdoor)
    alphaRadWall = params["alpha_rad_wall"]
    alphaWallOut = params["alpha_wall_out"]
    eWin = params["eWin"] # coefficient of emission of windows
    aWin = params["aWin"] # coefficient of absorption of windows
    alphaRadWin = params["alpha_rad_win"]
    alphaWinOut = params["alpha_win_out"]
    wfWall = params["wfWall"] # weight factors of the walls
    wfWin = params["wfWin"] # weight factors of the windows
    wfGro = params["wfGro"] # weight factor of the ground (0 if not considered)
    TGro = params["T_Gro"] # ground temperature
    n = len(wfWall) # orientations
    
    # Compute equivalent long wave and short wave temperatures
    delTEqLW = (TBlaSky - TDryBul) * (eExt * alphaRadWall / (alphaRadWall + alphaWallOut))
    delTEqLWWin = (TBlaSky - TDryBul) * (eWin * alphaRadWin / (alphaRadWin + alphaWinOut))
    delTEqSW = HSol * aExt / (alphaRadWall + alphaWallOut)
    delTEqSWWin = HSol * aWin / (alphaRadWin + alphaWinOut)
    
    # Compute equivalent window and wall temperatures
    if params["withLongwave"]:
        TEqWin = np.array([TDryBul + (delTEqSWWin + delTEqLWWin) * (1 - sunblind[:,i]) for i in range(n)]).T
        TEqWall = np.array([TDryBul + delTEqLW[:,i] + delTEqSW[:,i] for i in range(n)]).T
    else:
        TEqWin = np.array([TDryBul + delTEqSWWin * (1 - sunblind[:,i]) for i in range(n)]).T
        TEqWall = np.array([TDryBul + delTEqSW[:,i] for i in range(n)]).T
    
    # Compute equivalent air temperature
    TEqAir = np.dot(TEqWall, wfWall) + TGro * wfGro
    TEqAirWin = np.dot(TEqWin, wfWin)
    
    # Return result
    return TEqAir, TEqAirWin

if __name__ == "__main__":
    times_per_hour = 60
    no_tile = 60
    
    t_outside_raw = np.loadtxt("inputs/case08_t_amb.csv", delimiter=",")
    t_outside = np.array([t_outside_raw[2*i,1] for i in range(24)])
    t_outside_adj = np.repeat(t_outside, times_per_hour)
    t_outside_tiled = np.tile(t_outside_adj, no_tile)
    
    q_sol_rad_win_raw = np.loadtxt("inputs/case08_q_sol_win.csv", usecols=(1,2))
    solarRad_win = q_sol_rad_win_raw[0:24,:]
    
    sunblind_in = np.zeros_like(solarRad_win)
    sunblind_in[solarRad_win > 100] = 0.85
    sunblind_in_adj = np.repeat(sunblind_in, times_per_hour, axis=0)
    
    q_sol_rad_wall_raw = np.loadtxt("inputs/case08_q_sol_wall.csv", usecols=(1,2))
    solarRad_wall = q_sol_rad_wall_raw[0:24,:]
    solarRad_wall_adj = np.repeat(solarRad_wall, times_per_hour, axis=0)
    solarRad_wall_tiled = np.tile(solarRad_wall_adj.T, no_tile).T
    
    t_black_sky = np.zeros_like(t_outside) + 273.15
    
    params = {"aExt": 0.7,
              "eExt": 0.9,
              "wfWall": [0.05796831135677373, 0.13249899738691134],
              "wfWin": [0.4047663456281575, 0.4047663456281575],
              "wfGro": 0,
              "T_Gro": 273.15+12,
              "alpha_wall_out": 20,
              "alpha_rad_wall": 5,
              "withLongwave": False}
    
    t_equal_air = equal_air_temp(solarRad_wall, t_black_sky, t_outside, sunblind_in, params)
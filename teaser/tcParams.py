#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 14:19:07 2016

@author: Markus
"""
from __future__ import division
import numpy as np

#%% Load a standard result file
def load_res(filename):
    res = np.loadtxt(filename, delimiter=",", skiprows=1) # Skip time step 0
    
    # ignore time
    result = res[:,1:res.shape[1]]
    
    day1 = result[0:24, :]
    day2 = result[24:48, :]
    day3 = result[48:72, :]
    
    return (day1, day2, day3)

#%% Common house inputs
def get_house_data(case=1):

    """
    Meaning of Parameters
    ---------------------
    R1i: resistor of inner wall
    C1i: capacity of inner wall
    Ai: inner wall area
    RRest: remaining Resistor between exterior wall and the ambient air
    R1o: resistor of exterior wall
    C1o: capacity of exterior wall
    Ao: area of exterior wall for each orientation
    Aw: window area for each orientation
    At: solar radiation transmitting area of windows for each orientation
    Vair: air volume
    rhoair: air density
    cair: thermal heat capacity of the air
    splitfac: factor splitting convective and radiative gains from solar radiation
    g: energy transmittance
    alphaiwi: heat transfer coefficient from inner wall to air zone
    alphaowi: heat transfer coefficient from outer wall to air zone
    alphaWall: heat transfer coefficient from ambient air to outer wall
    withInnerwalls: True, if inner walls are considered
    """    
    
    if case in (1,2):
        return {"R1i": 0.000595693407511, 
                "C1i": 14836354.6282, 
                "Ai": 75.5, 
                "RRest": 0.03895919557, 
                "R1o": 0.00436791293674, 
                "C1o": 1600848.94,
                "Ao": [10.5], 
                "Aw": np.zeros(1), 
                "At": np.zeros(1), 
                "Vair": 52.5, 
                "rhoair": 1.19, 
                "cair": 0,
                "splitfac": 0.09,
                "g": 1,
                "alphaiwi": 2.24,
                "alphaowi": 2.7,
                "alphaWall": 25 * 10.5, # 25 * sum(Ao)
                "withInnerwalls": True}
    elif case in (3,4):
        return {"R1i": 0.003237138, 
                "C1i": 7297100, 
                "Ai": 75.5, 
                "RRest": 0.039330865, 
                "R1o": 0.00404935160802, 
                "C1o": 47900,
                "Ao": [10.5], 
                "Aw": np.zeros(1), 
                "At": np.zeros(1), 
                "Vair": 52.5, 
                "rhoair": 1.19, 
                "cair": 0,
                "splitfac": 0.09,
                "g": 1,
                "alphaiwi": 2.24,
                "alphaowi": 2.7,
                "alphaWall": 25 * 10.5, # 25 * sum(Ao)
                "withInnerwalls": True}
    elif case in (5,12):
        return {"R1i": 0.000595693407511, 
                "C1i": 14836354.6282, 
                "Ai": 75.5, 
                "RRest": 0.03895919557, 
                "R1o": 0.00436791293674, 
                "C1o": 1600848.94,
                "Ao": [10.5], 
                "Aw": np.zeros(1), 
                "At": [7], 
                "Vair": 0, 
                "rhoair": 1.19, 
                "cair": 1007,
                "splitfac": 0.09,
                "g": 1,
                "alphaiwi": 2.24,
                "alphaowi": 2.7,
                "alphaWall": 25 * 10.5, # 25 * sum(Ao)
                "withInnerwalls": True}
    elif case in (6,):
        return {"R1i": 0.000595515, 
                "C1i": 14836200, 
                "Ai": 75.5, 
                "RRest": 0.038959197, 
                "R1o": 0.004367913, 
                "C1o": 1600800,
                "Ao": [10.5], 
                "Aw": np.zeros(1), 
                "At": [0], 
                "Vair": 52.5, 
                "rhoair": 1.19, 
                "cair": 0,
                "splitfac": 0.09,
                "g": 1,
                "alphaiwi": 2.24,
                "alphaowi": 2.7,
                "alphaWall": 25 * 10.5, # 25 * sum(Ao)
                "withInnerwalls": True}
    elif case in (7,):
        return {"R1i": 0.000595693407511, 
                "C1i": 14836354.6282, 
                "Ai": 75.5, 
                "RRest": 0.03895919557, 
                "R1o": 0.00436791293674, 
                "C1o": 1600848.94,
                "Ao": [10.5], 
                "Aw": np.zeros(1), 
                "At": [0], 
                "Vair": 52.5, 
                "rhoair": 1.19, 
                "cair": 0,
                "splitfac": 0.09,
                "g": 1,
                "alphaiwi": 2.24,
                "alphaowi": 2.7,
                "alphaWall": 25 * 10.5, # 25 * sum(Ao)
                "withInnerwalls": True}
    elif case in (8,9):
        return {"R1i": 0.000668895639141, 
                "C1i": 12391363.8631, 
                "Ai": 60.5, 
                "RRest": 0.01913729904, 
                "R1o": 0.0017362530106, 
                "C1o": 5259932.23,
                "Ao": [10.5,15], 
                "Aw": np.zeros(2), 
                "At": [7,7], 
                "Vair": 52.5, 
                "rhoair": 1.19, 
                "cair": 0,
                "splitfac": 0.09,
                "g": 1,
                "alphaiwi": 2.12,
                "alphaowi": 2.7,
                "alphaWall": 25 * 25.5, # 25 * sum(Ao)
                "withInnerwalls": True}
    elif case in (10,):
        return {"R1i": 0.000779671554640369, 
                "C1i": 12333949.4129606, 
                "Ai": 58, 
                "RRest": 0.011638548, 
                "R1o": 0.00171957697767797, 
                "C1o": 4338751.41,
                "Ao": [28], 
                "Aw": np.zeros(1), 
                "At": [7,], 
                "Vair": 52.5, 
                "rhoair": 1.19, 
                "cair": 0,
                "splitfac": 0.09,
                "g": 1,
                "alphaiwi": 2.12,
                "alphaowi": 2.398,
                "alphaWall": 28 * 9.75, # 9.75 * sum(Ao)
                "withInnerwalls": True}
    elif case in (11,):
        return {"R1i": 0.000595693407511, 
                "C1i": 14836354.6282, 
                "Ai": 75.5, 
                "RRest": 0.03895919557, 
                "R1o": 0.00436791293674, 
                "C1o": 1600848.94,
                "Ao": [10.5], 
                "Aw": np.zeros(1), 
                "At": [0], 
                "Vair": 0, 
                "rhoair": 1.19, 
                "cair": 1007,
                "splitfac": 0.09,
                "g": 1,
                "alphaiwi": 3,
                "alphaowi": 2.7,
                "alphaWall": 25 * 10.5, # 25 * sum(Ao)
                "withInnerwalls": True}

#%%
def get_eqAirTemp_params(case=8):
    if case in (8,9):
        return {"aExt": 0.7,
                "eExt": 0.9,
                "wfWall": [0.05796831135677373, 0.13249899738691134],
                "wfWin": [0.4047663456281575, 0.4047663456281575],
                "wfGro": 0,
                "T_Gro": 273.15 + 12,
                "alpha_wall_out": 20,
                "alpha_rad_wall": 5,
                "withLongwave": False}
    elif case in (10,):
        return {"aExt": 0.7,
                "eExt": 0.9,
                "wfWall": [0.04646093176283288,],
                "wfWin": [0.32441554918476245,],
                "wfGro": 0.6291235190524047,
                "T_Gro": 273.15 + 15,
                "alpha_wall_out": 20,
                "alpha_rad_wall": 5,
                "withLongwave": False}
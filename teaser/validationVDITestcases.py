#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 14:19:07 2016

@author: Markus
"""
from __future__ import division
import numpy as np
import LowOrderModelVDIvalidation as LOM

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
                "alphaWall": 28 * 9.75, # 28 * sum(Ao)
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


#%%
def testCase1(timesteps, times_per_hour=60, n=4):
    
    timesteps_day = int(24 * times_per_hour)
    
    withInnerwalls = True
    R1i = 0.000595693407511
    C1i = 14836354.6282
    Ai = 75.5
    alphaiwi = 2.24
    epsi = 1
    alphaRad = np.zeros(timesteps) + 5    

    withWindows = False
    RWin = 0.00000001
    splitfac = 0.09
    Aw = np.zeros(n)
    A_win_tot = 0
    epsw = 1
    g = 1
    
    withOuterwalls = True
    RRest = 0.03895919557
    R1o = 0.00436791293674
    C1o = 1600848.94
    Ao = 10.5
    alphaowi = 2.7
    epso = 1
    alphaWall = 25 * Ao
    
    Vair = 52.5
    rhoair = 1.19
    cair = 0
    
    Tv = np.zeros(timesteps) + 295.15 # in K
    ventRate = np.zeros(timesteps)
    
    solarRad_in = np.zeros(timesteps)   
    source_igRad = np.zeros(timesteps)
    
    Q = {}
    Q["ig"] = np.zeros(timesteps_day)
    for q in range(int(6*timesteps_day/24), int(18*timesteps_day/24)):
        Q["ig"][q] = 1000
    Q["ig"] = np.tile(Q["ig"], 60)
        
    equalAirTemp = np.zeros(timesteps) + 295.15 # all temperatures in K
    
    (T_air1, T_air10, T_air60) = load_res("inputs/case01_res.csv")
    
    return   (R1i, C1i, Ai, RRest, R1o, C1o, Ao, RWin, Aw, A_win_tot,
             alphaiwi, epsi, alphaowi, epso, alphaWall, alphaRad,
             Tv, ventRate, solarRad_in, Vair, rhoair, cair,
             source_igRad, Q["ig"], equalAirTemp,
             withInnerwalls, withWindows, withOuterwalls,
             splitfac, epsw, g, T_air1[:,0], T_air10[:,0], T_air60[:,0])
             
#%%
def testCase2(houseData, n, timesteps, model, testcase):
    
    withInnerwalls = True
    R1i = 0.000595515
    C1i = 1.4836200e7/3600
    Ai = 75.5
    alphaiwi = 2.23642384
    epsi = 1
    alphaRad = np.zeros(timesteps)+5    

    withWindows = False
    RWin = 0
    splitfac = 0.09
    Aw = np.zeros(n)
    A_win_tot = 0
    epsw = 1
    g = 1
    
    withOuterwalls = True
    RRest = 0.042768721
    R1o = 0.004367913
    C1o = 1.6008e6/3600
    Ao = 10.5
    alphaowi = 2.7
    epso = 1
    
    Vair = 52.5
    rhoair = 1.19
    cair = 1007/3600    
    
    Tv = np.zeros(timesteps) + 295.15 # in K
    ventRate = 0
    
    solarRad_in = np.zeros(timesteps)    

    source_igRad = np.zeros(timesteps)
    for q in xrange(6,18):
        source_igRad[q] = 1e3
    
    Q = {}
    Q["ig"] = np.zeros(timesteps)
        
    equalAirTemp = np.zeros(timesteps)+295.15 # all temperatures in K
    C_HC = 0
    Tig = np.zeros(timesteps)
    
    return   (R1i, C1i, Ai, RRest, R1o, C1o, Ao, RWin, Aw, A_win_tot,
             alphaiwi, epsi, alphaowi, epso, alphaRad,
             Tv, ventRate, solarRad_in, Vair, rhoair, cair,
             source_igRad, Q["ig"], equalAirTemp,
             withInnerwalls, withWindows, withOuterwalls,
             splitfac, epsw, g, C_HC, Tig, model)
             
#%%
def testCase3(houseData, n, timesteps, model, testcase):
    
    withInnerwalls = True
    R1i = 0.003237138
    C1i = 7.297100e6/3600
    Ai = 75.5
    alphaiwi = 2.23642384
    epsi = 1
    alphaRad = np.zeros(timesteps)+5    

    withWindows = False
    RWin = 0
    splitfac = 0.09
    Aw = np.zeros(n)
    A_win_tot = 0
    epsw = 1
    g = 1
    
    withOuterwalls = True
    RRest = 0.043140385
    R1o = 0.004049352
    C1o = 4.79e4/3600
    Ao = 10.5
    alphaowi = 2.7
    epso = 1
    
    Vair = 52.5
    rhoair = 1.19
    cair = 1007/3600    
    
    Tv = np.zeros(timesteps) + 295.15 # in K
    ventRate = 0
    
    solarRad_in = np.zeros(timesteps)    
    source_igRad = np.zeros(timesteps)
    
    Q = {}
    Q["ig"] = np.zeros(timesteps)
    for q in xrange(6,18):
        Q["ig"][q] = 1e3
        
    equalAirTemp = np.zeros(timesteps)+295.15 # all temperatures in K
    C_HC = 0
    Tig = np.zeros(timesteps)
    
    return   (R1i, C1i, Ai, RRest, R1o, C1o, Ao, RWin, Aw, A_win_tot,
             alphaiwi, epsi, alphaowi, epso, alphaRad,
             Tv, ventRate, solarRad_in, Vair, rhoair, cair,
             source_igRad, Q["ig"], equalAirTemp,
             withInnerwalls, withWindows, withOuterwalls,
             splitfac, epsw, g, C_HC, Tig, model)
             
#%%
def testCase4(houseData, n, timesteps, model, testcase):
    
    withInnerwalls = True
    R1i = 0.003237138
    C1i = 7.297100e6/3600
    Ai = 75.5
    alphaiwi = 2.23642384
    epsi = 1
    alphaRad = np.zeros(timesteps)+5    

    withWindows = False
    RWin = 0
    splitfac = 0.09
    Aw = np.zeros(n)
    A_win_tot = 0
    epsw = 1
    g = 1
    
    withOuterwalls = True
    RRest = 0.043140385
    R1o = 0.004049352
    C1o = 4.79e4/3600
    Ao = 10.5
    alphaowi = 2.7
    epso = 1
    
    Vair = 52.5
    rhoair = 1.19
    cair = 1007/3600    
    
    Tv = np.zeros(timesteps) + 295.15 # in K
    ventRate = 0
    
    solarRad_in = np.zeros(timesteps)    
    source_igRad = np.zeros(timesteps)
    for q in xrange(6,18):
        source_igRad[q] = 1e3
        
    Q = {}
    Q["ig"] = np.zeros(timesteps)
        
    equalAirTemp = np.zeros(timesteps)+295.15 # all temperatures in K
    C_HC = 0
    Tig = np.zeros(timesteps)
    
    return   (R1i, C1i, Ai, RRest, R1o, C1o, Ao, RWin, Aw, A_win_tot,
             alphaiwi, epsi, alphaowi, epso, alphaRad,
             Tv, ventRate, solarRad_in, Vair, rhoair, cair,
             source_igRad, Q["ig"], equalAirTemp,
             withInnerwalls, withWindows, withOuterwalls,
             splitfac, epsw, g, C_HC, Tig, model)
             
#%%
def testCase5(houseData, n, timesteps, model, testcase):
    
    n = 5    
    
    withInnerwalls = True
    R1i = 0.000595515
    C1i = 1.48362e7/3600
    Ai = 75.5
    alphaiwi = 2.23642384
    epsi = 1
    alphaRad = np.zeros(timesteps)+5    

    withWindows = True
    RWin = 0
    splitfac = 0.09
    Aw = np.zeros(n)
    Aw[2] = 7
    A_win_tot = 7
    epsw = 1
    g = 1
    
    withOuterwalls = True
    RRest = 0.042768721
    R1o = 0.004367913
    C1o = 1.6008e6/3600
    Ao = 10.5
    alphaowi = 2.7
    epso = 1
    
    Vair = 52.5
    rhoair = 1.19
    cair = 1007/3600    
    
    Tv = np.zeros(timesteps) + 295.15 # in K
    ventRate = 0
    
    solarsunblinds = np.zeros((n,timesteps))
    solarsunblinds[2,5:20] = np.array([17,38,59,98,186,287,359,385,359,287,186,98,59,38,17])
    gsunblind = np.array([0,0,0.15,0,0])
    I_max = 100
    
    solarRad_in, sunblindsig = LOM.sundblinds(houseData,solarsunblinds,n,I_max,gsunblind,Aw,testcase)   
    source_igRad = np.zeros(timesteps)
    Q = {}
    Q["ig"] = np.zeros(timesteps)
    for q in xrange(8,18):
        source_igRad[q] = 80
        Q["ig"][q] = 280 
        
    equalAirTemp = np.zeros(timesteps)+295.15 # all temperatures in K
    equalAirTemp[:24] = np.array([291.95,290.25,289.65,289.25,289.65,290.95,293.45,295.95,297.95,299.85,301.25,302.15,302.85,303.55,304.05,304.15,303.95,303.25,302.05,300.15,297.85,296.05,295.05,294.05,])
    C_HC = 0
    Tig = np.zeros(timesteps)
    
    return   (R1i, C1i, Ai, RRest, R1o, C1o, Ao, RWin, Aw, A_win_tot,
             alphaiwi, epsi, alphaowi, epso, alphaRad,
             Tv, ventRate, solarRad_in, Vair, rhoair, cair,
             source_igRad, Q["ig"], equalAirTemp,
             withInnerwalls, withWindows, withOuterwalls,
             splitfac, epsw, g, C_HC, Tig, model)
             
#%%
def testCase6(houseData, n, timesteps, model, testcase):    
    
    withInnerwalls = True
    R1i = 0.000595515
    C1i = 1.48362e7/3600
    Ai = 75.5
    alphaiwi = 2.23642384
    epsi = 1
    alphaRad = np.zeros(timesteps)+5    

    withWindows = False
    RWin = 0
    splitfac = 0.03
    Aw = np.zeros(n)
    A_win_tot = 7
    epsw = 0.95
    g = 0.15
    
    withOuterwalls = True
    RRest = 0.042768721
    R1o = 0.004367913
    C1o = 1.6008e6/3600
    Ao = 10.5
    alphaowi = 2.7
    epso = 1
    
    Vair = 0.001
    rhoair = 1.19
    cair = 1007/3600    
    
    Tv = np.zeros(timesteps) + 295.15 # in K
    ventRate = 0
    
    solarRad_in = np.zeros(timesteps)   
    source_igRad = np.zeros(timesteps)
    Q = {}
    Q["ig"] = np.zeros(timesteps)
    for q in xrange(6,18):
        source_igRad[q] = 1000 
        
    equalAirTemp = np.zeros(timesteps)+295.15 # all temperatures in K
    C_HC = 1
    Tig = np.zeros(timesteps)
    Tig[:25] = np.array([295.1,295.1,295.1,295.1,295.1,295.1,295.1,300.1,300.1,300.1,300.1,300.1,300.1,300.1,300.1,300.1,300.1,300.1,300.1,295.1,295.1,295.1,295.1,295.1,295.1])
    
    Tair = {}
    for t in xrange(25):
        Tair[t] = model.addVar(vtype="C", name="Tair_"+str(t), lb=-100.)
    model.update()
    for t in xrange(25):
        model.addConstr(Tair[t] == Tig[t],
                        name="init_Tair_testcase6_"+str(t))
    model.update()
    
    return   (R1i, C1i, Ai, RRest, R1o, C1o, Ao, RWin, Aw, A_win_tot,
             alphaiwi, epsi, alphaowi, epso, alphaRad,
             Tv, ventRate, solarRad_in, Vair, rhoair, cair,
             source_igRad, Q["ig"], equalAirTemp,
             withInnerwalls, withWindows, withOuterwalls,
             splitfac, epsw, g, C_HC, Tig, model)
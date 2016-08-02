#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 14:19:07 2016

@author: Markus
"""
from __future__ import division
import numpy as np
import LowOrderModelVDIvalidation as LOM

def testCase1(houseData, n, timesteps, model, testcase):
    
    withInnerwalls = True
    R1i = 0.000595515
    C1i = 1.48362e7/3600
    Ai = 75.5
    alphaiwi = 2.236423594
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
    alphaowi = 2.699999358
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
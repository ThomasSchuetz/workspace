#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 02 13:59:12 2016

@author: tsz
"""

from __future__ import division

import numpy as np

def reducedOrderModelVDI(houseData, weatherTemperature, solarRad_in, equalAirTemp, alphaRad, ventRate,
                         Q_ig, source_igRad, krad, t_set_heating, t_set_cooling, dt=1,
                         T_air_init=293.15, T_iw_init=295.15, T_ow_init=295.15):
#%% partialReducedOrderModel                                          
    # parameters
    
    timesteps       = len(alphaRad)

    withInnerwalls  = houseData["withInnerwalls"]   # If inner walls are existent
#    withWindows     = houseData["withWindows"]      # If windows are existent
#    withOuterwalls  = houseData["withOuterwalls"]   # If outer walls (including windows) are existent
    
    R1i             = houseData["R1i"] # Resistor 1 inner wall
    C1i             = houseData["C1i"] # Capacity 1 inner wall in Wh/K
    Ai              = houseData["Ai"] # Inner wall area
    RRest           = houseData["RRest"] # Resistor Rest outer wall
    R1o             = houseData["R1o"] # Resistor 1 outer wall
    C1o             = houseData["C1o"] # Capacity 1 outer wall in Wh/K
    Ao              = houseData["Ao"] # Outer wall area
    Aw              = houseData["Aw"] # Window area
    
    Vair            = houseData["Vair"]     # Volume of the air in the zone
    rhoair          = houseData["rhoair"]   # Density of the air
    cair            = houseData["cair"] # Heat capacity of the air in Wh/KgK
    splitfac        = houseData["splitfac"] # Factor for conv. part of rad. through windows
    
    g               = houseData["g"]    # Total energy transmittance
    
    alphaiwi        = houseData["alphaiwi"] # Coefficient of heat transfer for inner walls
    alphaowi        = houseData["alphaowi"] # Outer wall's coefficient of heat transfer (inner side)
    
    A_win_tot = sum(Aw)

#%% Time variable inputs
    # convective heat entry from solar iradiation   
    Q_solar_conv = solarRad_in * splitfac * g * A_win_tot
    
    # splitters:
    # on each splitter: one output goes to outer wall, one goes to inner wall
    # therefor dimension is 2 if inner walls exist => 2 outgoing signals
    
    # therm. splitter solar radiative:
    Q_solar_rad = solarRad_in * (splitfac - 1) * g * A_win_tot
    
    # therm. splitter loads radiative:
    Q_loads_rad = - krad * source_igRad

    # each splitter goes into inner (if true) and outer walls
    # split heat transfer according to area ratio
    if withInnerwalls:
        dim = 2
        splitFacSolar = [(Ao - A_win_tot) / (Ao + Ai - A_win_tot), Ai/(Ao + Ai - A_win_tot)]
        splitFacLoads = [Ao / (Ao + Ai), Ai / (Ao + Ai)]
        Q_solarRadToInnerWall = -splitFacSolar[-dim+1] * Q_solar_rad
        Q_loadsToInnerWall    = -splitFacLoads[-dim+1] * Q_loads_rad
    else:
        dim = 1
        splitFacSolar = [(Ao - A_win_tot) / (Ao + Ai - A_win_tot)]
        splitFacLoads = [Ao / (Ao + Ai)]
        Q_solarRadToInnerWall = 0
        Q_loadsToInnerWall    = 0
        
    Q_solarRadToOuterWalli = -splitFacSolar[-dim] * Q_solar_rad
    Q_loadsToOuterWalli    = -splitFacLoads[-dim] * Q_loads_rad
    
#%%     Define system of linear equations: 
    # A * x = rhs
    # x = [T_ow, T_owi, T_iw, T_iwi, T_air, Q_air, Q_HC] (all at time t)
    
    # Results
    T_ow = []
    T_owi = []
    T_iw = []
    T_iwi = []
    T_air = []
    Q_air = []
    Q_HC = []
    
    # Dummy time step
    T_ow_prev = T_ow_init
    T_iw_prev = T_iw_init
    T_air_prev = T_air_init

    for t in range(timesteps):
        # Common equations
        A = np.zeros((7,7))
        rhs = np.zeros(A.shape[0])
    
        A[0,0] = C1o / dt + 1 / RRest + 1 / R1o
        A[0,1] = -1 / R1o
        A[1,0] = 1 / R1o
        A[1,1] = -Ao * (alphaRad[t] + alphaowi) - 1 / R1o
        A[1,3] = Ao * alphaRad[t]
        A[1,4] = Ao * alphaowi
        A[2,2] = C1i / dt + 1 / R1i
        A[2,3] = -1 / R1i
        A[3,1] = Ao * alphaRad[t]
        A[3,2] = 1 / R1i
        A[3,3] = -Ao * alphaRad[t] - Ai * alphaiwi - 1 / R1i
        A[3,4] = Ai * alphaiwi
        A[4,1] = Ao * alphaowi
        A[4,3] = Ai * alphaiwi
        A[4,4] = -Ao * alphaowi - Ai * alphaiwi - ventRate[t] * Vair * cair * rhoair
        A[4,5] = -1
        A[4,6] = 1
        A[5,4] = ventRate[t] * Vair * cair * rhoair
        A[5,5] = -1
        
        rhs[0] = equalAirTemp[t] / RRest + C1o * T_ow_prev / dt
        rhs[1] = -Q_solarRadToOuterWalli[t] - Q_loadsToOuterWalli[t]
        rhs[2] = C1i * T_iw_prev / dt
        rhs[3] = -Q_solarRadToInnerWall[t] - Q_loadsToInnerWall[t]
        rhs[4] = ventRate[t] * Vair * cair * rhoair * Tv[t] + Q_solar_conv[t] + Q_ig[t]
        rhs[5] = rhoair * cair * Vair * T_air_prev / dt
    
        x = _calc_timestep(A, rhs, t_set_heating[t], t_set_cooling[t])

        T_ow.append(x[0])
        T_owi.append(x[1])
        T_iw.append(x[2])
        T_iwi.append(x[3])
        T_air.append(x[4])
        Q_air.append(x[5])
        Q_HC.append(x[6])
        
        T_ow_prev = x[0]
        T_iw_prev = x[2]
        T_air_prev = x[4]
    
    return (T_air, Q_HC)

#%%
def _calc_timestep(A, rhs, t_set_heating=291.15, t_set_cooling=300.15):
    
    # Calculate without further heat inputs to determine if heating 
    # or cooling is needed
    # x = [T_ow, T_owi, T_iw, T_iwi, T_air, Q_air, Q_HC]
    x_noHeat = _calc_temperatue(A, rhs, q_hc_fix=0)
    
    if x_noHeat[4] < t_set_heating:
        return _calc_heatflow(A, rhs, t_air_set=t_set_heating)
    elif x_noHeat[4] > t_set_cooling:
        return _calc_heatflow(A, rhs, t_air_set=t_set_cooling)
    else:
        return x_noHeat

#%%
def _calc_temperatue(A, rhs, q_hc_fix=0):

    # Delete all entries in the final line of A:
    A[-1,:] = 0
    
    # Add Q_HC = q_hc_fix
    A[6,6] = 1
    rhs[6] = q_hc_fix
    
    result = np.linalg.solve(A, rhs)
    
    return result

#%%
def _calc_heatflow(A, rhs, t_air_set=293.15):
    # Delete all entries in the final line of A:
    A[-1,:] = 0
    
    # Add T_air = t_air_set
    A[6,4] = 1
    rhs[6] = t_air_set
    
    result = np.linalg.solve(A, rhs)
    
    return result

#%%
if __name__ == "__main__":
#    import pickle as pkl
#    with open("inputs_ROM_VDI.pkl", "rb") as fin:
#        houseData = pkl.load(fin)
#        weatherTemperature = pkl.load(fin)
#        solarRad_in = pkl.load(fin)
#        equalAirTemp = pkl.load(fin)
#        alphaRad = pkl.load(fin)
    
    time_steps = 24
    import validationVDITestcases as tc
    (R1i, C1i, Ai, RRest, R1o, C1o, Ao, RWin, Aw, A_win_tot,
     alphaiwi, epsi, alphaowi, epso, alphaRad,
     Tv, ventRate, solarRad_in, Vair, rhoair, cair,
     source_igRad, Q_ig, equalAirTemp,
     withInnerwalls, withWindows, withOuterwalls,
     splitfac, epsw, g) = tc.testCase1(time_steps)
     
    krad = 1
    
    houseData = {"R1i":R1i, "C1i":C1i, "Ai":Ai, "RRest":RRest, "R1o":R1o, "C1o":C1o,
                 "Ao":Ao, "Aw":Aw, "Vair":Vair, "rhoair":rhoair, "cair":cair,
                 "splitfac":splitfac, "g":g, "alphaiwi":alphaiwi, "alphaowi":alphaowi,
                 "withInnerwalls":True}
    
    weatherTemperature = Tv
    
    t_set_heating   = np.zeros(time_steps)# + 293.15
    t_set_cooling   = np.zeros(time_steps)+600# + 300.15
    
    T_air, Q_hc = reducedOrderModelVDI(houseData, weatherTemperature, solarRad_in,
                                       equalAirTemp, alphaRad, ventRate, Q_ig, source_igRad, krad,
                                       t_set_heating, t_set_cooling)
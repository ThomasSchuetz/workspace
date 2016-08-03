#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 02 13:59:12 2016

@author: tsz
"""

from __future__ import division

import numpy as np

def reducedOrderModelVDI(houseData, weatherTemperature, solarRad_in, equalAirTemp, alphaRad, dt=1,
                         T_air_init=293.15, T_iw_init=293.15, T_ow_init=293.15):
#%% partialReducedOrderModel                                          
    # parameters
    
    n               = len(houseData["orientationswallshorizontal"]) # Number of orientations (without ground)
    T               = range(1,len(alphaRad))
    timesteps       = len(alphaRad)

    withInnerwalls  = houseData["withInnerwalls"]   # If inner walls are existent
    withWindows     = houseData["withWindows"]      # If windows are existent
    withOuterwalls  = houseData["withOuterwalls"]   # If outer walls (including windows) are existent
    
    R1i             = houseData["R1i"] # Resistor 1 inner wall
    C1i             = houseData["C1i"] / 3600 # Capacity 1 inner wall in Wh/K
    Ai              = houseData["Ai"] # Inner wall area
    RRest           = houseData["RRest"] # Resistor Rest outer wall
    R1o             = houseData["R1o"] # Resistor 1 outer wall
    C1o             = houseData["C1o"] / 3600 # Capacity 1 outer wall in Wh/K
    Ao              = houseData["Ao"] # Outer wall area
    Aw              = houseData["Aw"] # Window area
    
    Vair            = houseData["Vair"]     # Volume of the air in the zone
    rhoair          = houseData["rhoair"]   # Density of the air
    cair            = houseData["cair"] / 3600 # Heat capacity of the air in Wh/KgK
    splitfac        = houseData["splitfac"] # Factor for conv. part of rad. through windows
    
    g               = houseData["g"]    # Total energy transmittance
    
    alphaiwi        = houseData["alphaiwi"] # Coefficient of heat transfer for inner walls
    alphaowi        = houseData["alphaowi"] # Outer wall's coefficient of heat transfer (inner side)
    
    krad            = 0.4
    ventRate        = 1
    Tv              = weatherTemperature + 273.15
    source_igRad    = np.zeros(timesteps) + 200

    Q = {}
    Q["ig"]         = np.zeros(timesteps) + 120
    
    A_win_tot = sum(Aw)
    
    equalAirTemp += 273.15
    C_HC = 0
    Tig = np.zeros(timesteps)
    
#%% convective heat entry from solar iradiation   
    k = splitfac * g * A_win_tot
    y = k*solarRad_in

    Q["solarConv"] = y
    
    Q_solar_conv = solarRad_in * splitfac * g * A_win_tot
    
#%% splitters:
    # on each splitter: one output goes to outer wall, one goes to inner wall
    # therefor dimension is 2 if inner walls exist => 2 outgoing signals
    
    # therm. splitter solar radiative:
    k = (1-splitfac)*g*A_win_tot
    y = k*solarRad_in
    
    Q["solarRad"] = -y
    
    Q_solar_rad = solarRad_in * (splitfac - 1) * g * A_win_tot
    
    # therm. splitter loads radiative:
    y = krad * source_igRad
    
    Q["loadsRad"] = -y
    
    Q_loads_rad = - krad * source_igRad

    # each splitter goes into inner (if true) and outer walls
    # split heat transfer according to area ratio
    if withInnerwalls:
        dim = 2
        splitFacSolar = [(Ao - A_win_tot) / (Ao + Ai - A_win_tot), Ai/(Ao + Ai - A_win_tot)]
        splitFacLoads = [Ao / (Ao + Ai), Ai / (Ao + Ai)]
        Q_solarRadToInnerWall   = -splitFacSolar[-dim+1] * Q["solarRad"]
        Q_loadsToInnerWall      = -splitFacLoads[-dim+1] * Q["loadsRad"]
    else:
        dim = 1
        splitFacSolar = [(Ao - A_win_tot) / (Ao + Ai - A_win_tot)]
        splitFacLoads = [Ao / (Ao + Ai)]
        Q_solarRadToInnerWall   = 0
        Q_loadsToInnerWall      = 0
        
    Q_solarRadToOuterWalli   = -splitFacSolar[-dim] * Q["solarRad"]
    Q_loadsToOuterWalli      = -splitFacLoads[-dim] * Q["loadsRad"]
    
    
    
    # Define system of linear equations: 
    # A * x = rhs
    # x = [T_ow, T_owi, T_iw, T_iwi, T_air, Q_air, Q_HC] (all at time t)
    
    # Dummy time step
    t = 0
    T_ow_prev = T_ow_init
    T_iw_prev = T_iw_init
    T_air_prev = T_air_init
    
    q_hc_fix = 0
    t_air_set = 294
    
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
    A[4,4] = -Ao * alphaowi - Ai * alphaiwi - ventRate * Vair * cair * rhoair
    A[4,5] = -1
    A[4,6] = 1
    A[5,4] = ventRate * Vair * cair * rhoair
    A[5,5] = -1
    
    rhs[0] = equalAirTemp[t] / RRest + C1o * T_ow_prev / dt
    rhs[1] = -Q_solarRadToOuterWalli[t] - Q_loadsToOuterWalli[t]
    rhs[2] = C1i * T_iw_prev / dt
    rhs[3] = -Q_solarRadToInnerWall[t] - Q_loadsToInnerWall[t]
    rhs[4] = ventRate * Vair * cair * rhoair * Tv[t] + Q["solarConv"][t] + Q["ig"][t]
    rhs[5] = rhoair * cair * Vair * T_air_prev / dt
    
    
    # Fixed heat flow, e.g. 0:    
    A[6,6] = 1
    rhs[6] = q_hc_fix
    
    # Fixed set temperature:
    A[6,4] = 1
    rhs[6] = t_air_set
    
    
#%% add variables to the model
    """
    variables:
        Tair - achieved room temperature
        Towi - inside temperature of the outer wall
        Tiwi - inside temperature of the inner wall
        Tow  - temperature of the outer wall (outside temperature is the equalAirTemp)
        Tiw  - temperature of the inner wall
        Q_HC - heat or cooling coming from storage
        Qair - heat/cooling transfer going into the room
    """
#    Tair = {}
#    Towi = {}
#    Tiwi = {}
#    Tow  = {}
#    Tiw  = {}
#    Q_HC = {}
#    Qair = {}
#    
#    Tair[0]   = model.getVarByName("Tair_0")
#    Tow[0]    = model.getVarByName("Tow_0")
#    Tiw[0]    = model.getVarByName("Tiw_0")
#    Q_HC[0]   = model.getVarByName("Q_HC_0")
#    
#    for t in T:    
#        if model.getVarByName("Tair_"+str(t)) in model.getVars():
#            Tair[t] = model.getVarByName("Tair_"+str(t))
#        else:
#            Tair[t] = model.addVar(vtype="C", name="Tair_"+str(t), lb=-100.)        
#        Tow[t]  = model.addVar(vtype="C", name="Tow_"+str(t),  lb=-100.)
#        Tiw[t]  = model.addVar(vtype="C", name="Tiw_"+str(t),  lb=-100.)
#        
#        Towi[t]   = model.addVar(vtype="C", name="Towi_"+str(t), lb=-100.)
#        Tiwi[t]   = model.addVar(vtype="C", name="Tiwi_"+str(t), lb=-100.)
#        
##        Q_HC[timestep]   = model.addVar(vtype="C", name="Q_HC_"+str(timestep), lb=-1e5)
#        Q_HC[t]   = model.getVarByName("Q_HC_"+str(t))
#        Qair[t]   = model.addVar(vtype="C", name="Qair_"+str(t), lb=-1e5)
#    
#    model.update()
#    
##%% add main contraints: balances
##    for t in T:
#    testperiod = range(1,49)
#    for t in testperiod:
#        # outer wall
#        model.addConstr(C1o * (Tow[t]-Tow[t-1])/dt == (equalAirTemp[t]-Tow[t])/RRest - (Tow[t]-Towi[t])/R1o,
#                        name="outerWall_transfer_"+str(t))
#        model.addConstr(Q_solarRadToOuterWalli[t] + Q_loadsToOuterWalli[t] - 
#                        Ao*alphaRad[t]*(Towi[t]-Tiwi[t]) + Ao*alphaowi*(Tair[t]-Towi[t]) + 
#                        (Tow[t]-Towi[t])/R1o == 0,
#                        name="outerWall_balance_"+str(t))
#        
#        # inner wall
#        model.addConstr(C1i * (Tiw[t]-Tiw[t-1])/dt == (Tiwi[t]-Tiw[t])/R1i,
#                        name="innerWall_transfer_"+str(t))
#        model.addConstr(Q_solarRadToInnerWall[t] + Q_loadsToInnerWall[t] +
#                        Ao*alphaRad[t]*(Towi[t]-Tiwi[t]) + Ai*alphaiwi*(Tair[t]-Tiwi[t]) - 
#                        (Tiwi[t]-Tiw[t])/R1i == 0,
#                        name="innerWall_balance_"+str(t))
#        
#        # room
#        model.addConstr(Qair[t] == -Ao*alphaowi*(Tair[t]-Towi[t]) - Ai*alphaiwi*(Tair[t]-Tiwi[t]) + 
#                        ventRate*Vair*cair*rhoair*(Tv[t]-Tair[t]) + Q["solarConv"][t] + Q["ig"][t] + Q_HC[t],
#                        name="room_balance_"+str(t))
#        
#        # additional constraints
#        model.addConstr(Qair[t] == (Tair[t]-Tair[t-1])/dt * rhoair*cair*Vair,
#                        name="room_transfer_"+str(t))
#        
#        # bounds
#    #    bounds = {}
#    #    bounds["Tair_up"]  = 300
#    #    bounds["Tair_low"] = 288
#    #    Tair[t] <= bounds["Tair_up"]
#    #    Tair[t] >= bounds["Tair_low"]
#                        
#    model.update()
#    model.write("lgs.lp")
#                        

#%%
if __name__ == "__main__":
    import pickle as pkl
    with open("inputs_ROM_VDI.pkl", "rb") as fin:
        houseData = pkl.load(fin)
        weatherTemperature = pkl.load(fin)
        solarRad_in = pkl.load(fin)
        equalAirTemp = pkl.load(fin)
        alphaRad = pkl.load(fin)
    
    reducedOrderModelVDI(houseData, weatherTemperature, solarRad_in, equalAirTemp, alphaRad)
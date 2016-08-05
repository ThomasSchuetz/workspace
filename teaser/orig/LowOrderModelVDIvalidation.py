#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 10:20:57 2016

@author: Markus
"""

from __future__ import division

import math
import parser_buildings
import numpy as np
import re
import gurobipy as gp
import validationVDITestcases as valid

def getTRYData(houseData):
    
    rad_sky     = []
    rad_earth   = []
    temp        = []
    sun_dir     = []
    sun_diff    = []
    
    # Parse all lines 
    with open("TRY2010_12_Jahr.dat") as f:
        # Read all lines at once
        all_lines = f.readlines()
        
        # Skip the first 38 lines
        for i in range(38, len(all_lines)):
            
            s = all_lines[i].split()

            rad_sky.append(float(s[16]))
            rad_earth.append(float(s[17]))
            temp.append(float(s[8]))
            sun_dir.append(float(s[13]))
            sun_diff.append(float(s[14]))
        
        # get location data
        if re.match("Lage", all_lines[2]) != None:
            i = all_lines[2].replace("<- B.","").replace("<- L.","").replace(" ","").split("N")
            j = filter(None, re.split("[° \']+",i[0].replace("Lage:","")))
            latitude = float(j[0])+float(j[1])/60
            i = i[1].split("O")
            j = filter(None, re.split("[° \']+",i[0]))
            longitude = float(j[0])+float(j[1])/60
            altitude = float(i[1].replace("Meter\xfcber",""))
            
#        location = (latitude, longitude)
        location = (49.5, 8.5)
        altitude = 0
        
        # calculate orientation depending iradiation
        import sun
        timeZone = 1
        albedo = 0.2
        
        beta = houseData["orientationswallshorizontal"]        
        n = len(beta)
        
        offset = 0
        
        gamma = np.array([0,90,180,270])
        gamma += offset
        if n == 4:
            pass
        elif n == 5:
            gamma.append(0)
        elif n == 6:
            # in the current Teaser data file: beta = [45,90,90,45,90,90]
            gamma = -np.array([0,0,90,0,180,270])
        
        # Sun radiation on surfaces
        SunRad  = sun.getSolarGains(0, 3600, 8760, 
                                    timeZone=timeZone,
                                    location=location,
                                    altitude=altitude,
                                    beta=beta,
                                    gamma=gamma,
                                    beam = np.array(sun_dir),
                                    diffuse = np.array(sun_diff),
                                    albedo = albedo)
    
    return np.array(rad_sky), np.array(rad_earth), np.array(temp), SunRad

def eqAirTempVDI(weatherData, houseData, solarRad_in, sunblindsig):      
#%% partialEqAirTemp 
    aowo            = houseData["aowo"] # Coefficient of absorption of the outer walls
    eowo            = houseData["epso"] # Coefficient of emission   of the outer walls
    n               = len(houseData["orientationswallshorizontal"]) # Number of orientations (without ground)
    T_ground        = houseData["temperatureground"] # Temperature of the ground in contact with ground slab
    withLongwave    = True # If longwave radiation exchange is considered
    
    wf_wall     = houseData["weightfactorswall"]    # Weight factors of the walls
    wf_win      = houseData["weightfactorswindow"]  # Weight factors of the windows
    wf_ground   = houseData["weightfactorground"]   # Weight factor of the ground (0 if not considered)
    
    unitvec     = np.ones(n)
    
    T_air   = weatherData[0] # outdoor air temperature
    E_sky   = weatherData[1] # Iradiation from sky
    E_earth = weatherData[2] # Iradiation from land surface
    
    #%% EqAirTempVDI
    alphaowo    = houseData["alphaowo"] # Outer wall's coefficient of heat transfer (outer side)
    orientationswallshorizontal = np.array(houseData["orientationswallshorizontal"]) # orientations of the walls against the vertical (wall,roof)    
    
    timesteps = np.size(E_earth)
    alpharad     = np.zeros(timesteps)
    equalAirTemp = np.zeros(timesteps)
    T_earth      = np.zeros(timesteps)
    T_sky        = np.zeros(timesteps)

    T_eqLW       = np.zeros((n,timesteps))
    T_eqSW       = np.zeros((n,timesteps))
    T_eqWin      = np.zeros((n,timesteps))
    T_eqWall     = np.zeros((n,timesteps))    
    
    for t in xrange(timesteps):
        # scalars
        T_earth[t] = ((-E_earth[t]/(0.93*5.67))**0.25)*100 -273.15
        T_sky[t]   = ((E_sky[t]/(0.93*5.67))**0.25)*100    -273.15    
        
        if abs(E_sky[t]+E_earth[t])<0.1:
            alpharad[t] = 5.0
        else:
            alpharad[t] = (E_sky[t]+(E_earth[t]/0.93))/(T_sky[t]-T_earth[t])        
        
        phiprivate = np.zeros(n)
        phiprivate = unitvec*0.5 + np.cos(orientationswallshorizontal*(math.pi/180))*0.5
        
        T_eqLW[:,t] = ((T_earth[t]-T_air[t])*(unitvec-phiprivate)+(T_sky[t]-T_air[t])*phiprivate)*(eowo*alpharad[t]/alphaowo)
        T_eqSW[:,t] = solarRad_in[:,t]*(aowo/alphaowo)
           
        if withLongwave == True:
            T_eqWin[:,t]  = T_air[t]*unitvec+T_eqLW[:,t]*abs(sunblindsig[:,t]-unitvec)
            T_eqWall[:,t] = T_air[t]*unitvec+T_eqLW[:,t]+T_eqSW[:,t]
        else:
            T_eqWin[:,t]  = T_air[t]*unitvec
            T_eqWall[:,t] = T_air[t]*unitvec+T_eqSW[:,t]
        
        # scalar products        
        equalAirTemp[t] = np.dot(T_eqWall[:,t],wf_wall) + np.dot(T_eqWin[:,t],wf_win) + (T_ground-273.15)*wf_ground
    
    return equalAirTemp, alpharad

def reducedOrderModelVDI(houseData, weatherTemperature, solarRad_in,
                         equalAirTemp, alphaRad, model, testcase=False):
#%% partialReducedOrderModel                                          
    # parameters
    
    Q = {}                    
    if testcase == True:
        n               = len(houseData["orientationswallshorizontal"]) # Number of orientations (without ground)
        T               = xrange(1,len(alphaRad))
        timesteps       = len(alphaRad)
        
        (R1i, C1i, Ai, RRest, R1o, C1o, Ao, RWin, Aw, A_win_tot,
         alphaiwi, epsi, alphaowi, epso, alphaRad,
         Tv, ventRate, solarRad_in, Vair, rhoair, cair,
         source_igRad, Q["ig"], equalAirTemp,
         withInnerwalls, withWindows, withOuterwalls,
         splitfac, epsw, g, C_HC, Tig, model) = valid.testCase1(houseData,n,timesteps,model,testcase)
         
        krad = 1
         
    else:
        n               = len(houseData["orientationswallshorizontal"]) # Number of orientations (without ground)
        T               = xrange(1,len(alphaRad))
        timesteps       = len(alphaRad)
    
        withInnerwalls  = houseData["withInnerwalls"]   # If inner walls are existent
        withWindows     = houseData["withWindows"]      # If windows are existent
        withOuterwalls  = houseData["withOuterwalls"]   # If outer walls (including windows) are existent
        
        R1i             = houseData["R1i"] # Resistor 1 inner wall
        C1i             = houseData["C1i"]/3600 # Capacity 1 inner wall in Wh/K
        Ai              = houseData["Ai"] # Inner wall area
        RRest           = houseData["RRest"] # Resistor Rest outer wall
        R1o             = houseData["R1o"] # Resistor 1 outer wall
        C1o             = houseData["C1o"]/3600 # Capacity 1 outer wall in Wh/K
        Ao              = houseData["Ao"] # Outer wall area
        Aw              = houseData["Aw"] # Window area
        RWin            = houseData["RWin"]
        
        Vair            = houseData["Vair"]     # Volume of the air in the zone
        rhoair          = houseData["rhoair"]   # Density of the air
        cair            = houseData["cair"]/3600 # Heat capacity of the air in Wh/KgK
        splitfac        = houseData["splitfac"] # Factor for conv. part of rad. through windows
        
        epsw            = houseData["epsw"] # Emissivity of the windows
        g               = houseData["g"]    # Total energy transmittance
        
        alphaiwi        = houseData["alphaiwi"] # Coefficient of heat transfer for inner walls
        alphaowi        = houseData["alphaowi"] # Outer wall's coefficient of heat transfer (inner side)
        alphaConvWinInner = houseData["alphaConvWinInner"]
        epsi            = houseData["epsi"]     # Emissivity of the inner walls
        epso            = houseData["epso"]     # Emissivity of the outer walls
        
        krad            = 0.4
        ventRate        = 1
        Tv              = weatherTemperature + 273.15
        source_igRad    = np.zeros(timesteps)+200
        Q["ig"]         = np.zeros(timesteps)+120
        
        A_win_tot = sum(Aw[i] for i in xrange(n))
        
        equalAirTemp += 273.15
        C_HC = 0
        Tig = np.zeros(timesteps)
    
#%% convective heat entry from solar iradiation   
    k = splitfac * g * A_win_tot
    y = k*solarRad_in

    Q["solarConv"] = y
    
#%% splitters:
    # on each splitter: one output goes to outer wall, one goes to inner wall
    # therefor dimension is 2 if inner walls exist => 2 outgoing signals
    
    # therm. splitter solar radiative:
    k = (1-splitfac)*g*A_win_tot
    y = k*solarRad_in
    
    Q["solarRad"] = -y
    
    # therm. splitter loads radiative:
    y = krad * source_igRad
    
    Q["loadsRad"] = -y

    # each splitter goes into inner (if true) and outer walls
    # split heat transfer according to area ratio
    if withInnerwalls == True:
        dim = 2
        splitFacSolar = [(Ao - A_win_tot)/(Ao + Ai - A_win_tot), Ai/(Ao + Ai - A_win_tot)]
        splitFacLoads = [Ao / (Ao + Ai), Ai / (Ao + Ai)]
        Q_solarRadToInnerWall   = -splitFacSolar[-dim+1]*Q["solarRad"]
        Q_loadsToInnerWall      = -splitFacLoads[-dim+1]*Q["loadsRad"]
    else:
        dim = 1
        splitFacSolar = [(Ao - A_win_tot)/(Ao + Ai - A_win_tot)]
        splitFacLoads = [Ao / (Ao + Ai)]
        Q_solarRadToInnerWall   = 0
        Q_loadsToInnerWall      = 0
        
    Q_solarRadToOuterWalli   = -splitFacSolar[-dim]*Q["solarRad"]
    Q_loadsToOuterWalli      = -splitFacLoads[-dim]*Q["loadsRad"]
    
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
    Tair = {}
    Towi = {}
    Tiwi = {}
    Tow  = {}
    Tiw  = {}
    Q_HC = {}
    Qair = {}
    
    Tair[0]   = model.getVarByName("Tair_0")
    Tow[0]    = model.getVarByName("Tow_0")
    Tiw[0]    = model.getVarByName("Tiw_0")
    Q_HC[0]   = model.getVarByName("Q_HC_0")
    
    for t in T:    
        if model.getVarByName("Tair_"+str(t)) in model.getVars():
            Tair[t] = model.getVarByName("Tair_"+str(t))
        else:
            Tair[t] = model.addVar(vtype="C", name="Tair_"+str(t), lb=-100.)        
        Tow[t]  = model.addVar(vtype="C", name="Tow_"+str(t),  lb=-100.)
        Tiw[t]  = model.addVar(vtype="C", name="Tiw_"+str(t),  lb=-100.)
        
        Towi[t]   = model.addVar(vtype="C", name="Towi_"+str(t), lb=-100.)
        Tiwi[t]   = model.addVar(vtype="C", name="Tiwi_"+str(t), lb=-100.)
        
#        Q_HC[timestep]   = model.addVar(vtype="C", name="Q_HC_"+str(timestep), lb=-1e5)
        Q_HC[t]   = model.getVarByName("Q_HC_"+str(t))
        Qair[t]   = model.addVar(vtype="C", name="Qair_"+str(t), lb=-1e5)
    
    model.update()
    
#%% add main contraints: balances
#    for t in T:
    testperiod = range(1,49)
    for t in testperiod:
        # outer wall
        model.addConstr(C1o * (Tow[t]-Tow[t-1])/dt == (equalAirTemp[t]-Tow[t])/RRest - (Tow[t]-Towi[t])/R1o,
                        name="outerWall_transfer_"+str(t))
        model.addConstr(Q_solarRadToOuterWalli[t] + Q_loadsToOuterWalli[t] - 
                        Ao*alphaRad[t]*(Towi[t]-Tiwi[t]) + Ao*alphaowi*(Tair[t]-Towi[t]) + 
                        (Tow[t]-Towi[t])/R1o == 0,
                        name="outerWall_balance_"+str(t))
        
        # inner wall
        model.addConstr(C1i * (Tiw[t]-Tiw[t-1])/dt == (Tiwi[t]-Tiw[t])/R1i,
                        name="innerWall_transfer_"+str(t))
        model.addConstr(Q_solarRadToInnerWall[t] + Q_loadsToInnerWall[t] +
                        Ao*alphaRad[t]*(Towi[t]-Tiwi[t]) + Ai*alphaiwi*(Tair[t]-Tiwi[t]) - 
                        (Tiwi[t]-Tiw[t])/R1i == 0,
                        name="innerWall_balance_"+str(t))
        
        # room
        model.addConstr(Qair[t] == -Ao*alphaowi*(Tair[t]-Towi[t]) - Ai*alphaiwi*(Tair[t]-Tiwi[t]) + 
                        ventRate*Vair*cair*rhoair*(Tv[t]-Tair[t]) + Q["solarConv"][t] + Q["ig"][t] + Q_HC[t],
                        name="room_balance_"+str(t))
        
        # additional constraints
        model.addConstr(Qair[t] == (Tair[t]-Tair[t-1])/dt * rhoair*cair*Vair,
                        name="room_transfer_"+str(t))
        
        # bounds
    #    bounds = {}
    #    bounds["Tair_up"]  = 300
    #    bounds["Tair_low"] = 288
    #    Tair[t] <= bounds["Tair_up"]
    #    Tair[t] >= bounds["Tair_low"]
                        
    model.update()
    model.write("lgs.lp")
                        
#%% get testing results
    
    model.setParam("Aggregate", 0)
                    
    # Check for infeasible
#    model.computeIIS()
#    model.write("model.ilp")
#    print('\nConstraints:')        
#    for c in model.getConstrs():
#        if c.IISConstr:
#            print('%s' % c.constrName)
#    print('\nBounds:')
#    for v in model.getVars():
#        if v.IISLB > 0 :
#            print('Lower bound: %s' % v.VarName)
#        elif v.IISUB > 0:
#            print('Upper bound: %s' % v.VarName)
                                 
    model.optimize()

    res = {}
    res["Tair"] = {}
    res["Tiwi"] = {}
    res["Towi"] = {}
    res["Tow"]  = {}
    res["Tiw"]  = {}
    res["Qair"] = {}
    res["Qsrc"] = {}
    res["Qiwi"] = {}
    res["Qowi"] = {}
    res["Qowc"] = {}
    res["Qiwc"] = {}
    res["QHT"]  = {}
    res["QHC"]  = {}
    res["Qv"]   = {}
    res["Qig"]  = {}
    res["Teq"]  = {}
    res["QTSW1"]= {}
    res["QTSL1"]= {}
    res["QTSW2"]= {}
    res["QTSL2"]= {}

    for t in testperiod:
#    for t in T:
        res["Tair"][t] = Tair[t].X-273.15
        res["Tiwi"][t] = Tiwi[t].X-273.15
        res["Towi"][t] = Towi[t].X-273.15
        res["Tow"][t]  = Tow[t].X-273.15
        res["Tiw"][t]  = Tiw[t].X-273.15
        
        res["Qair"][t] = Qair[t].X
        res["Qsrc"][t] = Q["solarConv"][t]
        res["Qiwi"][t] = (Tiwi[t].X-Tiw[t].X)/R1i 
        res["Qowi"][t] = (Tow[t].X-Towi[t].X)/R1o
        res["Qowc"][t] = (Tair[t].X-Towi[t].X)*alphaowi*Ao
        res["Qiwc"][t] = (Tair[t].X-Tiwi[t].X)*alphaiwi*Ai
        res["QHT"][t]  = (Towi[t].X-Tiwi[t].X)*alphaRad[t]*Ao
        res["QHC"][t]  = Q_HC[t].X
        res["Qv"][t]   = ventRate*rhoair*cair*Vair*(Tv[t]-Tair[t].X)
        res["Qig"][t]  = Q["ig"][t]
        res["Teq"][t]  = equalAirTemp[t]-273.15
        
        res["QTSW1"][t]= Q_solarRadToOuterWalli[t]
        res["QTSL1"][t]= Q_loadsToOuterWalli[t]
        res["QTSW2"][t]= Q_solarRadToInnerWall[t]
        res["QTSL2"][t]= Q_loadsToInnerWall[t]

#%% return model with added variables/constraints
    return model
    
def sundblinds(houseData, solarBeforeSunblinds,
               n = 5,
               I_max = 150,
               gsunblind = np.ones(5),
               Aw = np.zeros(5),
               testcase = False):
#%% calculate solar input for reduced order model and for eqAirTemp through sunblinds
    
    if testcase == False:    
        n = len(houseData["orientationswallshorizontal"]) # Number of orientations (without ground)
        I_max = houseData["Imax"] # Intensity at which the sunblind closes
        gsunblind = houseData["gsunblind"] # Total energy transmittances if sunblind is closed
        Aw = houseData["Aw"] # Window area
    
    solar = {}

    T = np.size(solarBeforeSunblinds,1)
    sunblindsig = np.zeros((n,T))
    solar["afterSunblind"] = np.zeros((n,T))

    for t in xrange(T):
        for i in xrange(n):
            if solarBeforeSunblinds[i,t] > I_max:
                sunblindsig[i,t] = 1-gsunblind[i]
                solar["afterSunblind"][i,t] = solarBeforeSunblinds[i,t] * gsunblind[i]
            else:
                sunblindsig[i,t] = 0
                solar["afterSunblind"][i,t] = solarBeforeSunblinds[i,t]
    solar["RadWinTrans"] = solar["afterSunblind"]
    
    # weight solar iradiation by window area
    weightfactors = Aw
    if sum(weightfactors[i] for i in xrange(n)) < 0.0001:
        temp = 0.0001
    else: 
        temp = sum(weightfactors)
    solar["weightedSum"] = np.zeros(T)
    for t in xrange(T):
        solar["weightedSum"][t] = sum([solar["RadWinTrans"][i,t] * weightfactors[i] for i in xrange(n)]) / temp
    
    solarRad_in_redOrdMod = solar["weightedSum"]

    return solarRad_in_redOrdMod, sunblindsig

def addConstraints(model, dt, testcase):
#%% 
    # get building inputs
    filename = "TEASER4_meine_Geo.mo"
    houseData = parser_buildings.parse_record(filename)
    
    # get weather inputs
    raw_inputs = {}
 
    (raw_inputs["solar_irrad_sky"],
    raw_inputs["solar_irrad_earth"],
    raw_inputs["temperature"], solarRad_in) = getTRYData(houseData)
    
    # calculate solar input for reducedOrderModel (ROM) and rad through sunblinds for equalAirTemp (EAT)
    solarRad_in_redOrdMod, sunblindsig = sundblinds(houseData, solarRad_in, testcase)

    # calculate equivalent air temperature on outer walls
    equalAirTemp, alphaRad = eqAirTempVDI([raw_inputs["temperature"],
                                          raw_inputs["solar_irrad_sky"],
                                          raw_inputs["solar_irrad_earth"]],
                                          houseData,
                                          solarRad_in,
                                          sunblindsig)
    
    model = reducedOrderModelVDI(houseData,
                                 raw_inputs["temperature"],
                                 solarRad_in_redOrdMod,
                                 equalAirTemp,
                                 alphaRad,
                                 model,
                                 testcase)
                                     
    return model, equalAirTemp

def createUseCase(model):
#%%    
    Tair = {}
    Tow  = {}
    Tiw  = {}
    Q_HC = {}
    
    Tair[0] = model.addVar(vtype="C", name="Tair_0", lb=-100.)
    Tow[0]  = model.addVar(vtype="C", name="Tow_0",  lb=-100.)
    Tiw[0]  = model.addVar(vtype="C", name="Tiw_0",  lb=-100.)
    for i in xrange(8760):
        Q_HC[i] = model.addVar(vtype="C", name="Q_HC_"+str(i), lb=-1e5)

    model.update()

    T0all = 295.15
    model.addConstr(Tair[0] == 293.15,
                    name="initial_Tair")
    model.addConstr(Tow[0]  == T0all,
                    name="initial_Tow")
    model.addConstr(Tiw[0]  == T0all,
                    name="initial_Tiw")
    for i in xrange(8760):
        model.addConstr(Q_HC[i] == 0,
                        name="initial_HC_"+str(i))
        
    model.update()
    
    return model
    
#%%
if __name__ == "__main__":
    model = gp.Model("test")
    
    model = createUseCase(model)    
    
    dt = 1
    model, temp = addConstraints(model, dt, testcase=True)
    
    import matplotlib.pyplot as plt
    plt.plot(temp)
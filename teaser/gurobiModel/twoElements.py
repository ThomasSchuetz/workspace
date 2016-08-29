#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 10:20:57 2016

@author: Markus
"""

from __future__ import division

import numpy as np
import gurobipy as gp

#%%
def splitFacVal(nCol, A_array, AExt, AWin):
    """
    This is a direct translation from RC --> BaseClasses --> splitFacVal
    
    Parameters
    ----------
    nCol : int
        Number of orientations
    A_array : list
        [ATotExt, ATotWin]
    AExt : list
        Vector of exterior wall areas
    AWin : list
        Vector of window areas
    
    Example
    -------
    >>> # Define areas
    >>> AExt = [10.5]
    >>> AWin = [0]
    >>> A_int = 75.5
    >>> A_ar = [sum(AExt), sum(AWin), AInt]
    >>> # Calculate split factors for inner walls and outside walls
    >>> splitFac_IW = splitFacVal(dim, 1, A_ar, [0], [0])
    >>> splitFac_OW = splitFacVal(dim, len(AExt), A_ar, AExt, AWin)
    """
    
    A_tot = sum(A_array) # total area

    rows = len(A_array)

    # Counters
    i = 0 # A_array
    j = 0 # Row
    k = 0 # Column

    result = np.zeros((rows, nCol))    
    
    for A in A_array:
        if A > 0:
            k = 0
            if i == 0:
                for A_wall in AExt:
                    result[j,k] = (A - A_wall) / (A_tot - A_wall - AWin[k])
                    k += 1
            elif i == 1:
                for A_wall in AExt:
                    result[j,k] = (A - AWin[k]) / (A_tot - A_wall - AWin[k])
                    k += 1
            else:
                for A_wall in AExt:
                    result[j,k] = A / (A_tot - A_wall - AWin[k])
                    k += 1
            j += 1
        i += 1
    
    return result


def twoElements(params, solRad, window, extWall, windowIndoorSurface, extWallIndoorSurface,
                intWallIndoorSurface, intGainsConv, intGainsRad, ports, model, dt):
#%% # parameters
    # Thermal zone                
    Vair                = params["Vair"] # Air volume in the zone
    alphaRad            = params["alphaRad"] # Coefficient of heat transfer for linearized radiation exchange between walls
    n                   = params["nOrientations"] # Number of orientations

    # Windows
    AWin                = params["AWin"] # Vector of areas of windows by orientations
    ATransparent        = params["ATransparent"] # Vector of areas of transparent elements by orientations
    alphaWin            = params["alphaWin"] # Convective coefficient of heat transfer of windows (indoor)
    RWin                = params["RWin"] # Resistor for windows
    gWin                = params["gWin"] # Total energy transmittance of windows
    ratioWinConRad      = params["ratioWinConRad"] # Ratio for windows between indoor convective and radiative heat emission
    indoorPortWin       = params["indoorPortWin"]
    
    # Exterior walls
    AExt                = params["AExt"] # Vector of areas of exterior walls by orientations
    alphaExt            = params["alphaExt"] # Convective coefficient of heat transfer of exterior walls (indoor)
    nExt                = params["nExt"] # Number of RC-elements of exterior walls
    RExt                = params["RExt"] # Vector of resistances of exterior walls, from inside to outside
    RExtRem             = params["RExtRem"] # Resistance of remaining resistor RExtRem between capacity n and outside
    CExt                = params["CExt"]/3600 # Vector of heat capacities of exterior walls, from inside to outside, in Wh/K
    indoorPortExtWalls  = params["indoorPortExtWalls"]
    
    # Interior walls
    AInt                = params["AInt"] # Area of interior walls
    alphaInt            = params["alphaInt"] # Convective coefficient of heat transfer of interior walls (indoor)
    nInt                = params["nInt"] # Number of RC-elements of interior walls
    RInt                = params["RInt"] # Vector of resistances of interior walls, from port to center
    CInt                = params["CInt"]/3600 # Vector of heat capacities of interior walls, from port to center
    indoorPortIntWalls  = params["indoorPortIntWalls"]
    
    rhoair              = ports["rhoair"] # in Kg/m3
    cair                = ports["cair"] # in Wh
    ventRate            = ports["ventRate"]
    Tv                  = ports["Tv"] # in K
    
    timesteps           = len(alphaRad)
    
#%% solar radiation   
    solar_radiation = {}
    # convective heat entry from solar radiation
    solar_radiation["eConvSol"] = np.zeros((timesteps, len(ATransparent)))
    for i in range(len(ATransparent)):
        solar_radiation["eConvSol"][:,i] = gWin*ratioWinConRad*ATransparent[i] * solRad[:,i]
    solar_radiation["sumSolRad"] = np.sum(solar_radiation["eConvSol"], axis=1)
    solar_radiation["convHeatSol"] = -solar_radiation["sumSolRad"]
    Q_solarConv = solar_radiation["convHeatSol"]

    # radiative heat entry from solar radiation
    # split radiation onto exterior walls, interior walls (and windows)
    splitFacSolar = splitFacVal(len(AExt), [sum(AExt), sum(AWin), AInt], AExt, AWin)
    solar_radiation["eRadSol"] = np.zeros((timesteps, len(ATransparent)))
    for i in range(len(ATransparent)):
        solar_radiation["eRadSol"] = gWin*(1-ratioWinConRad)*ATransparent[i] * solRad[:,i]
    solar_radiation["radHeatSol"] = -solar_radiation["eRadSol"]

    solar_radiation["thermSplit"] = np.zeros((timesteps, len(AExt), splitFacSolar.shape[0]))
    for i in range(len(AExt)):
        for j in range(splitFacSolar.shape[0]):
            solar_radiation["thermSplit"][:,i,j] = solar_radiation["radHeatSol"] * splitFacSolar[j,i]
    
    Q_solarRadExt = np.sum(solar_radiation["thermSplit"][:,:,0], axis=1)
    Q_solarRadInt = np.sum(solar_radiation["thermSplit"][:,:,1], axis=1)
    Q_solarRadWin = np.sum(solar_radiation["thermSplit"][:,:,2], axis=1)
    
#%% internal gains
    # convective heat entry from internal gains
    Q_igConv = intGainsConv
    
    # radiative heat entry from internal gains
    # split radiation onto exterior walls, interior walls (and windows)
    splitFacIg = splitFacVal(1, [sum(AExt), sum(AWin), AInt], [0], [0])
    
    Q_igRadExt = intGainsRad * splitFacIg[0,0]
    Q_igRadInt = intGainsRad * splitFacIg[1,0]
    Q_igRadWin = intGainsRad * splitFacIg[2,0]
    
#%% add variables to the model
    """
    variables:
        Tair - achieved room temperature
        Towi - inside temperature of the outer wall
        Tiwi - inside temperature of the inner wall
        Twin - inside temperature of the window
        Tow  - temperature of the outer wall (outside temperature is the equalAirTemp)
        Tiw  - temperature of the inner wall
        Q_HC - heat or cooling coming from storage
        Qair - heat/cooling transfer going into the room
        
    All temperatures in Kelvin
    """
    Tair = {}
    Towi = {}
    Tiwi = {}
    Twin = {}
    Tow  = {}
    Tiw  = {}
    Q_HC = {}
    Qair = {}
    
    for t in range(timesteps):    
        if t > 0:
            if ports["setAirTemp"]:
                Tair[t] = model.getVarByName("Tair_"+str(t))
            else:
                Tair[t] = model.addVar(vtype="C", name="Tair_"+str(t), lb=-100.)
            Tow[t]  = model.addVar(vtype="C", name="Tow_"+str(t),  lb=-100.)
            Tiw[t]  = model.addVar(vtype="C", name="Tiw_"+str(t),  lb=-100.)
        else:
            Tair[0]   = model.getVarByName("Tair_0")
#            if model.getVarByName("Tow_0") in model.getVars():
            Tow[0]    = model.getVarByName("Tow_0")
#            if model.getVarByName("Tiw_0") in model.getVars():
            Tiw[0]    = model.getVarByName("Tiw_0")            
        
        Towi[t] = model.addVar(vtype="C", name="Towi_"+str(t), lb=-100.)
        Tiwi[t] = model.addVar(vtype="C", name="Tiwi_"+str(t), lb=-100.)
        Twin[t] = model.addVar(vtype="C", name="Twin_"+str(t), lb=-100.)
        
        Qair[t]   = model.addVar(vtype="C", name="Qair_"+str(t), lb=-1e5)
        Q_HC[t]   = model.addVar(vtype="C", name="Q_HC_"+str(t), lb=-1e5)
    
    model.update()
    
#%% add main contraints: balances
#    testperiod = range(1,49)
    for t in range(1,timesteps+1):
        if ports["heaterCooler"]:
            pass
        else:
            model.addConstr(Q_HC[t] == 0, name="prevent_HC_"+str(t))
        
        # exterior wall
        model.addConstr(CExt * (Tow[t]-Tow[t-1])/dt == (extWall[t]-Tow[t])/RExtRem - (Tow[t]-Towi[t])/RExt,
                        name="exteriorWall_transfer_"+str(t))
        model.addConstr(Q_solarRadExt[t] + Q_igRadExt[t] - 
                        min(sum(AExt),AInt)*alphaRad[t]*(Towi[t]-Tiwi[t]) + 
                        sum(AExt)*alphaExt*(Tair[t]-Towi[t]) + 
                        (Tow[t]-Towi[t])/RExt + 
                        min(sum(AExt),sum(AWin))*alphaRad[t]*(Twin[t]-Towi[t]) +
                        extWallIndoorSurface[t] == 0,
                        name="exteriorWall_balance_"+str(t))
        
        # interior wall
        model.addConstr(CInt * (Tiw[t]-Tiw[t-1])/dt == (Tiwi[t]-Tiw[t])/RInt,
                        name="interiorWall_transfer_"+str(t))
        model.addConstr(Q_solarRadInt[t] + Q_igRadInt[t] +
                        min(sum(AExt),AInt)*alphaRad[t]*(Towi[t]-Tiwi[t]) + 
                        AInt*alphaInt*(Tair[t]-Tiwi[t]) + 
                        (Tiw[t]-Tiwi[t])/RInt +
                        min(sum(AWin),AInt)*alphaRad[t]*(Twin[t]-Tiwi[t]) +
                        intWallIndoorSurface[t] == 0,
                        name="interiorWall_balance_"+str(t))
                        
        model.update()
        model.write("lgs.lp")
        
        # window
        model.addConstr(Q_solarRadWin[t] + Q_igRadWin +
                        sum(AWin)*alphaWin*(Tair[t]-Twin[t]) +
                        (window[t]-Twin[t])/RWin-min(sum(AWin),sum(AExt))*alphaRad[t]*(Twin[t]-Towi[t]) -
                        min(sum(AWin),AInt)*alphaRad[t]*(Twin[t]-Tiwi[t]) +
                        windowIndoorSurface[t] == 0,
                        name="window_balance_"+str(t))        
        
        # air
        model.addConstr(Qair[t] == Q_solarConv[t] + Q_igConv[t] -
                        sum(AWin)*alphaWin*(Tair[t]-Twin[t]) -
                        sum(AExt)*alphaExt*(Tair[t]-Towi[t]) -
                        AInt*alphaInt*(Tair[t]-Tiwi[t]) +
                        ventRate*rhoair*cair*(Tv[t]-Tair[t]) +
                        Q_HC[t],
                        name="air_mass_balance_"+str(t))
        
        # additional constraints
        model.addConstr(Qair[t] == (Tair[t]-Tair[t-1])/dt * rhoair*cair*Vair,
                        name="room_transfer_"+str(t))
        
        # bounds
    #    bounds = {}
    #    bounds["Tair_up"]  = 300
    #    bounds["Tair_low"] = 288
    #    model.addConstr(Tair[t] <= bounds["Tair_up"])
    #    model.addConstr(Tair[t] >= bounds["Tair_low"])
                        
    model.update()
#    model.write("lgs.lp")
                        
#%% get testing results
    
    # "improve" feasibility for model with numerical issues
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

    for t in range(timesteps):
        res["Tair"][t] = Tair[t].X-273.15
        res["Tiwi"][t] = Tiwi[t].X-273.15
        res["Towi"][t] = Towi[t].X-273.15
        res["Tow"][t]  = Tow[t].X-273.15
        res["Tiw"][t]  = Tiw[t].X-273.15
        
        res["Qair"][t] = Qair[t].X
        res["Qsrc"][t] = Q_solarConv[t]
        res["Qiwi"][t] = (Tiwi[t].X-Tiw[t].X)/RInt 
        res["Qowi"][t] = (Tow[t].X-Towi[t].X)/RExt
        res["Qowc"][t] = (Tair[t].X-Towi[t].X)*alphaExt*sum(AExt)
        res["Qiwc"][t] = (Tair[t].X-Tiwi[t].X)*alphaInt*AInt
        res["QHT"][t]  = (Towi[t].X-Tiwi[t].X)*alphaRad[t]*sum(AExt)
        res["QHC"][t]  = Q_HC[t].X
        res["Qv"][t]   = ventRate*rhoair*cair*Vair*(Tv[t]-Tair[t].X)
        res["Qig"][t]  = intGainsConv[t] + intGainsRad[t]
        res["Teq"][t]  = extWall[t]-273.15
        
        res["QTSW1"][t]= Q_solarRadExt[t]
        res["QTSL1"][t]= Q_igRadExt[t]
        res["QTSW2"][t]= Q_solarRadInt[t]
        res["QTSL2"][t]= Q_igRadInt[t]

#%% return model with added variables/constraints
    return res["Tair"], res["QHC"], res["Qiwi"], res["Qowi"]

def addConstraints(model, dt, testcase):
#%% 
    # get building parameters
    filename = "TEASER4_meine_Geo.mo"
    houseData = parser_buildings.parse_record(filename)
    
    # get weather data
    filename = "TRY2010_12_Jahr.dat"
    raw_inputs = {}
 
    beta, gamma = weather.get_betaGamma(houseData["orientationswallshorizontal"])
    (raw_inputs["solar_irrad_sky"],
    raw_inputs["solar_irrad_earth"],
    raw_inputs["temperature"], solarRad_in) = weather.get_weather(filename, beta, gamma)
    
    # calculate solar input for reducedOrderModel (ROM) and rad through sunblinds for equalAirTemp (EAT)
    solarRad_in_redOrdMod, sunblindsig = sunblinds.sunblinds(houseData, solarRad_in, testcase)

    # calculate equivalent air temperature on outer walls
    equalAirTemp, alphaRad = eqAirTemp.eqAirTempVDI([raw_inputs["temperature"],
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
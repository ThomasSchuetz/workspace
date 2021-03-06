#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 01 10:35:56 2015

"""

from __future__ import division

import gurobipy as gp
import numpy as np
import datetime

def compute(pass_house, marginals, eco, devs, demands, params):
    """
    Parameters
    ----------
    eco : dict
        - b :
        - crf :
        - prChange : 
        - q : 
        - rate : 
        - sub_CHP : 
        - t_calc : 
        - tax : 
        - gas, Boiler : 
        - gas, CHP : 
        - gas, c_meter : 
        - pr, el : 
        - sell, CHP : 
        - sell, PV : 
    devs : dict
        - bat :
        - boiler :
        - chp :
        - eh :
        - hp :
        - pv :
        - stc :
        - tes :
    demands : dict
        - design_heat_load : 
        - electricity : 
        - heat : 
        - solar_irrad : 
        - temperature : 
        - weights : 
    params : dict
        - dt : time step length (s)
        - maximum roof area (m2)
        - mip_gap : Solver setting (-)
        - time_limit : Solver setting (s)
        - days : 
        - time_steps : 
    """
    
    # Extract parameters
    dt = params["dt"]
    sigma = marginals["sigma"]
    
    # Define subsets
    heater = ("boiler", "chp", "eh", "hp")
    storage = ("bat", "tes")
    solar = ("pv", "stc")
    device = heater+storage+solar
    
    time_steps = range(params["time_steps"])
    days = range(params["days"])
    
    try:
        model = gp.Model("Design computation")
        
        # Define variables
        # Costs and Revenues
        c_inv = {dev: model.addVar(vtype="C", name="c_inv_"+dev)
                 for dev in devs.keys()}
        c_om = {dev: model.addVar(vtype="C", name="c_om_"+dev)
                  for dev in devs.keys()}
        c_dem = {dev: model.addVar(vtype="C", name="c_dem_"+dev)
                  for dev in ("boiler", "chp", "grid")}
        c_meter = model.addVar(vtype="C", name="c_meter")
        
        revenue = {dev: model.addVar(vtype="C", name="revenue_"+dev)
                  for dev in ("chp", "pv")}
        chp_subsidy = model.addVar(vtype="C", name="chp_subsidy")
        
        # SOC, power, heat and energy
        soc = {}
        power = {}
        heat = {}
        energy = {}
        soc_nom = {}
        power_nom = {}
        heat_nom = {}
        for d in days: # All days
            for t in time_steps: # All time steps of all days
                timetag = "_"+str(d)+"_"+str(t)
                for dev in storage: # All storage devices
                    soc[dev,d,t] = model.addVar(vtype="C",
                                                name="SOC_"+dev+"_"+timetag)
                
                for dev in (heater+solar):
                    power[dev,d,t] = model.addVar(vtype="C",
                                                  name="P_"+dev+"_"+timetag)
                    heat[dev,d,t] = model.addVar(vtype="C",
                                                 name="Q_"+dev+"_"+timetag)

                for dev in heater:
                    energy[dev,d,t] = model.addVar(vtype="C",
                                                   name="E_"+dev+"_"+timetag)
        
                for dev in heater:
                    heat_nom[dev,d,t] = model.addVar(vtype="C",
                                              name="Q_nom_"+dev+"_"+timetag)
                dev = "hp"
                power_nom[dev,d,t] = model.addVar(vtype="C",
                                           name="P_nom_"+dev+"_"+timetag)
        
        for dev in storage:
            soc_nom[dev] = model.addVar(vtype="C", name="SOC_nom_"+dev)

        # Storage initial SOC's
        soc_init = {}
        for dev in storage:
            for d in days:
                tag = dev + "_" + "_" + str(d)
                soc_init[dev,d] = model.addVar(vtype="C", name="SOC_init_"+tag)

        # Storage volume
        volume = model.addVar(vtype="C", name="Volume")

        # Storage charging and discharging
        p_ch = {}
        p_dch = {}
        for d in days:
            for t in time_steps:
                timetag = "_" + str(d) + "_" + str(t)

                p_ch[d,t] = model.addVar(vtype="C", name="P_ch"+timetag)
                p_dch[d,t] = model.addVar(vtype="C", name="P_dch"+timetag)
                    
        
        # Electricity imports, sold and self-used electricity
        p_imp = {}
        p_use = {}
        p_sell = {}
        for d in days:
            for t in time_steps:
                timetag = "_"+str(d)+"_"+str(t)
                
                p_imp[d,t] = model.addVar(vtype="C", name="P_imp_"+timetag)
                for dev in ("chp", "pv"):
                    p_use[dev,d,t] = model.addVar(vtype="C", 
                                                  name="P_use_"+dev+timetag)
                    p_sell[dev,d,t] = model.addVar(vtype="C", 
                                                   name="P_sell_"+dev+timetag)

        # Variables for determining the variable costs (based on capacity)
        capacity = {dev : model.addVar(vtype="C", name="Capacity_"+dev) 
                    for dev in devs.keys()}
        
        # Activation and purchase decision variables
        x = {}  # Purchase (all devices)
        y = {}  # Acitivation (heaters)
        for dev in devs.keys():
            x[dev] = model.addVar(vtype="B", name="x_"+dev, lb=0.0, ub=1.0)
        
        for d in days:
            for t in time_steps:
                timetag = "_"+str(d)+"_"+str(t)
                for dev in heater: # All heating devices
                    y[dev,d,t] = model.addVar(vtype="B", lb=0.0, ub=1.0,
                                              name="y_"+dev+"_"+timetag)
        # PV and STC areas
        area = {}
        for dev in solar:
            area[dev] = model.addVar(vtype="C", name="Area_"+dev)
        
        # Update
        model.update()

        # Define constraints
        # Objective
        model.setObjective(sum(c_inv[key] for key in c_inv.keys())
                         + sum(c_om[key] for key in c_om.keys())
                         + sum(c_dem[key] for key in c_dem.keys())
                         + c_meter
                         - chp_subsidy
                         - sum(revenue[key] for key in revenue.keys())
                         - sigma, 
                         gp.GRB.MINIMIZE)
#        model.setObjective(x["eh"]+x["hp"]+x["boiler"]+x["chp"]+x["stc"] 
#                           + sum(sum(y["boiler",d,t] for d in days) for t in time_steps)
#                           , gp.GRB.MAXIMIZE)

        # Determine device capacity
        # Boiler, CHP, EH, HP: Based on nominal heat output
        # PV and STC: Installed area
        # TES: Volume
        # BAT: Capacity
        # INV: PV Power (AC) / eta (resulting DC power)
        for d in days:
            for t in time_steps:
                for dev in heater:
                    model.addConstr(capacity[dev] >= heat_nom[dev,d,t],
                                    name="Capacity_"+dev+"_"+str(d)+"_"+str(t))
                dev = "inv"
                model.addConstr(capacity[dev] >= 
                                1 / devs[dev]["eta"] * power["pv",d,t], 
                                name="Capacity_"+dev)
        for dev in solar:
            model.addConstr(capacity[dev] == area[dev],
                            name="Capacity_"+dev)
        dev = "tes"
        model.addConstr(capacity[dev] == volume, name="Capacity_"+dev)
        dev = "bat"
        model.addConstr(capacity[dev] == soc_nom[dev], name="Capacity_"+dev)
        
        # Capacity is bounded:
        for dev in heater:
            model.addConstr(x[dev] * devs[dev]["Q_nom_min"] <= capacity[dev],
                            name="Capacity_min_"+dev)
            model.addConstr(x[dev] * devs[dev]["Q_nom_max"] >= capacity[dev],
                            name="Capacity_max_"+dev)
        
        # Economic constraints
        # Investment costs
        for dev in devs.keys():
            model.addConstr(c_inv[dev] == eco["crf"] * eco["tax"] *
                (1-devs[dev]["rval"]) *
                (x[dev] * devs[dev]["c_inv_fix"] 
                 + capacity[dev] * devs[dev]["c_inv_var"]),
                name="Investment_costs_"+dev)
        
        # Operation and maintenance
        for dev in c_om.keys():
            model.addConstr(c_om[dev] == eco["b"]["infl"] * eco["crf"] 
                * eco["tax"] * devs[dev]["c_om_rel"] * 
                (x[dev] * devs[dev]["c_inv_fix"] +
                 capacity[dev] * devs[dev]["c_inv_var"]),
                name="O&M_costs_"+dev)

        # Demand related costs (gas)
        for dev in ("boiler", "chp"):
            model.addConstr(c_dem[dev] == eco["b"]["gas"] * eco["crf"] *
                eco["gas",dev] * dt * 
                sum(demands["weights"][d] * 
                    sum(energy[dev,d,t] for t in time_steps)
                    for d in days),
                name="Demand_costs_"+dev)
                
        # Demand related costs (electricity)    # eco["b"]["el"] * eco["crf"]
        dev = "grid"
        model.addConstr(c_dem[dev] == sum(#demands["weights"][d] * 
#                sum((demands["electricity"][d,t] + power["hp",d,t] +power["eh",d,t] -power["chp",d,t] -power["pv",d,t]) * marginals["pi"][d,t] for t in time_steps)
                 sum((p_imp[d,t]) * marginals["pi"][d,t] for t in time_steps)
                 for d in days),
            name="Demand_costs_"+dev)
        
        # Revenues for selling electricity to the grid / neighborhood
        for dev in ("chp", "pv"):
            model.addConstr(revenue[dev] == sum(#demands["weights"][d] * 
                    sum(p_sell[dev,d,t] * marginals["pi"][d,t] for t in time_steps) 
                    for d in days),
                name="Feed_in_rev_"+dev)
                
#        # Revenues for selling electricity to the grid
#        for dev in ("chp", "pv"):
#            model.addConstr(revenue[dev] == eco["b"]["eex"] * eco["crf"] * dt *
#                sum(demands["weights"][d] *
#                    sum(p_sell[dev,d,t] * eco["sell",dev] 
#                    for t in time_steps)
#                for d in days),
#                name="Feed_in_rev_"+dev)

        # Metering costs
        for dev in ("boiler", "chp"):
            model.addConstr(c_meter >= eco["b"]["infl"] * eco["crf"] 
                            * eco["gas","c_meter"] * x[dev],
                            name="Metering_costs_"+dev)

        # CHP subsidies
        model.addConstr(chp_subsidy == eco["b"]["eex"] * eco["crf"] * dt *
                sum(demands["weights"][d] *
                    sum(power["chp",d,t] * eco["sub_chp"] 
                    for t in time_steps)
                for d in days),
                name="Subsidies_chp")       

        # Technical constraints
        # Maximum area constraints
        for dev in solar:
            # Minimum area for each device
            model.addConstr(area[dev] >= x[dev] * devs[dev]["area_min"],
                            name="Minimum_area_"+dev)
            # Maximum area for each device
            model.addConstr(area[dev] <= x[dev] * params["A_max"],
                            name="Maximum_area_"+dev)
                            
        # Maximum available area
        model.addConstr(sum(area[dev] for dev in solar) <= params["A_max"],
                        name="Maximum_total_area")
        
        # Necessity of an inverter:
        model.addConstr(x["inv"] == x["pv"])
        
        # Devices can be switched on only if they have been purchased
        for dev in heater:
            model.addConstr(sum( sum(y[dev,d,t] for t in time_steps)
                            for d in days) 
                  <= params["time_steps"] * params["days"] * x[dev],
                  name="Activation_"+dev)
        
        # Devices nominal values (heat_nom = y * capacity)
        for dev in heater:
            for d in days:
                for t in time_steps:
                    # Abbreviations
                    timetag = "_" + str(d) + "_" + str(t)
                    
                    q_nom_min = devs[dev]["Q_nom_min"]# * devs[dev]["mod_lvl"]
                    q_nom_max = devs[dev]["Q_nom_max"]
                    
                    model.addConstr(heat_nom[dev,d,t]
                        <= q_nom_max * y[dev,d,t],
                        name="Max_heat_1_"+dev+"_"+timetag)
                    model.addConstr(heat_nom[dev,d,t]
                        >= q_nom_min * y[dev,d,t],
                        name="Min_heat_1_"+dev+"_"+timetag)
                        
                    model.addConstr(capacity[dev] - heat_nom[dev,d,t]
                        <= q_nom_max * (x[dev] - y[dev,d,t]),
                        name="Max_heat_2_"+dev+"_"+timetag)
                    model.addConstr(capacity[dev] - heat_nom[dev,d,t]
                        >= q_nom_min * (x[dev] - y[dev,d,t]),
                        name="Min_heat_2_"+dev+"_"+timetag)

        # Compute nominal power consumption of HP:
        dev = "hp"
        for d in days:
            for t in time_steps:
                model.addConstr(power_nom[dev,d,t] * devs[dev]["cop_a2w35"]
                                == heat_nom[dev,d,t],
                                name="Power_nom_"+dev+"_"+str(d)+"_"+str(t))

        # Devices operation
        # Heat output between mod_lvl*Q_nom and Q_nom (P_nom for heat pumps)
        # Power and Energy directly result from Heat output
        for dev in heater:
            for d in days:
                for t in time_steps:
                    # Abbreviations
                    timetag = "_" + str(d) + "_" + str(t)
                    
                    mod_lvl = devs[dev]["mod_lvl"]
                    eta     = devs[dev]["eta"][d,t]
                    omega   = devs[dev]["omega"][d,t]
                    
                    if dev == "hp":
                        model.addConstr(power[dev,d,t] <= power_nom[dev,d,t],
                            name="Max_pow_operation_"+dev+"_"+timetag)
                        model.addConstr(power[dev,d,t] >= 
                            power_nom[dev,d,t] * mod_lvl,
                            name="Min_pow_operation_"+dev+"_"+timetag)
                    else:
                        model.addConstr(heat[dev,d,t] <= heat_nom[dev,d,t],
                            name="Max_heat_operation_"+dev+"_"+timetag)
                        model.addConstr(heat[dev,d,t] >= 
                            heat_nom[dev,d,t] * mod_lvl,
                            name="Min_heat_operation_"+dev+"_"+timetag)
                    
                    model.addConstr(power[dev,d,t] == 1/eta * heat[dev,d,t],
                          name="Power_equation_"+dev+"_"+timetag)
                            
                    model.addConstr(energy[dev,d,t] == 
                          1/omega * (heat[dev,d,t]+power[dev,d,t]),
                          name="Energy_equation_"+dev+"_"+timetag)

        # Solar components
        for dev in solar:
            for d in days:
                for t in time_steps:
                    timetag = "_" + str(d) + "_" + str(t)
                    eta_th = devs[dev]["eta_th"][d][t]
                    eta_el = devs["inv"]["eta"] * devs[dev]["eta_el"][d][t]
                    solar_irrad = demands["solar_irrad"][d][t]

                    model.addConstr(heat[dev,d,t] == 
                          eta_th * area[dev] * solar_irrad,
                          name="Solar_thermal_"+dev+timetag)
                    model.addConstr(power[dev,d,t] == 
                          eta_el * area[dev] * solar_irrad,
                          name="Solar_electrical_"+dev+timetag)
        
        # Storage equations
        # TES soc[dev,d,t] soc_init[dev,d] soc_nom[dev]
        dev = "tes"
        k_loss = devs[dev]["k_loss"]
        eta_ch = devs[dev]["eta_ch"]
        eta_dch = devs[dev]["eta_dch"]
        for d in days:
            for t in time_steps:
                if t == 0:
                    if np.max(demands["weights"]) == 1:
                        if d == 0:
                           soc_prev = soc_init[dev,d]
                        else:
                           soc_prev = soc[dev,d-1,params["time_steps"]-1]
                    else:
                        soc_prev = soc_init[dev,d]
                else:
                    soc_prev = soc[dev,d,t-1]
                
                timetag = "_" + str(d) + "_" + str(t)
                
                charge = eta_ch * sum(heat[dv,d,t] for dv in (heater+solar))
                discharge = 1 / eta_dch * demands["heat"][d,t]
                
                model.addConstr(soc[dev,d,t] == (1-k_loss) * soc_prev + 
                                dt * (charge - discharge),
                                name="Storage_bal_"+dev+timetag)
        
        # Nominal storage content (SOC)
        for dev in storage:
            for d in days:
                # Inits
                model.addConstr(soc_nom[dev] >= soc_init[dev,d], 
                                name="SOC_nom_inits_"+dev+"_"+str(d))
                for t in time_steps:
                    # Regular storage loads
                    model.addConstr(soc_nom[dev] >= soc[dev,d,t],
                                    name="SOC_nom_"+dev+"_"+str(d)+"_"+str(t))
        
        # Compute volume and ensure bounds
        dev = "tes"
        model.addConstr(soc_nom[dev] == params["rho_w"] * params["c_w"] * 
                        devs[dev]["dT_max"] / 3600000 * volume, 
                        name="Storage_Volume")
        model.addConstr(volume >= x[dev] * devs[dev]["volume_min"], 
                        name="Storage_Volume_min")
        model.addConstr(volume <= x[dev] * devs[dev]["volume_max"], 
                        name="Storage_Volume_max")

        dev = "bat"
        k_loss = devs[dev]["k_loss"]
        for d in days:
            for t in time_steps:
                if t == 0:
                    if np.max(demands["weights"]) == 1:
                        if d == 0:
                           soc_prev = soc_init[dev,d]
                        else:
                           soc_prev = soc[dev,d-1,params["time_steps"]-1]
                    else:
                        soc_prev = soc_init[dev,d]
                else:
                    soc_prev = soc[dev,d,t-1]

                timetag = "_" + str(d) + "_" + str(t)
    
                model.addConstr(soc[dev,d,t] == (1-k_loss) * soc_prev + 
                    dt * (devs[dev]["eta"] * p_ch[d,t] - 
                          1 / devs[dev]["eta"] * p_dch[d,t]),
                    name="Storage_balance_"+dev+timetag)
    
                model.addConstr(p_ch[d,t] <= x[dev] * devs[dev]["P_ch_fix"] + 
                                capacity[dev] * devs[dev]["P_ch_var"],
                                name="P_ch_max"+timetag)
    
                model.addConstr(p_dch[d,t] <= x[dev] * devs[dev]["P_dch_fix"] + 
                                capacity[dev] * devs[dev]["P_dch_var"],
                                name="P_dch_max"+timetag)
        
        # Ensure bounds
        dev = "bat"
        model.addConstr(soc_nom[dev] >= x[dev] * devs[dev]["cap_min"],
                        name="Battery_capacity_min")
        model.addConstr(soc_nom[dev] <= x[dev] * devs[dev]["cap_max"],
                        name="Battery_capacity_max")
        
        # SOC repetitions
        for dev in storage:
            for d in range(params["days"]):
                if np.max(demands["weights"]) > 1:
                    model.addConstr(soc[dev,d,params["time_steps"]-1] ==
                                    soc_init[dev,d],
                                    name="repetitions_"+dev+"_"+str(d))

        # Electricity balance (house)
        for d in days:
            for t in time_steps:
                model.addConstr(demands["electricity"][d,t]
                    + p_ch[d,t] - p_dch[d,t]
                    + power["hp",d,t] + power["eh",d,t]
                    - p_use["chp",d,t] - p_use["pv",d,t]
                    == p_imp[d,t],
                    name="Electricity_balance_"+str(d)+"_"+str(t))
        
        # Split CHP and PV generation into self-consumed and sold powers
        for dev in ("chp", "pv"):
            for d in days:
                for t in time_steps:
                    model.addConstr(p_sell[dev,d,t] + p_use[dev,d,t] 
                            == power[dev,d,t],
                            name="power=sell+use_"+dev+"_"+str(d)+"_"+str(t))
                            
        # Heat pump's operation depends on storage temperature
        for d in days:
            for t in time_steps:
                # Abbreviations
                dT_relative = devs["hp"]["dT_max"] / devs["tes"]["dT_max"]
                # Residual storage content
                resSC = (devs["tes"]["volume_max"] * devs["tes"]["dT_max"]
                         * params["rho_w"] * params["c_w"] * (1 - dT_relative)
                         / 3600000)
                
                model.addConstr(soc["tes",d,t] <= soc_nom["tes"] * dT_relative 
                      + (1 - y["hp",d,t]) * resSC,
                      name="Heat_pump_act_"+str(d)+"_"+str(t))
        
        # Design heat load has to be covered
        cop_relative = devs["hp"]["cop_design"] / devs["hp"]["cop_a2w35"]
        model.addConstr(sum(capacity[dev] for dev in ("boiler","chp","eh"))
                        + capacity["hp"] * cop_relative 
                        >= demands["design_heat_load"],
                        name="Design_heat_load")
       
        # FORCE CERTAIN EVENTS
#        model.addConstr(x["eh"]+x["hp"]+x["chp"]+x["stc"] == 0)
#        model.addConstr(x["pv"] == 1)
        
        # Set solver parameters
        model.Params.TimeLimit = params["time_limit"]
        model.Params.MIPGap = params["mip_gap"]
        model.Params.MIPFocus = 3
        
        # Execute calculation
        model.optimize()
#        model.computeIIS()
#        model.write("model.ilp")
        
        # Retrieve results
        res_x = {dev : x[dev].X  for dev in devs}
        res_y = {}
        res_energy = {}
        res_power = {}
        res_heat = {}
        res_soc = {}
        
        for dev in heater:
            res_y[dev] = {(d,t): y[dev,d,t].X
                                 for d in days
                                 for t in time_steps}
            res_energy[dev] = {(d,t) : energy[dev,d,t].X  
                                       for d in days
                                       for t in time_steps}
        for dev in (heater + solar):
            res_power[dev] = {(d,t) : power[dev,d,t].X  
                                      for d in days
                                      for t in time_steps}
        
            res_heat[dev] = {(d,t) : heat[dev,d,t].X  
                                     for d in days
                                     for t in time_steps}
        
        for dev in storage:
            res_soc[dev] = {(d,t): soc[dev,d,t].X 
                                   for d in days
                                   for t in time_steps}
    
        res_p_imp = {(d,t) : p_imp[d,t].X for d in days for t in time_steps}
        res_p_ch = {(d,t) : p_ch[d,t].X for d in days for t in time_steps}
        res_p_dch = {(d,t) : p_dch[d,t].X for d in days for t in time_steps}

        res_area = {dev : area[dev].X for dev in area.keys()}
        
        res_cap = {dev : capacity[dev].X for dev in capacity.keys()}
        
        res_volume = max(0.001, volume.X)
        res_temperature = {(d,t): res_soc["tes"][d,t]/(res_volume*1000*4180) for d in days for t in time_steps}
        
        obj = model.ObjVal
        
        
        res_c_inv = {dev: c_inv[dev].X for dev in c_inv.keys()}
        res_c_om  = {dev: c_om[dev].X for dev in c_om.keys()}
        res_c_dem = {dev: c_dem[dev].X for dev in c_dem.keys()}
        res_c_met = c_meter.X
        
        # res_rev = {dev: revenue[dev].X for dev in revenue.keys()}
        res_chp_sub = chp_subsidy.X
        
        res_soc_nom = {dev: soc_nom[dev].X for dev in storage}
        res_power_nom = {}
        res_heat_nom = {}
        for dev in heater:
            res_heat_nom[dev] = {(d,t): heat_nom[dev,d,t].X for d in days for t in time_steps}
        dev = "hp"
        res_power_nom[dev] = {(d,t): power_nom[dev,d,t].X for d in days for t in time_steps}
        
        res_soc_init = {}
        for dev in storage:
            res_soc_init[dev] = {d: soc_init[dev,d].X for d in days}
        
        res_p_use = {}
        res_p_sell = {}
        for dev in ("chp", "pv"):
            res_p_use[dev] = {(d,t): p_use[dev,d,t].X for d in days for t in time_steps}
            res_p_sell[dev] = {(d,t): p_sell[dev,d,t].X for d in days for t in time_steps}
        
        for dev in devs:
            print dev + ": " + str(res_x[dev])
        
        print
        for dev in devs:
            print "Capacity " + dev + ": " + str(res_cap[dev])
        
        print
        for dev in solar:
            print dev + ": " + str(res_area[dev])
                
        print
        print "Obj: " + str(model.ObjVal)
        
        # Compute "costs" of the proposal:
        cost = {}
        cost["inv"] = []
        cost["om"] = []
        cost["dem"] = []
        cost["grid"] = []
        cost["inv"] = res_c_inv["bat"]+res_c_inv["boiler"]+res_c_inv["chp"]+res_c_inv["eh"]+res_c_inv["hp"]+res_c_inv["tes"]+res_c_inv["pv"]+res_c_inv["stc"]+res_c_inv["inv"]
        cost["om"] = res_c_om["bat"]+res_c_om["boiler"]+res_c_om["chp"]+res_c_om["eh"]+res_c_om["hp"]+res_c_om["tes"]+res_c_om["pv"]+res_c_om["stc"]+res_c_om["inv"]
        cost["dem"] = res_c_dem["boiler"]+res_c_dem["chp"]
        cost["grid"] = res_c_dem["grid"]
        costs = (res_c_met - res_chp_sub + cost["inv"] + cost["om"] +cost["dem"])
            
        # Compute "proposals" of the devices:
        proposals = {}
        proposals["house"] = {}
        for dev in ("chp", "pv"):
            proposals[dev] = {(d,t): power[dev,d,t].X  for d in days for t in time_steps}
        proposals["hp"] = res_power["hp"]
        proposals["eh"] = res_power["eh"]
        for d in days:
            for t in time_steps:
                proposals["house"][d,t] = res_p_imp[d,t] - res_p_sell["chp"][d,t] - res_p_sell["pv"][d,t]
    
        #res_p_imp
        #res_p_imp = {}
        
        objVal = obj
        runtime = model.getAttr("Runtime")
#        mipgap = model.getAttr("MIPgap") 
        mipgap = 0
        datetime.datetime.now()   
#        model.computeIIS()
#        model.write("model.ilp")
#        print('\nConstraints:')        
#        for c in model.getConstrs():
#            if c.IISConstr:
#                print('%s' % c.constrName)
#        print('\nBounds:')
#        for v in model.getVars():
#            if v.IISLB > 0 :
#                print('Lower bound: %s' % v.VarName)
#            elif v.IISUB > 0:
#                print('Upper bound: %s' % v.VarName) 
        
        # Return results
        return (res_x, res_y, res_energy, res_power, res_heat, res_soc, res_p_imp,
                res_p_ch, res_p_dch, res_p_use, res_p_sell, res_area, res_cap,
                res_volume, res_temperature, obj, res_c_inv, res_c_om, res_c_dem,
                res_c_met, res_chp_sub, res_soc_nom, res_power_nom,
                res_heat_nom, res_soc_init, devs, costs, proposals, cost, objVal, runtime, mipgap)
        
    except gp.GurobiError:
        print("Error")
    
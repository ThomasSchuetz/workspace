# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 10:31:11 2017

@author: tsz
"""

from __future__ import division

import numpy as np
import pickle as pkl
import xlsxwriter

def read_demands(building):
    solar_irrad = np.loadtxt("raw_inputs/building_" + building + "/solar_rad_35deg.csv") / 1000
    sh = np.loadtxt("raw_inputs/building_" + building + "/space_heating.csv")
    dhw = np.loadtxt("raw_inputs/building_" + building + "/dhw.csv")
    elec = np.loadtxt("raw_inputs/building_" + building + "/electricity.csv")
    heat = sh + dhw
    
    heat_split = heat.reshape((24, 365), order="F")
    solar_irrad_split = solar_irrad.reshape((24, 365), order="F")
    elec_split = elec.reshape((24, 365), order="F")
    
    return (solar_irrad, elec, heat, heat_split, solar_irrad_split, elec_split)

def read_results(building, method, days):
    with open("results/res_" + building + "_" + method + "_" + str(days) + ".pkl", "rb") as f_in:
        opti_res = pkl.load(f_in)
        eco = pkl.load(f_in)
        devs = pkl.load(f_in)
        clustered = pkl.load(f_in)
        par = pkl.load(f_in)
        comp_time = pkl.load(f_in)

        # (x, y, z, energy, power,                  0-4
        # heat, soc, p_imp, ch, dch,                5-9
        # p_use, p_sell, area, c_inv, c_om,         10-14
        # c_dem, c_met, rev, soc_init, emissions,   15-19
        # ObjVal, Runtime, MIPGap)                  20-22
        # = opti_res

    return (opti_res, devs, clustered, comp_time)

def evaluate(opti_res, comp_time, clustered, heat_split, elec_split, solar_irrad_split, heat, elec, 
             solar_irrad, devs, method):
    obj_val = opti_res[20]
    calc_time = comp_time["time optimization"]
    
    sse = 0
    sse_scaled = 0
    
    z = clustered["z"]
    
    if method in ("center", "medoid"):
        k = 0
        for i in range(365):
            if np.sum(z[i,:]) >= 1:
                for j in range(365):
                    if round(z[i,j]) >= 0.5:
                        sse += np.sum((heat_split[:,j] - clustered["heat"][k]) ** 2)
                        sse += np.sum((elec_split[:,j] - clustered["electricity"][k]) ** 2)
                        sse += np.sum((solar_irrad_split[:,j] - clustered["solar_irrad"][k]) ** 2)
                        
                        sse_scaled += np.sum(((heat_split[:,j] - clustered["heat"][k]) / np.max(heat)) ** 2)
                        sse_scaled += np.sum(((elec_split[:,j] - clustered["electricity"][k]) / np.max(elec)) ** 2)
                        sse_scaled += np.sum(((solar_irrad_split[:,j] - clustered["solar_irrad"][k]) / np.max(solar_irrad)) ** 2)
        
                k += 1
    
    elif method in ("mean", "median"):
        for j in range(365):
            sse += np.sum((heat_split[:,j] - clustered["heat"][z[j]]) ** 2)
            sse += np.sum((elec_split[:,j] - clustered["electricity"][z[j]]) ** 2)
            sse += np.sum((solar_irrad_split[:,j] - clustered["solar_irrad"][z[j]]) ** 2)
            
            sse_scaled += np.sum(((heat_split[:,j] - clustered["heat"][z[j]]) / np.max(heat)) ** 2)
            sse_scaled += np.sum(((elec_split[:,j] - clustered["electricity"][z[j]]) / np.max(elec)) ** 2)
            sse_scaled += np.sum(((solar_irrad_split[:,j] - clustered["solar_irrad"][z[j]]) / np.max(solar_irrad)) ** 2)
    
    elif method == "monthly":
        k = 0
        limits = np.cumsum(z)
        for j in range(365):
            if j == limits[k]:
                k += 1
            sse += np.sum((heat_split[:,j] - clustered["heat"][k,:]) ** 2)
            sse += np.sum((elec_split[:,j] - clustered["electricity"][k,:]) ** 2)
            sse += np.sum((solar_irrad_split[:,j] - clustered["solar_irrad"][k,:]) ** 2)
            
            sse_scaled += np.sum(((heat_split[:,j] - clustered["heat"][k,:]) / np.max(heat)) ** 2)
            sse_scaled += np.sum(((elec_split[:,j] - clustered["electricity"][k,:]) / np.max(elec)) ** 2)
            sse_scaled += np.sum(((solar_irrad_split[:,j] - clustered["solar_irrad"][k,:]) / np.max(solar_irrad)) ** 2)

    elif method == "seasonal":
        days_month = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
        days_month_cumsum = np.cumsum(days_month)
        
        def split_seasonal(orig_split, dmc):
            temp_split = np.split(orig_split, [dmc[1], dmc[4], dmc[7], dmc[10]], axis=1)
            return [np.hstack((temp_split[0], temp_split[4])), # Winter
                               temp_split[2], # Summer
                               np.hstack((temp_split[1], temp_split[3]))] # Transition

        heat_seasonal = split_seasonal(heat_split, days_month_cumsum)
        elec_seasonal = split_seasonal(elec_split, days_month_cumsum)
        solar_seasonal = split_seasonal(solar_irrad_split, days_month_cumsum)

        for k in range(len(z)):
            for j in range(z[k]):
                sse += np.sum((heat_seasonal[k][:,j] - clustered["heat"][k,:]) ** 2)
                sse += np.sum((elec_seasonal[k][:,j] - clustered["electricity"][k,:]) ** 2)
                sse += np.sum((solar_seasonal[k][:,j] - clustered["solar_irrad"][k,:]) ** 2)
                
                sse_scaled += np.sum(((heat_seasonal[k][:,j] - clustered["heat"][k,:]) / np.max(heat)) ** 2)
                sse_scaled += np.sum(((elec_seasonal[k][:,j] - clustered["electricity"][k,:]) / np.max(elec)) ** 2)
                sse_scaled += np.sum(((solar_seasonal[k][:,j] - clustered["solar_irrad"][k,:]) / np.max(solar_irrad)) ** 2)
    
    r_x = opti_res[0]
    r_z = opti_res[2]
    
    selection = {}
    caps = {}
    module_count = {}
    
    for dev in r_x.keys():
        selection[dev] = 0
        caps[dev] = 0
        
        for n in r_x[dev].keys():
            if r_x[dev][n] >= 0.5:
                selection[dev] = n
                if dev == "bat":
                    caps[dev] = devs[dev][n]["cap"]
                elif dev == "tes":
                    caps[dev] = devs[dev][n]["volume"]
                elif dev in ("boiler", "chp"):
                    caps[dev] = max(devs[dev][n]["Q_heat"])
                elif dev in ("eh", "hp"):
                    caps[dev] = devs[dev][n]["Q_nom"]
                elif dev == "inv":
                    caps[dev] = devs[dev][n]["P_nom_DC"]
                elif dev in ("pv", "stc"):
                    caps[dev] = devs[dev][n]["area"] * r_z[dev][n]
    
    for dev in r_z.keys():
        if selection[dev] == 0:
            module_count[dev] = 0
        else:
            module_count[dev] = r_z[dev][selection[dev]]    

    # C design --> investments, o&m, metering
    # C demand --> demands
    # Revenue --> rev
    # Total --> design + demand - revenue
    c_inv = opti_res[13]
    c_om = opti_res[14]
    c_dem = opti_res[15]
    c_met = opti_res[16]
    c_rev = opti_res[17]
    
    c_design = sum(c_inv.values()) + sum(c_om.values()) + c_met
    c_demand = sum(c_dem.values())
    c_revenue = sum(c_rev.values())

    c_total = c_design + c_demand - c_revenue #opti_res[20]
    
    costs = {"design": c_design, "demand": c_demand, "revenue": c_revenue, "total": c_total}

    return {"calc_time": calc_time, "sse":sse, "sse_scaled":sse_scaled,
            "selection": selection, "caps": caps, "modules": module_count,
            "obj_val": obj_val, "costs": costs}

# Layers: 1. House, 2. Method, 3. Number days
# Content: calc_time, sse, energy-sys, objective value
results = {}

for building in ("03", "07", "10"):
    results[building] = {}

    (solar_irrad, elec, heat, heat_split, solar_irrad_split, elec_split) = read_demands(building)
    
    for method in ("mean", "median", "center", "medoid"):
        results[building][method] = {}
        
        for d in range(1,26):
            (opti_res, devs, clustered, comp_time) = read_results(building, method, d)
            results[building][method][d] = evaluate(opti_res, comp_time, clustered, 
                    heat_split, elec_split, solar_irrad_split, heat, elec, solar_irrad, devs, method)
    
    (opti_res, devs, clustered, comp_time) = read_results(building, "monthly", 12)
    results[building]["monthly"] = evaluate(opti_res, comp_time, clustered, 
                    heat_split, elec_split, solar_irrad_split, heat, elec, solar_irrad, devs, "monthly")
                    
    (opti_res, devs, clustered, comp_time) = read_results(building, "seasonal", 3)
    results[building]["seasonal"] = evaluate(opti_res, comp_time, clustered, 
                    heat_split, elec_split, solar_irrad_split, heat, elec, solar_irrad, devs, "seasonal")

with open("results/summary.pkl", "wb") as fout:
    pkl.dump(results, fout, pkl.HIGHEST_PROTOCOL)

def write_col(method, d, col_counter, res):
    sheet.write(0,col_counter, method + " " + str(d))
    sheet.write(1,col_counter, res["obj_val"])
    sheet.write(3,col_counter, res["caps"]["boiler"])
    sheet.write(4,col_counter, res["caps"]["chp"])
    sheet.write(5,col_counter, res["caps"]["eh"])
    sheet.write(6,col_counter, res["caps"]["hp"])
    sheet.write(8,col_counter, res["caps"]["bat"])
    sheet.write(9,col_counter, res["caps"]["tes"])
    sheet.write(11,col_counter, res["caps"]["pv"])
    sheet.write(12,col_counter, res["caps"]["stc"])
    sheet.write(14,col_counter, res["sse_scaled"])
    sheet.write(15,col_counter, res["calc_time"])
    
    sheet.write(17,col_counter, res["costs"]["design"])
    sheet.write(18,col_counter, res["costs"]["demand"])
    sheet.write(19,col_counter, res["costs"]["revenue"])
    sheet.write(20,col_counter, res["costs"]["total"])

book = xlsxwriter.Workbook("results/results.xlsx")
for building in ("03", "07", "10"):
    sheet = book.add_worksheet(name="building_"+building)
    
    # Create outline
    sheet.write(1, 0, "Objective")
    
    sheet.write(3, 0, "BOI")
    sheet.write(4, 0, "CHP")
    sheet.write(5, 0, "EH")
    sheet.write(6, 0, "HP")

    sheet.write(8, 0, "BAT")
    sheet.write(9, 0, "TES")

    sheet.write(11, 0, "PV")
    sheet.write(12, 0, "STC")
    
    sheet.write(14, 0, "SSE scaled")
    sheet.write(15, 0, "Calc time")
    
    sheet.write(17, 0, "C design")
    sheet.write(18, 0, "C demand")
    sheet.write(19, 0, "C revenue")
    sheet.write(20, 0, "C total")

    col_counter = 1
    for method in ("mean", "median", "center", "medoid"):
        for d in range(1,26):
            write_col(method, d, col_counter, results[building][method][d])
            col_counter += 1
    
    write_col("monthly", 12, col_counter, results[building]["monthly"])
    write_col("seasonal", 3, col_counter+1, results[building]["seasonal"])
    
book.close()

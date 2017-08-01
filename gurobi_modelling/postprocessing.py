#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 02 10:49:10 2015

@author: tsz
"""

from __future__ import division


import pickle
import numpy as np
import os

def read_single_entry(filename):

    result = {}    
    
    # Load results
    with open(filename, "rb") as f_out:
        opti_res = pickle.load(f_out)
        eco = pickle.load(f_out)
        devs = pickle.load(f_out)
        clustered = pickle.load(f_out)
        par = pickle.load(f_out)
        comp_time = pickle.load(f_out)
    
    (x, y, energy, power, heat, soc, p_imp, p_ch, p_dch, p_use, p_sell, area, cap,
    volume, c_inv, c_om, c_dem, c_met, rev, chp_sub, soc_nom, power_nom, heat_nom,
    soc_init, p_use, p_sell, objVal, runtime, mipgap) = opti_res
    
    result["capacity"] = cap
    result["costs_invest"] = c_inv
    result["costs_om"] = c_om
    result["costs_demand"] = c_dem
    result["costs_metering"] = c_met
    result["rev_feed_in"] = rev
    result["rev_chp_sub"] = chp_sub
    result["objective"] = objVal
    result["mipgap"] = mipgap
    
    comp_time["time clustering"] = clustered["time used"]
    result["comp_time"] = comp_time

    result["total_invest"] = sum(c_inv.values()) + sum(c_om.values()) + c_met
    result["total_demand"] = sum(c_dem.values())
    result["total_revenue"] = sum(rev.values()) + chp_sub
    
    days = range(len(clustered["weights"]))
    result["ldc"] = [np.flipud(np.sort(np.hstack(
                        [np.tile(clustered[j][i,:], clustered["weights"][i]) 
                         for i in days]))) 
                     for j in ("electricity", "heat", "solar_irrad")]
    
    result["time_optimization"] = comp_time["time optimization"]
    result["time_clustering"] = comp_time["time clustering"]
    
    ###################################################
    # Determine structure
    result["structure"] = {dev: x[dev] for dev in x.keys()}
    
    # Characteristics
    # Total heat and electricity demand (building, not energy system)
    demand_electrical_building = 0
    demand_thermal_building = 0
    
    # Solar fractions
    
    # Total solar generation
    solar_thermal_total = 0
    solar_electrical_total = 0
    for d in days:
        demand_thermal_building += np.sum(clustered["heat"][d,:] * clustered["weights"][d])
        demand_electrical_building += np.sum(clustered["electricity"][d,:] * clustered["weights"][d])
        
        solar_thermal_total += np.sum(heat["stc"][d,:] * clustered["weights"][d])
        solar_electrical_total += np.sum(power["pv"][d,:] * clustered["weights"][d])
    
    # Transform to kWh/a
    to_kWh = lambda j: j * par["dt"] / (3600 * 1000)
    demand_electrical_building = to_kWh(demand_electrical_building)
    demand_thermal_building = to_kWh(demand_thermal_building)
    solar_thermal_total = to_kWh(solar_thermal_total)
    solar_electrical_total = to_kWh(solar_electrical_total)
    
    solar_fraction_electrical_building = solar_electrical_total / demand_electrical_building
    solar_fraction_thermal_building = solar_thermal_total / demand_thermal_building
#    result["solar_fraction_electrical_building"] = solar_fraction_electrical_building
#    result["solar_fraction_thermal_building"] = solar_fraction_thermal_building
    
    # Operating and full load hours:
    full_load_hours = {}
    operating_hours = {}
    for dev in ("chp", "boiler", "hp", "eh"):
        full_load_hours[dev] = 0
        operating_hours[dev] = 0
        if cap[dev] > 0:
            for d in days:
                operating_hours[dev] += np.sum(clustered["weights"][d] * np.sum(heat[dev][d,:]>0))
                full_load_hours[dev] += np.sum(clustered["weights"][d] * np.sum(heat[dev][d,:]/cap[dev]))
        
        full_load_hours[dev] *= par["dt"] / 3600
        operating_hours[dev] *= par["dt"] / 3600
    result["full_load_hours"] = full_load_hours
    result["operating_hours"] = operating_hours
    
    # Investments (not annualized)
    costs_total_inv = {}
    for dev in devs.keys():
        costs_total_inv[dev] = x[dev] * devs[dev]["c_inv_fix"] + cap[dev] * devs[dev]["c_inv_var"]
    
    # Electricity generation/consumption/grid
    electricity_generation = {}
    electricity_consumption = {}
    electricity_consumption["grid"] = 0
    electricity_generation["chp", "total"] = 0
    electricity_generation["chp", "sell"] = 0
    electricity_generation["chp", "use"] = 0
    electricity_generation["pv", "total"] = 0
    electricity_generation["pv", "sell"] = 0
    electricity_generation["pv", "use"] = 0
    for d in days:
        electricity_consumption["grid"] += np.sum(clustered["weights"][d] * p_imp[d,:])
        for dev in ("chp", "pv"):
            electricity_generation[dev, "total"] += np.sum(clustered["weights"][d] * (p_sell[dev][d,:] + p_use[dev][d,:]))
            electricity_generation[dev, "sell"] += np.sum(clustered["weights"][d] * p_sell[dev][d,:])
            electricity_generation[dev, "use"] += np.sum(clustered["weights"][d] * p_use[dev][d,:])
    
    electricity_consumption["grid"] = to_kWh(electricity_consumption["grid"])
    for dev in ("chp", "pv"):
        electricity_generation[dev, "total"] = to_kWh(electricity_generation[dev, "total"])
        electricity_generation[dev, "sell"] = to_kWh(electricity_generation[dev, "sell"])
        electricity_generation[dev, "use"] = to_kWh(electricity_generation[dev, "use"])
        
    for dev in ("hp", "eh"):
        electricity_consumption[dev] = 0
        for d in days:
            electricity_consumption[dev] += np.sum(clustered["weights"][d] * power[dev][d,:])
    
        electricity_consumption[dev] = to_kWh(electricity_consumption[dev])
    
    # Heat generation / consumption
    heat_total_consumption = np.sum([clustered["weights"][d] * clustered["heat"][d,:] for d in days])
    heat_total_consumption = to_kWh(heat_total_consumption)
    
    heat_total_generation = {}
    for dev in ("boiler", "chp", "hp", "eh", "stc"):
        heat_total_generation[dev] = 0
        for d in days:
            heat_total_generation[dev] += np.sum(clustered["weights"][d] * heat[dev][d,:])
        
        heat_total_generation[dev] = to_kWh(heat_total_generation[dev])
    result["heat_total_generation"] = heat_total_generation
    
    # Compute average and max. storage usage
    soc_average = {dev: np.mean(soc[dev]) for dev in ("tes", "bat")}
#    result["soc_average"] = soc_average
    soc_max = {dev: np.max(soc[dev]) for dev in ("tes", "bat")}
#    result["soc_max"] = soc_max
    
    # Compute storage losses
    sto_loss_tes = (sum([heat_total_generation[dev] for dev in heat_total_generation.keys()]) - heat_total_consumption)
    sto_loss_bat = (sum(electricity_generation[dev, "use"] for dev in ("chp", "pv"))
                    + electricity_consumption["grid"]
                    - sum(electricity_consumption[dev] for dev in ("eh", "hp"))
                    - demand_electrical_building)
#    result["sto_loss_tes"] = sto_loss_tes
#    result["sto_loss_bat"] = sto_loss_bat
    
    heat_fractions = {dev: heat_total_generation[dev] / (heat_total_consumption + sto_loss_tes)
                      for dev in heat_total_generation.keys()}
    result["heat_fractions"] = heat_fractions
    
    power_fractions = {}
    power_fractions["chp"] = p_use["chp"].sum() / (p_use["chp"].sum() + p_use["pv"].sum() + p_imp.sum())
    power_fractions["pv"] = p_use["pv"].sum() / (p_use["chp"].sum() + p_use["pv"].sum() + p_imp.sum())
    power_fractions["grid"] = p_imp.sum() / (p_use["chp"].sum() + p_use["pv"].sum() + p_imp.sum())
    result["power_fractions"] = power_fractions
    
    return result


def cleanup(dic, house_type, numbers, clus_types = ("mean", "center", "median", "medoid")):
    for number in numbers:
        for clus_type in clus_types:
            key = "res_"+house_type+"_"+clus_type+"_"+str(number)+".pkl"
            dic.pop(key)
    return dic

results_all = {}
results = {"AB":{}, "MFH":{}, "SFH":{}}

filenames = [file for file in os.listdir("results") if file.endswith(".pkl")]

for filename in filenames:
    
    fname = "results/" + filename

    results_all[filename] = read_single_entry(fname)

for k_res in results.keys():
    for k_all in results_all.keys():
        if k_res in k_all:
            results[k_res][k_all] = results_all[k_all]

#list_del = [1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17,18,19]
list_del = []
res_AB = cleanup(results["AB"], "AB", list_del)
res_MFH = cleanup(results["MFH"], "MFH", list_del)
res_SFH = cleanup(results["SFH"], "SFH", list_del)
eco_SFH = {key: results["SFH"][key]["objective"] for key in results["SFH"].keys()}
eco_MFH = {key: results["MFH"][key]["objective"] for key in results["MFH"].keys()}
eco_AB = {key: results["AB"][key]["objective"] for key in results["AB"].keys()}

def write_line(arg):
    if type(arg) == type({}):
        return (str(arg).replace("{","").replace("}","").replace(": -0.0", ": 0")
                        .replace(": 0.0", ": 0"). replace(": 1.0", ": 1").replace("u'", "").replace("'", ""))
    else:
        return arg

def write_sheet(workbook, name, key):
    ws = workbook.add_worksheet(name)
    headings = ["SFH", "MFH", "AB"]
    for i in range(len(headings)):
        ws.write(0, i+1, headings[i])
    
    categories = ["Reference", "365 uncon",
                  "Center 6", "Center 12", "Center 20",
                  "Medoid 6", "Medoid 12", "Medoid 20",
                  "Mean 6", "Mean 12", "Mean 20",
                  "Median 6", "Median 12", "Median 20",
                  "Monthly"]
    kSFH = ['orig_res_SFH_reference_365.pkl', 'res_SFH_reference_365_not_connected.pkl',
            'recalc_SFH_center_6.pkl', 'recalc_SFH_center_12.pkl', 'recalc_SFH_center_20.pkl', 
            'recalc_SFH_medoid_6.pkl', 'recalc_SFH_medoid_12.pkl', 'recalc_SFH_medoid_20.pkl', 
            'recalc_SFH_mean_6.pkl', 'recalc_SFH_mean_12.pkl', 'recalc_SFH_mean_20.pkl',
            'recalc_SFH_median_6.pkl', 'recalc_SFH_median_12.pkl', 'recalc_SFH_median_20.pkl',
            'recalc_SFH_monthly_12.pkl']
#    kMFH = ['orig_res_MFH_reference_365.pkl', 'res_MFH_reference_365_not_connected.pkl',
#            'res_MFH_center_6.pkl', 'res_MFH_center_12.pkl', 'res_MFH_center_20.pkl', 
#            'res_MFH_medoid_6.pkl', 'res_MFH_medoid_12.pkl', 'res_MFH_medoid_20.pkl', 
#            'res_MFH_mean_6.pkl', 'res_MFH_mean_12.pkl', 'res_MFH_mean_20.pkl',
#            'res_MFH_median_6.pkl', 'res_MFH_median_12.pkl', 'res_MFH_median_20.pkl',
#            'res_MFH_monthly_12.pkl']
    kAB = ['orig_res_AB_reference_365.pkl', 'res_AB_reference_365_not_connected.pkl',
           'recalc_AB_center_6.pkl', 'recalc_AB_center_12.pkl', 'recalc_AB_center_20.pkl', 
           'recalc_AB_medoid_6.pkl', 'recalc_AB_medoid_12.pkl', 'recalc_AB_medoid_20.pkl', 
           'recalc_AB_mean_6.pkl', 'recalc_AB_mean_12.pkl', 'recalc_AB_mean_20.pkl',
           'recalc_AB_median_6.pkl', 'recalc_AB_median_12.pkl', 'recalc_AB_median_20.pkl',
           'recalc_AB_monthly_12.pkl']
            
    for i in range(len(categories)):
        ws.write(i+1, 0, categories[i])
        if type(key) == type(""):
            ws.write(i+1, 1, write_line(res_SFH[kSFH[i]][key]))
#            ws.write(i+1, 2, write_line(res_MFH[kMFH[i]][key]))
            ws.write(i+1, 3, write_line(res_AB[kAB[i]][key]))
        else:
            ws.write(i+1, 1, write_line(res_SFH[kSFH[i]][key[0]][key[1]]))
#            ws.write(i+1, 2, write_line(res_MFH[kMFH[i]][key[0]][key[1]]))
            ws.write(i+1, 3, write_line(res_AB[kAB[i]][key[0]][key[1]]))

import xlsxwriter
wb = xlsxwriter.Workbook("results.xlsx")
write_sheet(wb, "total costs", "objective")
write_sheet(wb, "time optimization", ["comp_time", "time optimization"])
write_sheet(wb, "time clustering", ["comp_time", "time clustering"])
write_sheet(wb, "capacities", "capacity")
write_sheet(wb, "structure", "structure")
write_sheet(wb, "investments", "costs_invest")
write_sheet(wb, "full load hours", "full_load_hours")
write_sheet(wb, "heat fractions", "heat_fractions")
write_sheet(wb, "power fractions", "power_fractions")
wb.close()

import post_process_helpers as pph
days = 25
obj_SFH = pph.get_results_ordered(res_SFH, days=days, key="objective", prefix="recalc")
#obj_MFH = pph.get_results_ordered(res_MFH, days=days, key="objective")
obj_AB = pph.get_results_ordered(res_AB, days=days, key="objective", prefix="recalc")

#time_SFH = pph.get_results_ordered(res_SFH, days=days, key="time_optimization", prefix="recalc")
#time_MFH = pph.get_results_ordered(res_MFH, days=days, key="time_optimization")
#time_AB = pph.get_results_ordered(res_AB, days=days, key="time_optimization", prefix="recalc")

#timecls_SFH = pph.get_results_ordered(res_SFH, days=days, key="time_clustering")
#timecls_MFH = pph.get_results_ordered(res_MFH, days=days, key="time_clustering")
#timecls_AB = pph.get_results_ordered(res_AB, days=days, key="time_clustering")

#tot_dem_SFH = pph.get_results_ordered(res_SFH, days=days, key="total_demand")
#tot_dem_MFH = pph.get_results_ordered(res_MFH, days=days, key="total_demand")
#tot_dem_AB = pph.get_results_ordered(res_AB, days=days, key="total_demand")

#tot_inv_SFH = pph.get_results_ordered(res_SFH, days=days, key="total_invest")
#tot_inv_MFH = pph.get_results_ordered(res_MFH, days=days, key="total_invest")
#tot_inv_AB = pph.get_results_ordered(res_AB, days=days, key="total_invest")

#tot_rev_SFH = pph.get_results_ordered(res_SFH, days=days, key="total_revenue")
#tot_rev_MFH = pph.get_results_ordered(res_MFH, days=days, key="total_revenue")
#tot_rev_AB = pph.get_results_ordered(res_AB, days=days, key="total_revenue")


#pph.plot_ordered(obj_SFH, ytitle="Deviation SFH annual costs", filename="sfh_annual_costs.png", pos_legend=1)
#pph.plot_ordered(obj_MFH, ytitle="Deviation MFH annual costs", filename="mfh_annual_costs.png", pos_legend=1)
#pph.plot_ordered(obj_AB, ytitle="Deviation AB annual costs", filename="ab_annual_costs.png", pos_legend=1)
#
#pph.plot_ordered(time_SFH, ytitle="Time optimization in s",
#                 filename="sfh_time_optimization.png", relative=False, pos_legend=9)
#pph.plot_ordered(time_MFH, ytitle="Time optimization in s",
#                 filename="mfh_time_optimization.png", relative=False, pos_legend=2)
#pph.plot_ordered(time_AB, ytitle="Time optimization in s",
#                 filename="ab_time_optimization.png", relative=False, pos_legend=2)
#pph.plot_times(time_SFH, ytitle="Time optimization SFH in s",
#                 filename="sfh_time_optimization.png", pos_legend=2, threshold=600)
#pph.plot_times(time_MFH, ytitle="Time optimization MFH in s",
#                 filename="mfh_time_optimization.png", pos_legend=2, threshold=200)
#pph.plot_times(time_AB, ytitle="Time optimization AB in s",
#                 filename="ab_time_optimization.png", pos_legend=2, threshold=600)

#pph.plot_ordered(tot_dem_SFH, ytitle="Demand related costs in EUR",
#                 filename="sfh_tot_dem_optimization.png", relative=True, pos_legend=1)
#pph.plot_ordered(tot_dem_MFH, ytitle="Demand related costs in EUR",
#                 filename="mfh_tot_dem_optimization.png", relative=True, pos_legend=1)
#pph.plot_ordered(tot_dem_AB, ytitle="Demand related costs in EUR",
#                 filename="ab_tot_dem_optimization.png", relative=True, pos_legend=1)

#pph.plot_ordered(tot_inv_SFH, ytitle="Investments in EUR",
#                 filename="sfh_tot_inv_optimization.png", relative=True, pos_legend=1)
#pph.plot_ordered(tot_inv_MFH, ytitle="Investments in EUR",
#                 filename="mfh_tot_inv_optimization.png", relative=True, pos_legend=1)
#pph.plot_ordered(tot_inv_AB, ytitle="Investments in EUR",
#                 filename="ab_tot_inv_optimization.png", relative=True, pos_legend=1)

#pph.plot_ordered(tot_rev_SFH, ytitle="Revenues in EUR",
#                 filename="sfh_tot_rev_optimization.png", relative=True, pos_legend=1)
#pph.plot_ordered(tot_rev_MFH, ytitle="Revenues in EUR",
#                 filename="mfh_tot_rev_optimization.png", relative=True, pos_legend=4)
#pph.plot_ordered(tot_rev_AB, ytitle="Revenues in EUR",
#                 filename="ab_tot_rev_optimization.png", relative=True, pos_legend=1)
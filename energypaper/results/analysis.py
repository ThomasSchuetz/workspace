#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 08:29:09 2016

@author: tsz
"""

from __future__ import division

import pickle as pkl
import numpy as np

with open("inputs.pkl", "rb") as fin:
    pkl.load(fin)
    pkl.load(fin)
    eco = pkl.load(fin)
    devs = pkl.load(fin)
    clustered = pkl.load(fin)
    params = pkl.load(fin)

def evaluate(results, devs, demands):
    """
    Capacities, coverages, run times, full load hours, partial costs, costs, emissions, gap, calc times
    """
    res = {}
    res["gap"] = results["gap"]
    res["run time"] = results["runtime"]
    res["emissions"] = results["emissions"]
    res["c_total"] = results["c_total"]
    res["c_inv"] = results["c_inv"]
    res["c_om"] = results["c_om"]
    res["c_dem"] = results["c_dem"]
    res["c_fix"] = results["c_fix"]
    res["rev"] = results["rev"]
    res["sub"] = results["sub"]
    
    # Capacities
    capacities = {dev: 0 for dev in devs.keys()}
    for dev in devs.keys():
        for n in devs[dev].keys():
            if round(results["x"][dev][n]) == 1:
                if dev in ("pv", "stc"):
                    capacities[dev] = results["z"][dev][n] * devs[dev][n]["area"]
                elif dev == "inv":
                    capacities[dev] = devs[dev][n]["P_nom_DC"]
                elif dev == "tes":
                    capacities[dev] = devs[dev][n]["volume"]
                elif dev == "bat":
                    capacities[dev] = devs[dev][n]["cap"]
                elif dev == "eh":
                    capacities[dev] = devs[dev][n]["Q_nom"]
                elif dev in ("boiler", "chp"):
                    capacities[dev] = max(devs[dev][n]["Q_heat"])
                elif dev == "hp":
                    capacities[dev] = devs[dev][n]["c_inv"]

    res["cap"] = capacities
    
    (num_days, num_time_steps) = np.shape(demands["heat"])
    days = range(num_days)
    time_steps = range(num_time_steps)
    
    
    # Coverages (heat)
    # Coverages (electrical)
    heat_total = {}
    power_total = {}
    for dev in ("eh", "chp", "hp", "boiler", "stc", "pv"):
        heat_total[dev] = 0
        power_total[dev] = 0
        for n in devs[dev].keys():
            if round(results["x"][dev][n]) == 1:
                for d in days:
                    for t in time_steps:
                        heat_total[dev] += demands["weights"][d] * results["heat"][dev][n][d,t]
                        power_total[dev] += demands["weights"][d] * results["power"][dev][n][d,t]
    
    res["heat_total"] = heat_total
    res["power_total"] = power_total

    sum_heat = sum([heat_total[key] for key in heat_total.keys()])
    heat_partials = {key: heat_total[key]/sum_heat for key in heat_total.keys()}
    res["heat_partials"] = heat_partials
    
    grid_total = {}    
    for key in ("house", "hp"):
        grid_total[key] = 0
        for d in days:
            for t in time_steps:
                grid_total[key] += demands["weights"][d] * results["p_grid"][key][d,t]
    
    res["grid_total"] = grid_total
    
    use_total = {}
    sell_total = {}
    hp_total = {}
    for dev in ("chp", "bat", "pv"):
        use_total[dev] = 0
        sell_total[dev] = 0
        hp_total[dev] = 0
        if dev == "chp":
            for n in devs[dev].keys():
                if round(results["x"][dev][n]) == 1:
                    for d in days:
                        for t in time_steps:
                            use_total[dev] += demands["weights"][d] * results["p_use"][dev][n-1,d,t]
                            sell_total[dev] += demands["weights"][d] * results["p_sell"][dev][n-1,d,t]
                            hp_total[dev] += demands["weights"][d] * results["p_hp"][dev][n-1,d,t]
                            
        else:
            for d in days:
                for t in time_steps:
                    use_total[dev] += demands["weights"][d] * results["p_use"][dev][d,t]
                    sell_total[dev] += demands["weights"][d] * results["p_sell"][dev][d,t]
                    hp_total[dev] += demands["weights"][d] * results["p_hp"][dev][d,t]
    
    res["use_total"] = use_total
    res["sell_total"] = sell_total
    res["hp_total"] = hp_total

    if sum([power_total[dev] for dev in ("chp", "pv")]) == 0:
        self_consumption = 0
    else:
        self_consumption = (sum([use_total[dev] + hp_total[dev] for dev in ("chp", "bat", "pv")])
                            / (sum([power_total[dev] for dev in ("chp", "pv")])))
    res["self_consumption"] = self_consumption
    
    # Full load hours
    # Operating hours
    flh = {}
    oh = {}
    for dev in ("eh", "chp", "hp", "boiler"):
        flh[dev] = 0
        oh[dev] = 0
        for n in devs[dev].keys():
            if round(results["x"][dev][n]) == 1:
                for d in days:
                    for t in time_steps:
                        if results["heat"][dev][n][d,t] > 0:
                            oh[dev] += demands["weights"][d]
                            if dev == "eh":
                                flh[dev] += demands["weights"][d] * results["heat"][dev][n][d,t] / devs[dev][n]["Q_nom"]
                            elif dev in ("chp", "boiler"):
                                flh[dev] += demands["weights"][d] * results["heat"][dev][n][d,t] / max(devs[dev][n]["Q_heat"])
                            elif dev == "hp":
                                flh[dev] += demands["weights"][d] * results["heat"][dev][n][d,t] / max(devs[dev][n]["Q_heat"][d][t])
    res["full load hours"] = flh
    res["operating hours"] = oh
    res["z"] = results["z"]
    
    return res

def load_file(filename, devs, clustered):
    res = {}

    if filename == "kfw_with_restrictions.pkl":
        pass
    
    with open(filename, "rb") as fin:
        res["x"] = pkl.load(fin)
        res["y"] = pkl.load(fin)
        res["z"] = pkl.load(fin)
        res["x_tar"] = pkl.load(fin)
        res["x_gas"] = pkl.load(fin)
        res["x_el"] = pkl.load(fin)
        res["power"] = pkl.load(fin)
        res["heat"] = pkl.load(fin)
        res["energy"] = pkl.load(fin)
        res["p_grid"] = pkl.load(fin)
        res["G"] = pkl.load(fin)
        res["G_total"] = pkl.load(fin)
        res["El"] = pkl.load(fin)
        res["El_total"] = pkl.load(fin)
        res["soc"] = pkl.load(fin)
        res["soc_init"] = pkl.load(fin)
        res["charge"] = pkl.load(fin)
        res["discharge"] = pkl.load(fin)
        res["p_use"] = pkl.load(fin)
        res["p_sell"] = pkl.load(fin)
        res["p_hp"] = pkl.load(fin)
        res["c_inv"] = pkl.load(fin)
        res["c_om"] = pkl.load(fin)
        res["c_dem"] = pkl.load(fin)
        res["c_fix"] = pkl.load(fin)
        res["c_total"] = pkl.load(fin)
        res["rev"] = pkl.load(fin)
        res["sub"] = pkl.load(fin)
        res["emissions"] = pkl.load(fin)
        res["emissions_max"] = pkl.load(fin)
        res["objval"] = pkl.load(fin)
        res["runtime"] = pkl.load(fin)
        res["gap"] = pkl.load(fin)
    
    return evaluate(res, devs, clustered)

filenames = ["chp_with_kwkg.pkl", "chp_without_kwkg.pkl", "eeg_with_restrictions.pkl",
             "eeg_without_restrictions.pkl", "reference.pkl", "hp_with_tar.pkl", "hp_without_tar.pkl",
             "kfw_with_restrictions.pkl", "kfw_without_restrictions.pkl"]
keys = ["CHP (KWKG)", "CHP (without KWKG)", "PV (EEG)", "PV (without EEG)", 
        "Traditional", "HP (tariff)", "HP (without tariff)", "BAT (KfW)", "BAT (without KfW)"]
result = {}

#for filename in filenames:
for i in range(len(filenames)):
    filename = filenames[i]
    key = keys[i]
    result[key] = load_file(filename, devs, clustered)

fnames_moo = ["free_"+str(i)+".pkl" for i in range(10)]
res_moo = {}
for filename in fnames_moo:
    res_moo[filename] = load_file(filename, devs, clustered)


def write_first_col(sheet):
    sheet.write(0,0, "Scenario")
    sheet.write(1,0, "costs")
    sheet.write(2,0, "emissions")
    
    sheet.write(4,0, "c_inv")
    sheet.write(5,0, "c_om")
    sheet.write(6,0, "c_dem")
    sheet.write(7,0, "c_met")
    sheet.write(8,0, "rev")
    sheet.write(9,0, "sub")
    
    sheet.write(11,0, "cap boiler")
    sheet.write(12,0, "cap chp")
    sheet.write(13,0, "cap eh")
    sheet.write(14,0, "cap hp")
    sheet.write(15,0, "cap bat")
    sheet.write(16,0, "cap tes")
    sheet.write(17,0, "cap pv")
    sheet.write(18,0, "cap stc")
    sheet.write(19,0, "cap inv")
    
    sheet.write(21,0, "flh boiler")
    sheet.write(22,0, "flh chp")
    sheet.write(23,0, "flh eh")
    sheet.write(24,0, "flh hp")
    
    sheet.write(26,0, "oh boiler")
    sheet.write(27,0, "oh chp")
    sheet.write(28,0, "oh eh")
    sheet.write(29,0, "oh hp")
    
    sheet.write(31,0, "self_consumption")
    sheet.write(32,0, "heat_partial boiler")
    sheet.write(33,0, "heat_partial chp")
    sheet.write(34,0, "heat_partial eh")
    sheet.write(35,0, "heat_partial hp")
    sheet.write(36,0, "heat_partial stc")
    
    sheet.write(38,0, "heat total boiler")
    sheet.write(39,0, "heat total chp")
    sheet.write(40,0, "heat total eh")
    sheet.write(41,0, "heat total hp")
    sheet.write(42,0, "heat total stc")
    
    sheet.write(44,0, "power total pv")
    sheet.write(45,0, "power total chp")
    sheet.write(46,0, "power total hp")
    sheet.write(47,0, "power total grid")
    sheet.write(48,0, "power sell pv")
    sheet.write(49,0, "power use pv")
    sheet.write(50,0, "power sell chp")
    sheet.write(51,0, "power use chp")
    
    sheet.write(53,0, "run time")
    sheet.write(54,0, "gap")

def write_col(sheet, result, col, name):
    sheet.write(0,col, name)
    sheet.write(1,col, result["c_total"])
    sheet.write(2,col, result["emissions"])
    
    sheet.write(4,col, sum([result["c_inv"][key] for key in result["c_inv"].keys()]))
    sheet.write(5,col, sum([result["c_om"][key] for key in result["c_om"].keys()]))
    sheet.write(6,col, sum([result["c_dem"][key] for key in result["c_dem"].keys()]))
    sheet.write(7,col, sum([result["c_fix"][key] for key in result["c_fix"].keys()]))
    sheet.write(8,col, sum([result["rev"][key] for key in result["rev"].keys()]))
    sheet.write(9,col, sum([result["sub"][key] for key in result["sub"].keys()]))
    
    sheet.write(11,col, result["cap"]["boiler"])
    sheet.write(12,col, result["cap"]["chp"])
    sheet.write(13,col, result["cap"]["eh"])
    sheet.write(14,col, result["cap"]["hp"])
    sheet.write(15,col, result["cap"]["bat"])
    sheet.write(16,col, result["cap"]["tes"])
    sheet.write(17,col, result["cap"]["pv"])
    sheet.write(18,col, result["cap"]["stc"])
    sheet.write(19,col, result["cap"]["inv"])
    
    sheet.write(21,col, result["full load hours"]["boiler"])
    sheet.write(22,col, result["full load hours"]["chp"])
    sheet.write(23,col, result["full load hours"]["eh"])
    sheet.write(24,col, result["full load hours"]["hp"])
    
    sheet.write(26,col, result["operating hours"]["boiler"])
    sheet.write(27,col, result["operating hours"]["chp"])
    sheet.write(28,col, result["operating hours"]["eh"])
    sheet.write(29,col, result["operating hours"]["hp"])
    
    sheet.write(31,col, result["self_consumption"])
    sheet.write(32,col, result["heat_partials"]["boiler"])
    sheet.write(33,col, result["heat_partials"]["chp"])
    sheet.write(34,col, result["heat_partials"]["eh"])
    sheet.write(35,col, result["heat_partials"]["hp"])
    sheet.write(36,col, result["heat_partials"]["stc"])
    
    sheet.write(38,col, result["heat_total"]["boiler"])
    sheet.write(39,col, result["heat_total"]["chp"])
    sheet.write(40,col, result["heat_total"]["eh"])
    sheet.write(41,col, result["heat_total"]["hp"])
    sheet.write(42,col, result["heat_total"]["stc"])
    
    sheet.write(44,col, result["power_total"]["pv"])
    sheet.write(45,col, result["power_total"]["chp"])
    sheet.write(46,col, result["power_total"]["hp"])
    sheet.write(47,col, result["grid_total"]["house"] + result["grid_total"]["hp"])
    sheet.write(48,col, result["sell_total"]["pv"])
    sheet.write(49,col, result["use_total"]["pv"] + result["hp_total"]["pv"])
    sheet.write(50,col, result["sell_total"]["chp"])
    sheet.write(51,col, result["use_total"]["chp"] + result["hp_total"]["chp"])
    
    sheet.write(53,col, result["run time"])
    sheet.write(54,col, result["gap"])

if True:
    import xlsxwriter
    book = xlsxwriter.Workbook(filename="results.xlsx")
    
    sheet = book.add_worksheet("refs")
    write_first_col(sheet)
    for i in range(len(filenames)):
        filename = filenames[i]
        write_col(sheet, result[keys[i]], i+1, filename)
    
    moo_sheet = book.add_worksheet("moo")
    write_first_col(moo_sheet)
    for i in range(len(fnames_moo)):
        filename = fnames_moo[i]
        write_col(moo_sheet, res_moo[filename], i+1, filename)
    
    book.close()
    
    import matplotlib.pyplot as plt
    plt.figure()
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0,
                     box.width, box.height * 0.75])
    
    cos = []
    emi = []
    for i in range(len(fnames_moo)):
        cos.append(res_moo[fnames_moo[i]]["c_total"])
        emi.append(res_moo[fnames_moo[i]]["emissions"])
    plt.plot(cos, emi, label="Pareto set", color="black", linewidth=2)
    
    #markers = ["o", "D", "*", "s", "+", "h", "p", "1", "4", "3"]
    markers = ["o", "o", "D", "D", "+", "*", "*", "s", "s", "h", "p", "1", "4", "3"]
    #for i in range(len(result.keys())):
    #    name = (result.keys())[i]
    #    entry = result[name]
    #    plt.plot(entry["c_total"], entry["emissions"], linestyle="None", markersize=10, label=name, marker=markers[i])
    
    colors = ["#D7191C", "#ABD9E9", "#D7191C", "#ABD9E9", "black", "#D7191C", "#ABD9E9", "#D7191C", "#ABD9E9"]
    
    for i in range(len(keys)):
        name = keys[i]
        entry = result[name]
        color = colors[i]
        
        if name == "Traditional":
            plt.plot(entry["c_total"], entry["emissions"], linestyle="None", 
                     markersize=10, label=name, marker=markers[i], mfc=color, mec=color, mew=2)
        else:
            plt.plot(entry["c_total"], entry["emissions"], linestyle="None", 
                     markersize=10, label=name, marker=markers[i], mfc=color, mec=color)
    
    plt.legend(numpoints=1, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode="expand", borderaxespad=0.)
    
    plt.xlabel("Costs in Euro/a", fontsize=14)
    plt.ylabel("Emissions in tCO2/a", fontsize=14)
    
    plt.xlim([3000, 5000])
    plt.ylim([-4, 10])
    
    plt.savefig("pareto_set.pdf", bbox_inches='tight')
    plt.savefig("pareto_set.png", dpi=400, bbox_inches='tight')
    #plt.close()
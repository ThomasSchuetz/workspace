#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 05 11:04:37 2016

@author: tsz
"""

from __future__ import division

import matplotlib.pyplot as plt
import numpy as np

import xlrd

book = xlrd.open_workbook("formated_results.xlsx")
sheet_moo = book.sheet_by_name("moo")
sheet_other = book.sheet_by_name("refs")

keys_other = ["Traditional", "CHP with KWKG", "CHP w/o KWKG", "PV with EEG", "PV w/o EEG",
              "HP with tariff", "HP w/o tariff", "BAT with KfW", "BAT w/o KfW"]
indexes_moo = [1, 2, 3, 4, 7, 9, 10]
indexes_other = [5, 1, 2, 3, 4, 6, 7, 8, 9]
res_other = {}
res_moo = {}

#%% Read results for MOO
for i in indexes_moo:
    res_moo[i] = {"cap_heat":{}, "cap_el":{}, "cap_sto":{}, "self_cons":0, "heat_cov":{}}
    
    col = i
    
    res_moo[i]["cap_heat"]["boiler"] = sheet_moo.cell_value(11, col)
    res_moo[i]["cap_heat"]["chp"] = sheet_moo.cell_value(12, col)
    res_moo[i]["cap_heat"]["eh"] = sheet_moo.cell_value(13, col)
    res_moo[i]["cap_heat"]["hp"] = sheet_moo.cell_value(14, col)
    res_moo[i]["cap_heat"]["stc"] = sheet_moo.cell_value(18, col)
    
    res_moo[i]["cap_el"]["pv"] = sheet_moo.cell_value(17, col)
    res_moo[i]["cap_el"]["chp"] = sheet_moo.cell_value(17, col)
    
    res_moo[i]["cap_sto"]["bat"] = sheet_moo.cell_value(15, col)
    res_moo[i]["cap_sto"]["tes"] = sheet_moo.cell_value(16, col)
    
    res_moo[i]["self_cons"] = sheet_moo.cell_value(31, col)
    res_moo[i]["heat_cov"]["boiler"] = sheet_moo.cell_value(32, col)
    res_moo[i]["heat_cov"]["chp"] = sheet_moo.cell_value(33, col)
    res_moo[i]["heat_cov"]["eh"] = sheet_moo.cell_value(34, col)
    res_moo[i]["heat_cov"]["hp"] = sheet_moo.cell_value(35, col)
    res_moo[i]["heat_cov"]["stc"] = sheet_moo.cell_value(36, col)

#%% Read results for other scenarios
for i in range(len(keys_other)):
    res_other[keys_other[i]] = {"cap_heat":{}, "cap_el":{}, "cap_sto":{},
                                "self_cons":0, "heat_cov":{}}
    
    col = indexes_other[i]
    
    res_other[keys_other[i]]["cap_heat"]["boiler"] = sheet_other.cell_value(11, col)
    res_other[keys_other[i]]["cap_heat"]["chp"] = sheet_other.cell_value(12, col)
    res_other[keys_other[i]]["cap_heat"]["eh"] = sheet_other.cell_value(13, col)
    res_other[keys_other[i]]["cap_heat"]["hp"] = sheet_other.cell_value(14, col)
    res_other[keys_other[i]]["cap_heat"]["stc"] = sheet_other.cell_value(18, col)
    
    res_other[keys_other[i]]["cap_el"]["pv"] = sheet_other.cell_value(17, col)
    res_other[keys_other[i]]["cap_el"]["chp"] = sheet_other.cell_value(17, col)
    
    res_other[keys_other[i]]["cap_sto"]["bat"] = sheet_other.cell_value(15, col)
    res_other[keys_other[i]]["cap_sto"]["tes"] = sheet_other.cell_value(16, col)
    
    res_other[keys_other[i]]["self_cons"] = sheet_other.cell_value(31, col)
    res_other[keys_other[i]]["heat_cov"]["boiler"] = sheet_other.cell_value(32, col)
    res_other[keys_other[i]]["heat_cov"]["chp"] = sheet_other.cell_value(33, col)
    res_other[keys_other[i]]["heat_cov"]["eh"] = sheet_other.cell_value(34, col)
    res_other[keys_other[i]]["heat_cov"]["hp"] = sheet_other.cell_value(35, col)
    res_other[keys_other[i]]["heat_cov"]["stc"] = sheet_other.cell_value(36, col)

#%% Plot coverages (bars) and self-consumption rate (marker) for other MOO
fig, ax = plt.subplots(figsize=(8,8))

width = 0.5
index = np.arange(len(res_moo.keys()))+width/2

hc_boi_moo = np.array([res_moo[i]["heat_cov"]["boiler"] for i in indexes_moo])
hc_chp_moo = np.array([res_moo[i]["heat_cov"]["chp"] for i in indexes_moo])
hc_eh_moo = np.array([res_moo[i]["heat_cov"]["eh"] for i in indexes_moo])
hc_hp_moo = np.array([res_moo[i]["heat_cov"]["hp"] for i in indexes_moo])
hc_stc_moo = np.array([res_moo[i]["heat_cov"]["stc"] for i in indexes_moo])

self_cons_moo = np.array([res_moo[i]["self_cons"] for i in indexes_moo])

plt.bar(index, hc_boi_moo, width, color="#b30000", label="Cov. BOI")
plt.bar(index, hc_hp_moo, width, color="#e34a33", label="Cov. HP", bottom=hc_boi_moo)
plt.bar(index, hc_chp_moo, width, color="#fc8d59", label="Cov. CHP", bottom=hc_boi_moo+hc_hp_moo)
plt.bar(index, hc_stc_moo, width, color="#fdd49e", label="Cov. STC", bottom=hc_boi_moo+hc_hp_moo+hc_chp_moo)
plt.bar(index, hc_eh_moo, width, color="#fef0d9", label="Cov. EH", bottom=hc_boi_moo+hc_hp_moo+hc_chp_moo+hc_stc_moo)

# Self consumption rate
plt.plot(index+0.5*width, self_cons_moo, label="Self-consumption",
         ms=10, marker="+", mec="black", mfc="black", mew=2, linestyle="None")
#plt.plot(index[1:2]+0.5*width, self_cons_moo[1:2],
#         ms=10, marker="+", mec="black", mfc="black", mew=2, linestyle="None")

#plt.xlabel("Scenarios", fontsize=14)
labels_moo = [str(i) for i in range(len(indexes_moo))]
labels_moo[0] = "Min. costs"
labels_moo[-1] = "Min. emissions"
plt.xticks(index+0.5*width, labels_moo)

plt.ylabel("Rate in %", fontsize=14)
offset = 0.03
plt.ylim([-offset, 1+offset])
plt.yticks([0.2 * i for i in range(6)], [str(20 * i) for i in range(6)])

box = ax.get_position()
ax.set_position([box.x0, box.y0+0.1,
                 box.width, box.height * 0.75])
plt.legend(numpoints=1, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=3, mode="expand", borderaxespad=0.)

plt.savefig("plot_coverage_moo.pdf", bbox_inches='tight')
plt.savefig("plot_coverage_moo.png", dpi=400, bbox_inches='tight')

#%% Plot coverages (bars) and self-consumption rate (marker) for other scenearios
fig, ax = plt.subplots(figsize=(8,8))

width = 0.5
index = np.arange(len(res_other.keys()))+width/2

hc_boi_other = np.array([res_other[keys_other[i]]["heat_cov"]["boiler"] for i in range(len(keys_other))])
hc_chp_other = np.array([res_other[keys_other[i]]["heat_cov"]["chp"] for i in range(len(keys_other))])
hc_eh_other = np.array([res_other[keys_other[i]]["heat_cov"]["eh"] for i in range(len(keys_other))])
hc_hp_other = np.array([res_other[keys_other[i]]["heat_cov"]["hp"] for i in range(len(keys_other))])
hc_stc_other = np.array([res_other[keys_other[i]]["heat_cov"]["stc"] for i in range(len(keys_other))])

self_cons_other = np.array([res_other[keys_other[i]]["self_cons"] for i in range(len(keys_other))])

plt.bar(index, hc_boi_other, width, color="#b30000", label="Cov. BOI")
plt.bar(index, hc_hp_other, width, color="#e34a33", label="Cov. HP", bottom=hc_boi_other)
plt.bar(index, hc_chp_other, width, color="#fc8d59", label="Cov. CHP", bottom=hc_boi_other+hc_hp_other)
plt.bar(index, hc_stc_other, width, color="#fdd49e", label="Cov. STC", bottom=hc_boi_other+hc_hp_other+hc_chp_other)
plt.bar(index, hc_eh_other, width, color="#fef0d9", label="Cov. EH", bottom=hc_boi_other+hc_hp_other+hc_chp_other+hc_stc_other)

# Self consumption rate
plt.plot(index+0.5*width, self_cons_other, label="Self-consumption",
         ms=10, marker="+", mec="black", mfc="black", mew=2, linestyle="None")
#plt.plot(index[1:2]+0.5*width, self_cons_other[1:2],
#         ms=10, marker="+", mec="black", mfc="black", mew=2, linestyle="None")

#plt.xlabel("Scenarios", fontsize=14)
plt.xticks(index+0.5*width, keys_other, rotation="vertical")

plt.ylabel("Rate in %", fontsize=14)
offset = 0.03
plt.ylim([-offset, 1+offset])
plt.yticks([0.2 * i for i in range(6)], [str(20 * i) for i in range(6)])

box = ax.get_position()
ax.set_position([box.x0, box.y0+0.1,
                 box.width, box.height * 0.75])
plt.legend(numpoints=1, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=3, mode="expand", borderaxespad=0.)

plt.savefig("plot_coverage_other.pdf", bbox_inches='tight')
plt.savefig("plot_coverage_other.png", dpi=400, bbox_inches='tight')

##%% Plot TES and BAT capacities - MOO
fig, ax = plt.subplots(figsize=(8,8))
ax2 = ax.twinx()

width = 0.3
index = np.arange(len(res_moo.keys()))+width/2

cap_bat_moo = np.array([res_moo[i]["cap_sto"]["bat"] for i in indexes_moo])
cap_tes_moo = np.array([res_moo[i]["cap_sto"]["tes"] for i in indexes_moo])

#"#D7191C", "#ABD9E9"
bars1 = ax2.bar(index+width, cap_bat_moo, width, color="#ABD9E9", label="BAT")
bars2 = ax.bar(index, cap_tes_moo, width, color="#D7191C", label="TES")

ax2.set_ylabel("Battery capacity in kWh", fontsize=14)
ax.set_ylabel("TES water volume in m$^3$", fontsize=14)

plt.legend((bars2, bars1), ("TES", "BAT"), "upper left")

labels_moo = [str(i) for i in range(len(indexes_moo))]
labels_moo[0] = "Min. costs"
labels_moo[-1] = "Min. emissions"
plt.xticks(index+width, labels_moo)
ax.set_ylim([0, 2.05])
ax.set_xlim([0, index[-1]+2.5*width])

box = ax.get_position()
ax.set_position([box.x0, box.y0+0.1,
                 box.width, box.height * 0.75])
ax2.set_position([box.x0, box.y0+0.1,
                 box.width, box.height * 0.75])

plt.savefig("plot_cap_sto_moo.pdf", bbox_inches='tight')
plt.savefig("plot_cap_sto_moo.png", dpi=400, bbox_inches='tight')

#%% Plot TES and BAT capacities - other
fig, ax = plt.subplots(figsize=(8,8))
ax2 = ax.twinx()

width = 0.3
index = np.arange(len(res_other.keys()))+width/2

cap_bat_other = np.array([res_other[keys_other[i]]["cap_sto"]["bat"] for i in range(len(keys_other))])
cap_tes_other = np.array([res_other[keys_other[i]]["cap_sto"]["tes"] for i in range(len(keys_other))])

#"#D7191C", "#ABD9E9"
bars2 = ax.bar(index, cap_tes_other, width, color="#D7191C", label="TES")
bars1 = ax2.bar(index+width, cap_bat_other, width, color="#ABD9E9", label="BAT")

ax2.set_ylabel("Battery capacity in kWh", fontsize=14)
ax.set_ylabel("TES water volume in m$^3$", fontsize=14)

plt.legend((bars2, bars1), ("TES", "BAT"), "upper left")

ax.set_xticks(index+width)
ax.set_xticklabels(keys_other, rotation="vertical")
ax.set_ylim([0, 2.05])
ax2.set_ylim([0, 7.15])
ax.set_xlim([0, index[-1]+2.5*width])

box = ax.get_position()
ax.set_position([box.x0, box.y0+0.1,
                 box.width, box.height * 0.75])
ax2.set_position([box.x0, box.y0+0.1,
                 box.width, box.height * 0.75])

plt.savefig("plot_cap_sto_other.pdf", bbox_inches='tight')
plt.savefig("plot_cap_sto_other.png", dpi=400, bbox_inches='tight')

#%% Plot Heat generator's capacity - MOO
fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()

width = 0.12
index = np.arange(len(res_moo.keys()))+width/2

cap_boi_moo = np.array([res_moo[i]["cap_heat"]["boiler"] for i in indexes_moo])
cap_chp_moo = np.array([res_moo[i]["cap_heat"]["chp"] for i in indexes_moo])
cap_eh_moo = np.array([res_moo[i]["cap_heat"]["eh"] for i in indexes_moo])
cap_hp_moo = np.array([res_moo[i]["cap_heat"]["hp"] for i in indexes_moo])
cap_pv_moo = np.array([res_moo[i]["cap_el"]["pv"] for i in indexes_moo])
cap_stc_moo = np.array([res_moo[i]["cap_heat"]["stc"] for i in indexes_moo])

#"#D7191C", "#ABD9E9"
bars_boi = ax.bar(index+0*width, cap_boi_moo, width, color="#b30000", label="BOI")
bars_chp = ax.bar(index+1*width, cap_chp_moo, width, color="#e34a33", label="CHP")
bars_eh = ax.bar(index+2*width, cap_eh_moo, width, color="#fc8d59", label="EH")
bars_hp = ax.bar(index+3*width, cap_hp_moo, width, color="#fdd49e", label="HP")
bars_pv = ax2.bar(index+4*width, cap_pv_moo, width, color="#2171b5", label="PV", hatch="/")
bars_stc = ax2.bar(index+5*width, cap_stc_moo, width, color="#ABD9E9", label="STC", hatch="//")

ax.set_ylabel("Installed heating power in kW", fontsize=14)
ax2.set_ylabel("Covered roof area in m$^2$", fontsize=14)


labels_moo = [str(i) for i in range(len(indexes_moo))]
labels_moo[0] = "Min. costs"
labels_moo[-1] = "Min. emissions"
ax.set_xticks(index+3*width)
ax.set_xticklabels(labels_moo)
ax.set_xlim([0, index[-1]+6.5*width])
ax.set_ylim([0, 15])
ax.set_yticks([3*i for i in range(6)])
ax2.set_ylim([0, 75])

box = ax.get_position()
ax.set_position([box.x0, box.y0+0.1,
                 box.width, box.height * 0.75])
ax2.set_position([box.x0, box.y0+0.1,
                 box.width, box.height * 0.75])

plt.legend((bars_boi, bars_chp, bars_eh, bars_hp, bars_pv, bars_stc), 
           ("BOI", "CHP", "EH", "HP", "PV", "STC"),
           numpoints=1, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=3, mode="expand", borderaxespad=0.)

plt.savefig("plot_cap_gen_moo.pdf", bbox_inches='tight')
plt.savefig("plot_cap_gen_moo.png", dpi=400, bbox_inches='tight')

#%% Plot Heat generator's capacity - MOO
fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()

width = 0.12
index = np.arange(len(res_other.keys()))+width/2

cap_boi_moo = np.array([res_other[keys_other[i]]["cap_heat"]["boiler"] for i in range(len(keys_other))])
cap_chp_moo = np.array([res_other[keys_other[i]]["cap_heat"]["chp"] for i in range(len(keys_other))])
cap_eh_moo = np.array([res_other[keys_other[i]]["cap_heat"]["eh"] for i in range(len(keys_other))])
cap_hp_moo = np.array([res_other[keys_other[i]]["cap_heat"]["hp"] for i in range(len(keys_other))])
cap_stc_moo = np.array([res_other[keys_other[i]]["cap_heat"]["stc"] for i in range(len(keys_other))])
cap_pv_moo = np.array([res_other[keys_other[i]]["cap_el"]["pv"] for i in range(len(keys_other))])

#"#D7191C", "#ABD9E9"
bars_boi = ax.bar(index+0*width, cap_boi_moo, width, color="#b30000", label="BOI")
bars_chp = ax.bar(index+1*width, cap_chp_moo, width, color="#e34a33", label="CHP")
bars_eh = ax.bar(index+2*width, cap_eh_moo, width, color="#fc8d59", label="EH")
bars_hp = ax.bar(index+3*width, cap_hp_moo, width, color="#fdd49e", label="HP")
bars_pv = ax2.bar(index+4*width, cap_pv_moo, width, color="#2171b5", label="PV", hatch="/")
bars_stc = ax2.bar(index+5*width, cap_stc_moo, width, color="#ABD9E9", label="STC", hatch="//")

ax.set_ylabel("Installed heating power in kW", fontsize=14)
ax2.set_ylabel("Covered roof area in m$^2$", fontsize=14)

ax.set_xticks(index+3*width)
ax.set_xticklabels(keys_other, rotation="vertical")
ax.set_xlim([0, index[-1]+6.5*width])
ax.set_ylim([0, 15])
ax.set_yticks([3*i for i in range(6)])
ax2.set_ylim([0, 75])

box = ax.get_position()
ax.set_position([box.x0, box.y0+0.1,
                 box.width, box.height * 0.75])
ax2.set_position([box.x0, box.y0+0.1,
                 box.width, box.height * 0.75])

plt.legend((bars_boi, bars_chp, bars_eh, bars_hp, bars_pv, bars_stc),
           ("BOI", "CHP", "EH", "HP", "PV", "STC"),
           numpoints=1, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=3, mode="expand", borderaxespad=0.)

plt.savefig("plot_cap_gen_other.pdf", bbox_inches='tight')
plt.savefig("plot_cap_gen_other.png", dpi=400, bbox_inches='tight')


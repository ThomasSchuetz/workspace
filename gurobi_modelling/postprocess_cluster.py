# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 17:06:00 2016

@author: tsz
"""

from __future__ import division


import pickle
import numpy as np
import os

def read_single_entry(filename):

    # Load results
    with open(filename, "rb") as f_out:
        clustered = pickle.load(f_out)

    if "center" in filename or "medoid" in filename:
        z = clustered["z"]
        indexes = [i for i in range(365) if np.sum(z[i,:])>=1]
        clustered["z"] = [np.argmax(np.array(indexes)==np.argmax(z[:,i])) for i in range(365)]
    elif "365" in filename:
        clustered["z"] = range(365)
    elif "monthly" in filename:
        days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        clustered["z"] = np.hstack([np.tile(i, days[i]) for i in range(len(days))])
    
    # Compute full year load curve
    keys = ("electricity", "heat", "solar_irrad")
    clustered["ldc"] = [np.hstack([clustered[k][clustered["z"][i],:] for i in range(365)]) for k in keys]
    
    return clustered


def cleanup(dic, house_type, numbers, clus_types = ("mean", "center", "median", "medoid")):
    for number in numbers:
        for clus_type in clus_types:
            key = "res_"+house_type+"_"+clus_type+"_"+str(number)+".pkl"
            dic.pop(key)
    return dic

scale = lambda x, maxes: [x[i] / maxes[i] for i in range(3)]

def get_ref_prof(house="AB"):
    inputs = {}

    dhw  = np.loadtxt("raw_inputs/dhw_"+house+".csv")
    sh   = np.loadtxt("raw_inputs/heat_"+house+".csv")
    
    inputs["electricity"] = np.loadtxt("raw_inputs/elec_"+house+".csv")
    inputs["heat"]        = dhw + sh
    inputs["solar_irrad"] = np.loadtxt("raw_inputs/solar_rad_35deg.csv")
    inputs["temperature"] = np.loadtxt("raw_inputs/temperature.csv")
    
    
    ins = [inputs["electricity"], inputs["heat"], inputs["solar_irrad"]]
    maxes = [np.max(i) for i in ins]
    scaled_ins = scale(ins, maxes)
    
    return scaled_ins, maxes

(scaled_ins_AB, max_AB) = get_ref_prof("AB")
(scaled_ins_MFH, max_MFH) = get_ref_prof("MFH")
(scaled_ins_SFH, max_SFH) = get_ref_prof("SFH")

results_all = {}
results = {"AB":{}, "MFH":{}, "SFH":{}}

filenames = [file for file in os.listdir("res_clus") if file.endswith(".pkl")]

for filename in filenames:
    
    fname = "res_clus/" + filename

    results_all[filename] = read_single_entry(fname)

for k_res in results.keys():
    for k_all in results_all.keys():
        if k_res in k_all:
            results[k_res][k_all] = results_all[k_all]

import post_process_helpers as pph

days = 25

ldc_AB = pph.get_results_ordered(results["AB"], days=days, key="ldc", prefix="clus")
fdid_AB = pph.clus_fdid(ldc_AB, ref_max=max_AB, ref_prof=scaled_ins_AB, days=days)
pph.plot_ordered(fdid_AB, ytitle="SSE AB", relative=False, filename="ab_sse_900.png", pos_legend=0)

ldc_SFH = pph.get_results_ordered(results["SFH"], days=days, key="ldc", prefix="clus")
fdid_SFH = pph.clus_fdid(ldc_SFH, ref_max=max_SFH, ref_prof=scaled_ins_SFH, days=days)
pph.plot_ordered(fdid_SFH, ytitle="SSE SFH", relative=False, filename="sfh_sse_900.png", pos_legend=0)

#
#center = results_all["clus_AB_center_7.pkl"]
#
#scaled_center = scale(center["ldc"], maxes)
#
#deviation = np.sum([np.square(np.array(scaled_ins) - np.array(scaled_center))])
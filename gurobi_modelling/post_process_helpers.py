# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 14:17:02 2016

@author: tsz
"""

from __future__ import division

import matplotlib.pyplot as plt
import numpy as np

def clus_fdid(results, ref_max, ref_prof, days=20):
    res = {}
    
    f = [results["monthly"][j] / ref_max[j] for j in range(3)]
    res["monthly"] = _sse(ref_prof, f)
    f = [results["reference"][j] / ref_max[j] for j in range(3)]
    res["reference"] = _sse(ref_prof, f)
    
    types = ["center", "mean", "median", "medoid"]
#    types = ["center", "medoid"]
    for t in types:
        temp = []
        for i in range(days):
            f = [results[t][i,j,:] / ref_max[j] for j in range(3)]
            temp.append(_sse(ref_prof, f))
        res[t] = np.array(temp)
        
    return res
    
def _sse(data1, data2):
    
    
    return np.sum([np.square(np.array(data1) - np.array(data2))])

def fdid(results, days=25):
    """
    """
    scale = lambda x: [k / np.max(k) for k in x["ldc"]]
    
    if "AB" in results.keys()[0]:
        house = "AB"
    elif "SFH" in results.keys()[0]:
        house = "SFH"
    else:
        house = "MFH"
    
    types = ["center", "mean", "median", "medoid"]
    
    ref_scaled = scale(results["orig_res_"+house+"_reference_365.pkl"])
    
    res = {}
    res["reference"] = _sse(ref_scaled, scale(results["orig_res_"+house+"_reference_365.pkl"]))
    res["unconnected"] = _sse(ref_scaled, scale(results["res_"+house+"_reference_365_not_connected.pkl"]))
    res["monthly"] = _sse(ref_scaled, scale(results["res_"+house+"_monthly_12.pkl"]))
    
    for t in types:
        temp = []
        for i in range(days):
            temp.append(_sse(ref_scaled, scale(results["res_"+house+"_"+t+"_"+str(i+1)+".pkl"])))
        res[t] = np.array(temp)
        
    return res

def get_results_ordered(results, days=25, key="objective", prefix="res"):
    """
    """

    if "AB" in results.keys()[0]:
        house = "AB"
    elif "SFH" in results.keys()[0]:
        house = "SFH"
    else:
        house = "MFH"
    
    types = ["center", "mean", "median", "medoid"]
#    types = ["center", "medoid"]
    
    res = {}
    if prefix == "res_clus":
        res["reference"] = results["clus_"+house+"_reference_365.pkl"][key]
    else:
        res["reference"] = results["orig_res_"+house+"_reference_365.pkl"][key]
        res["unconnected"] = results["res_"+house+"_reference_365_not_connected.pkl"][key]
        
    
    res["monthly"] = results[prefix+"_"+house+"_monthly_12.pkl"][key]
    
    for t in types:
        temp = []
        for i in range(days):
            temp.append(results[prefix+"_"+house+"_"+t+"_"+str(i+1)+".pkl"][key])
        res[t] = np.array(temp)
        
    return res

def plot_ordered(res, ytitle="", relative=True, filename="result.pdf", pos_legend=0):
    if relative:
        reference = res["reference"]
        for key in res.keys():
            res[key] = (reference-res[key]) / reference

    days = len(res["center"])
    colors = ["#d7191c", "#fdae61", "#abd9e9", "#2c7bb6"]
    i = 0
    
#    keys = ["center", "mean", "median", "medoid", "monthly", "reference"]
    keys = ["center", "mean", "median", "medoid", "monthly"]
#    keys = ["center", "medoid", "monthly"]
    
    plt.figure()
    for key in keys:
        if key in ("reference", "unconnected", "monthly"):
            if key == "reference":
                plt.plot([1, days], [res["reference"], res["reference"]],
                         linestyle = ":", color = "grey")#, label="reference")
            elif key == "monthly":
                plt.plot(12, res["monthly"], linestyle="None",
                         marker = "o", color = "black", label="monthly")
            else:
                pass
        else:
            plt.plot(np.arange(1,days+1), res[key], 
                     linestyle="-", linewidth=2, 
                     color=colors[i], label="k-"+key+"s")
            i += 1
    plt.legend(numpoints=1, loc=pos_legend)
    plt.xlim([1, days])
    plt.ylim([-0.05, 0.2])
    plt.xticks([1] + [(j+1) * 5 for j in range(5)])
    
    plt.xlabel("Demand days", fontsize=12)
    plt.ylabel(ytitle, fontsize=12)

    plt.savefig(filename, dpi=400)

def plot_times(res, ytitle="", filename="result.pdf", pos_legend=0, threshold=300):
    times = []
    days = len(res["center"])
    keys = ["center", "mean", "median", "medoid"]
    for i in range(days):
        if i == 11:
            temp = res["monthly"]
            counter = 1
        else:
            temp = 1
            counter = 0
        
        for k in keys:
            if res[k][i] < threshold:
                temp *= res[k][i]
                counter += 1
        
        times.append(np.power(temp, 1/counter))
    
    plt.figure()
    plt.plot(np.arange(1,days+1), times, linestyle="-", linewidth=2, color="black")
    plt.xlim([1, days])
    plt.xticks([1] + [(j+1) * 5 for j in range(5)])
    plt.ylim([0,100])
    
    plt.xlabel("Demand days", fontsize=12)
    plt.ylabel(ytitle, fontsize=12)

    plt.savefig(filename, dpi=400)
        
        
        
        
        
        
        
        
        
        
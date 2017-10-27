#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import pickle
import clustering

# Bounds for the number of typical days
min_days = 1
max_days = 12

def cluster(number_clusters, google_or=False):
    raw_inputs = {}
    raw_inputs["electricity"] = np.loadtxt("raw_inputs/electricity.csv")
    raw_inputs["dhw"]         = np.loadtxt("raw_inputs/dhw.csv")
    raw_inputs["sh"]          = np.loadtxt("raw_inputs/space_heating.csv")
    raw_inputs["heat"]        = raw_inputs["dhw"] + raw_inputs["sh"]
    raw_inputs["solar_irrad"] = np.loadtxt("raw_inputs/solar_rad_35deg.csv") / 1000
    raw_inputs["solar_irrad"] = np.maximum(raw_inputs["solar_irrad"], 0)
    raw_inputs["temperature"] = np.loadtxt("raw_inputs/temperature.csv")
    
###############################################################################
    # Clustering
    inputs_clustering = np.array([raw_inputs["electricity"], 
                                  raw_inputs["heat"], 
                                  raw_inputs["solar_irrad"]])
    inputs_additional = np.array([raw_inputs["temperature"]])
    
    clus_res = clustering.cluster(inputs_clustering, 
                                  inputs_additional=inputs_additional,
                                  method="medoid",
                                  number_clusters=number_clusters,
                                  time_limit=600,
                                  mip_gap=0.0,
                                  google_or=google_or)
    (inputs, typ_inputs_add, nc, scaling_factors, z, times, obj, gap) = clus_res
    
    filename = ("results/res_googleOR_" + str(google_or) + "_" + str(number_clusters) + "days.pkl")
    with open(filename, "wb") as f_in:
        pickle.dump(times, f_in, pickle.HIGHEST_PROTOCOL)
        pickle.dump(obj, f_in, pickle.HIGHEST_PROTOCOL)
        pickle.dump(gap, f_in, pickle.HIGHEST_PROTOCOL)

for number_clusters in range(min_days, max_days+1):
    cluster(number_clusters, google_or=True)
    cluster(number_clusters, google_or=False)
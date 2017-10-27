#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 10:09:31 2017

@author: Thomas
"""

from __future__ import division
import pickle

# Bounds for the number of typical days
min_days = 1
max_days = 12

def read_res(filename):
    with open(filename, "rb") as f_in:
        times = pickle.load(f_in)
        obj = pickle.load(f_in)
        gap = pickle.load(f_in)
    
    return times, obj, gap


res = {"google_or": {}, "gurobi": {}}

for number_clusters in range(min_days, max_days+1):
    res["google_or"][number_clusters] = read_res("res_googleOR_True_"+str(number_clusters)+"days.pkl")
    res["gurobi"][number_clusters] = read_res("res_googleOR_False_"+str(number_clusters)+"days.pkl")
    
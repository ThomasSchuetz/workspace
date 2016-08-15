#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 12:16:05 2016

@author: tsz
"""
import pickle as pkl
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
    print res["z"]
    


filename = "chp_with_kwkg.pkl"
print filename
load_file(filename, 0, 0)

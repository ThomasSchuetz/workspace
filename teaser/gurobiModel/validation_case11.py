#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 12:17:52 2016

@author: tsz
"""
from __future__ import division

import numpy as np
import gurobipy as gp
import matplotlib.pyplot as plt

import twoElements
import testcases as tc

import time

model = gp.Model("tc11")

# Definition of time horizon
times_per_hour = 1
timesteps = 24 * 60 * times_per_hour # 60 days
timesteps_day = int(24 * times_per_hour)

# Zero inputs    
ventRate = np.zeros(timesteps)
solRad   = np.zeros((timesteps,1))
intGainsConv = np.zeros(timesteps)
window   = np.zeros(timesteps)
windowIndoorSurface     = np.zeros(timesteps)
extWallIndoorSurface    = np.zeros(timesteps)

# Constant inputs
extWall = np.zeros(timesteps) + 295.15 # all temperatures in K
Tv      = np.zeros(timesteps) + 295.15 # in K

# Initial Values
T_start = 273.15 + 22

# Variable inputs
intGainsRad = np.zeros(timesteps_day)
for q in range(int(6*timesteps_day/24), int(18*timesteps_day/24)):
    intGainsRad[q] = 1000
intGainsRad = np.tile(intGainsRad, 60)

# Load constant house parameters
params = tc.get_house_data(timesteps, case=11)

# Ventilation and air parameters
ports = {}
ports["Tv"] = Tv
ports["ventRate"] = ventRate
ports["rhoair"] = 1.19 # kg/m3
ports["cair"] = 1007 #
ports["heaterCooler"] = False
ports["setAirTemp"] = False

# Load initial values into the model
Tair = model.addVar(vtype="C", name="Tair_start", lb=-100.)
Tow  = model.addVar(vtype="C", name="Tow_start", lb=-100.)
Tiw  = model.addVar(vtype="C", name="Tiw_start", lb=-100.)
model.update()
model.addConstr(Tair == T_start)
model.addConstr(Tow == T_start)
model.addConstr(Tiw == T_start)
model.update()

# Load initial values into the model
Tair = {}
Q_HC = {}
Tpre = {}
dT   = {}
intWallIndoorSurface = {}
for t in range(timesteps):
    Tair[t] = model.addVar(vtype="C", name="Tair_"+str(t), lb=-100.)
    Q_HC[t] = model.addVar(vtype="C", name="Q_HC_"+str(t), lb=0., ub=500.)
    Tpre[t] = model.addVar(vtype="C", name="Tpre_"+str(t), lb=-100.)
    dT[t]   = model.addVar(vtype="C", name="dT_"+str(t))
    intWallIndoorSurface[t] = model.addVar(vtype="C", name="intWallIndoorSurface_"+str(t),
                                           lb=-500., ub=0.)

# preset temperature
t_set = np.zeros(timesteps_day) + 273.15 + 22
for q in range(int(6*timesteps_day/24), int(18*timesteps_day/24)):
    t_set[q] = 273.15 + 27
t_set = np.tile(t_set, 60)

Tair_start = model.addVar(vtype="C", name="Tair_start", lb=-100.)
Tow        = model.addVar(vtype="C", name="Tow_start", lb=-100.)
Tiw        = model.addVar(vtype="C", name="Tiw_start", lb=-100.)
model.update()

for t in range(timesteps):
    model.addConstr(Tpre[t] == t_set[t])
    model.addConstr(dT[t] >= Tair[t] - Tpre[t]) # Tair > Tpre
    model.addConstr(dT[t] >= Tpre[t] - Tair[t]) # Tair < Tpre
model.addConstr(Tair_start == T_start)
model.addConstr(Tow        == T_start)
model.addConstr(Tiw        == T_start)

# objective: minimize sum of temperature deviations between preset and actual temperature
model.setObjective(sum(dT[t] for t in range(timesteps)),gp.GRB.MINIMIZE)
model.update()

now = time.time()

# Calculate indoor air temperature
model = twoElements.twoElements(params, solRad, window, extWall, 
                                                  windowIndoorSurface, 
                                                  extWallIndoorSurface,
                                                  intWallIndoorSurface, 
                                                  intGainsConv, intGainsRad, ports,
                                                  model,
                                                  dt=int(3600/times_per_hour))

dt=int(3600/times_per_hour)                                                 
T_air = {}
Q_HC = {}
intWallIndoorSurface = {}
Q_iw = {}
Q_ow = {}
Tiw  = {}
Tow  = {}
Tow_start = model.getVarByName("Tow_start")
Tiw_start = model.getVarByName("Tiw_start")
for t in range(timesteps):
    T_air[t] = model.getVarByName("Tair_"+str(t)).X
    Q_HC[t] = (model.getVarByName("Q_HC_"+str(t))).X
    intWallIndoorSurface[t] = (model.getVarByName("intWallIndoorSurface_"+str(t))).X
    Tiw[t] = model.getVarByName("Tiw_"+str(t))
    Tow[t] = model.getVarByName("Tow_"+str(t))
    
    if t > 0:
        Q_ow[t] = (Tow[t].X-Tow[t-1].X)*params["CExt"]/dt
        Q_iw[t] = (Tiw[t].X-Tiw[t-1].X)*params["CInt"]/dt
    else:
        Q_ow[t] = (Tow[t].X-Tow_start.X)*params["CExt"]/dt
        Q_iw[t] = (Tiw[t].X-Tiw_start.X)*params["CInt"]/dt
Q_hc = np.array([Q_HC[t] for t in range(timesteps)]) + np.array([intWallIndoorSurface[t] for t in range(timesteps)])
Q_iw = np.array([Q_iw[t] for t in range(timesteps)])
Q_ow = np.array([Q_ow[t] for t in range(timesteps)])
T_air = np.array([T_air[t] for t in range(timesteps)])

print 
print "Time used: " +str(time.time()-now)
print

# Compute averaged results
Q_hc_mean = np.array([np.mean(Q_hc[i*times_per_hour:(i+1)*times_per_hour]) for i in range(24*60)])
Q_iw_mean = np.array([np.mean(Q_iw[i*times_per_hour:(i+1)*times_per_hour]) for i in range(24*60)])
Q_ow_mean = np.array([np.mean(Q_ow[i*times_per_hour:(i+1)*times_per_hour]) for i in range(24*60)])

Q_hc_1 = Q_hc_mean[0:24] + Q_iw_mean[0:24] + Q_ow_mean[0:24]
Q_hc_10 = Q_hc_mean[216:240] + Q_iw_mean[216:240] + Q_ow_mean[216:240]
Q_hc_60 = Q_hc_mean[1416:1440] + Q_iw_mean[1416:1440] + Q_ow_mean[1416:1440]

T_air_c = T_air - 273.15
T_air_mean = np.array([np.mean(T_air_c[i*times_per_hour:(i+1)*times_per_hour]) for i in range(24*60)])

T_air_1 = T_air_mean[0:24]
T_air_10 = T_air_mean[216:240]
T_air_60 = T_air_mean[1416:1440]

# Load reference results    
(load_res_1, load_res_10, load_res_60) = tc.load_res("inputs/case11_res.csv")
Q_hc_ref_1 = load_res_1[:,1]
Q_hc_ref_10 = load_res_10[:,1]
Q_hc_ref_60 = load_res_60[:,1]

T_air_ref_1 = load_res_1[:,0]
T_air_ref_10 = load_res_10[:,0]
T_air_ref_60 = load_res_60[:,0]


# Plot comparisons
def plot_result(res, ref, title="Results day 1"):
    plt.figure()
    ax_top = plt.subplot(211)
    plt.plot(ref, label="Reference", color="black", linestyle="--")
    plt.plot(res, label="Simulation", color="blue", linestyle="-")
    plt.legend()
    plt.ylabel("Heat load in W")
    
    plt.title(title)

    plt.subplot(212, sharex=ax_top)
    plt.plot(res-ref, label="Ref. - Sim.")
    plt.legend()
    plt.ylabel("Heat load difference in W")
    plt.xticks([4*i for i in range(7)])
    plt.xlim([1,24])
    plt.xlabel("Time in h")

plot_result(T_air_1, T_air_ref_1, "Results temperatures day 1")
plot_result(T_air_10, T_air_ref_10, "Results temperatures day 10")
plot_result(T_air_60, T_air_ref_60, "Results temperatures day 60")

plot_result(Q_hc_1, Q_hc_ref_1, "Results heating/cooling day 1")
plot_result(Q_hc_10, Q_hc_ref_10, "Results heating/cooling day 10")
plot_result(Q_hc_60, Q_hc_ref_60, "Results heating/cooling day 60")

print("Deviations temperature in K:")
print("Max. deviation day 1: " + str(np.max(np.abs(T_air_1 - T_air_ref_1))))
print("Max. deviation day 10: " + str(np.max(np.abs(T_air_10 - T_air_ref_10))))
print("Max. deviation day 60: " + str(np.max(np.abs(T_air_60 - T_air_ref_60))))
print("")
print("Deviations heating/cooling in W:")
print("Max. deviation day 1: " + str(np.max(np.abs(Q_hc_1 - Q_hc_ref_1))))
print("Max. deviation day 10: " + str(np.max(np.abs(Q_hc_10 - Q_hc_ref_10))))
print("Max. deviation day 60: " + str(np.max(np.abs(Q_hc_60 - Q_hc_ref_60))))
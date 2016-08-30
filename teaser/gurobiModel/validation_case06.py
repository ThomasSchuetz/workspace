#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 13:26:01 2016

@author: tsz
"""
from __future__ import division

import numpy as np
import gurobipy as gp
import matplotlib.pyplot as plt

import twoElements
import testcases as tc

model = gp.Model("tc6")

# Definition of time horizon
times_per_hour = 10
timesteps = 24 * 60 * times_per_hour # 60 days
timesteps_day = int(24 * times_per_hour)

# Zero inputs    
ventRate     = np.zeros(timesteps)
solRad       = np.zeros((timesteps,1))
intGainsConv = np.zeros(timesteps)
window       = np.zeros(timesteps)
windowIndoorSurface     = np.zeros(timesteps)
extWallIndoorSurface    = np.zeros(timesteps)
intWallIndoorSurface    = np.zeros(timesteps)

# Constant inputs
extWall     = np.zeros(timesteps) + 295.15 # all temperatures in K
Tv          = np.zeros(timesteps) + 295.15 # in K

# Initial Values
T_start = 273.15 + 22

# Variable inputs
intGainsRad = np.zeros(timesteps_day)
for q in range(int(6*timesteps_day/24), int(18*timesteps_day/24)):
    intGainsRad[q] = 1000
intGainsRad = np.tile(intGainsRad, 60)

# Load constant house parameters
params = tc.get_house_data(timesteps, case=6)

# Ventilation and air parameters
ports = {}
ports["Tv"] = Tv
ports["ventRate"] = ventRate
ports["rhoair"] = 1.19 # kg/m3
ports["cair"] = 0 #
ports["heaterCooler"] = True
ports["setAirTemp"] = True

# Load initial values into the model
Tair = {}
Q_HC = {}
#y    = {}
z    = {}
for t in range(timesteps):
    Tair[t] = model.addVar(vtype="C", name="Tair_"+str(t), lb=-100.)
    Q_HC[t] = model.addVar(vtype="C", name="Q_HC_"+str(t), lb=-1e5)
#    y[t]    = model.addVar(vtype="B", name="y_"+str(t)) # =1 if abs(Q_HC[t] is positive)
    z[t]    = model.addVar(vtype="C", name="z_"+str(t))

t_set = np.zeros(timesteps_day) + 273.15 + 22
for q in range(int(6*timesteps_day/24), int(18*timesteps_day/24)):
    t_set[q] = 273.15 + 27
t_set = np.tile(t_set, 60)

Tow  = model.addVar(vtype="C", name="Tow_0", lb=-100.)
Tiw  = model.addVar(vtype="C", name="Tiw_0", lb=-100.)
model.update()

for t in range(timesteps):
    model.addConstr(Tair[t] == t_set[t])
    model.addConstr(z[t] >= Q_HC[t])
    model.addConstr(z[t] >= -Q_HC[t])
model.addConstr(Tow == T_start)
model.addConstr(Tiw == T_start)

model.setObjective(sum(z[t] for t in range(timesteps)),gp.GRB.MINIMIZE)
model.update()

# Calculate indoor air temperature
T_air, Q_HC, Q_iw, Q_ow = twoElements.twoElements(params, solRad, window, extWall, 
                                                  windowIndoorSurface, 
                                                  extWallIndoorSurface,
                                                  intWallIndoorSurface, 
                                                  intGainsConv, intGainsRad, ports,
                                                  model,
                                                  dt=int(3600/times_per_hour))

# Compute averaged results
Q_hc_mean = np.array([np.mean(Q_HC[i*times_per_hour:(i+1)*times_per_hour]) for i in range(24*60)])

Q_hc_1 = Q_hc_mean[0:24]
Q_hc_10 = Q_hc_mean[216:240]
Q_hc_60 = Q_hc_mean[1416:1440]

# Load reference results    
(Q_hc_ref_1, Q_hc_ref_10, Q_hc_ref_60) = tc.load_res("inputs/case06_res.csv")
Q_hc_ref_1 = -Q_hc_ref_1[:,0]
Q_hc_ref_10 = -Q_hc_ref_10[:,0]
Q_hc_ref_60 = -Q_hc_ref_60[:,0]


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

plot_result(Q_hc_1, Q_hc_ref_1, "Results day 1")
plot_result(Q_hc_10, Q_hc_ref_10, "Results day 10")
plot_result(Q_hc_60, Q_hc_ref_60, "Results day 60")

print("Max. deviation day 1: " + str(np.max(np.abs(Q_hc_1 - Q_hc_ref_1))))
print("Max. deviation day 10: " + str(np.max(np.abs(Q_hc_10 - Q_hc_ref_10))))
print("Max. deviation day 60: " + str(np.max(np.abs(Q_hc_60 - Q_hc_ref_60))))
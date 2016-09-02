#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 11:46:59 2016

@author: tsz
"""
from __future__ import division

import numpy as np
import gurobipy as gp
import matplotlib.pyplot as plt

import twoElements
import testcases as tc
import eqAirTemp

model = gp.Model("tc9")

# Definition of time horizon
times_per_hour = 15
timesteps = 24 * 60 * times_per_hour # 60 days
timesteps_day = int(24 * times_per_hour)

# Zero inputs    
ventRate = np.zeros(timesteps)
window   = np.zeros(timesteps)
windowIndoorSurface     = np.zeros(timesteps)
extWallIndoorSurface    = np.zeros(timesteps)
intWallIndoorSurface    = np.zeros(timesteps)

# Initial Values
T_start = 273.15 + 22

# Variable inputs
intGainsConv = np.zeros(timesteps_day)
intGainsRad  = np.zeros(timesteps_day)
for q in range(int(7*timesteps_day/24), int(17*timesteps_day/24)):
    intGainsConv[q] = 80 + 200
    intGainsRad [q] = 80
intGainsConv = np.tile(intGainsConv, 60)
intGainsRad  = np.tile(intGainsRad , 60)

q_sol_rad_win_raw = np.loadtxt("inputs/case08_q_sol_win.csv", usecols=(1,2))
solarRad_win = q_sol_rad_win_raw[0:24,:]
solarRad_win[solarRad_win > 100] = solarRad_win[solarRad_win > 100] * 0.15
solarRad_win_adj = np.repeat(solarRad_win, times_per_hour, axis=0)
solRad = np.tile(solarRad_win_adj.T, 60).T

sunblind_in = np.zeros_like(solarRad_win)
sunblind_in[solarRad_win > 100] = 0.85
sunblind_in_adj = np.repeat(sunblind_in, times_per_hour, axis=0)
sunblind_in_tiled = np.tile(sunblind_in_adj.T, 60).T

q_sol_rad_wall_raw = np.loadtxt("inputs/case08_q_sol_wall.csv", usecols=(1,2))
solarRad_wall = q_sol_rad_wall_raw[0:24,:]
solarRad_wall_adj = np.repeat(solarRad_wall, times_per_hour, axis=0)
solarRad_wall_tiled = np.tile(solarRad_wall_adj.T, 60).T

t_outside_raw = np.loadtxt("inputs/case08_t_amb.csv", delimiter=",")
t_outside = ([t_outside_raw[2*i,1] for i in range(24)])
t_outside_adj = np.repeat(t_outside, times_per_hour)
Tv = np.tile(t_outside_adj, 60)

H_sky_raw = np.loadtxt("inputs/case09_h_sky.csv", usecols=(1,))
H_sky = H_sky_raw[0:24]
t_black_sky_in = 65.99081593 * (H_sky ** 0.25)
t_black_sky_adj = np.repeat(t_black_sky_in, times_per_hour)
t_black_sky = np.tile(t_black_sky_adj, 60)

extWall = eqAirTemp.equal_air_temp(solarRad_wall_tiled, 
                                   t_black_sky, 
                                   Tv, 
                                   sunblind_in_tiled, 
                                   tc.get_eqAirTemp_params(case=9))

# Load constant house parameters
params = tc.get_house_data(timesteps, case=9)

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

# Calculate indoor air temperature
T_air, Q_HC, Q_iw, Q_ow = twoElements.twoElements(params, solRad, window, extWall, 
                                                  windowIndoorSurface, 
                                                  extWallIndoorSurface,
                                                  intWallIndoorSurface, 
                                                  intGainsConv, intGainsRad, ports,
                                                  model,
                                                  dt=int(3600/times_per_hour))

# Compute averaged results
T_air_c = T_air - 273.15
T_air_mean = np.array([np.mean(T_air_c[i*times_per_hour:(i+1)*times_per_hour]) for i in range(24*60)])

T_air_1 = T_air_mean[0:24]
T_air_10 = T_air_mean[216:240]
T_air_60 = T_air_mean[1416:1440]

# Load reference results    
(T_air_ref_1, T_air_ref_10, T_air_ref_60) = tc.load_res("inputs/case09_res.csv")
T_air_ref_1 = T_air_ref_1[:,0]
T_air_ref_10 = T_air_ref_10[:,0]
T_air_ref_60 = T_air_ref_60[:,0]


# Plot comparisons
def plot_result(res, ref, title="Results day 1"):
    plt.figure()
    ax_top = plt.subplot(211)
    plt.plot(res, label="Reference", color="black", linestyle="--")
    plt.plot(ref, label="Simulation", color="blue", linestyle="-")
    plt.legend()
    plt.ylabel("Temperature in degC")
    
    plt.title(title)

    plt.subplot(212, sharex=ax_top)
    plt.plot(res-ref, label="Ref. - Sim.")
    plt.legend()
    plt.ylabel("Temperature difference in K")
    plt.xticks([4*i for i in range(7)])
    plt.xlim([1,24])
    plt.xlabel("Time in h")

plot_result(T_air_1, T_air_ref_1, "Results day 1")
plot_result(T_air_10, T_air_ref_10, "Results day 10")
plot_result(T_air_60, T_air_ref_60, "Results day 60")

print("Max. deviation day 1: " + str(np.max(np.abs(T_air_1 - T_air_ref_1))))
print("Max. deviation day 10: " + str(np.max(np.abs(T_air_10 - T_air_ref_10))))
print("Max. deviation day 60: " + str(np.max(np.abs(T_air_60 - T_air_ref_60))))
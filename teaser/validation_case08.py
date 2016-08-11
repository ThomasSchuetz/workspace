#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 17:41:23 2016

@author: tsz
"""
from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

import low_order_VDI
import validationVDITestcases as tc
import equal_air_temperature as eq_air_temp

# Definition of time horizon
times_per_hour = 60
timesteps = 24 * 60 * times_per_hour # 60 days
timesteps_day = int(24 * times_per_hour)

# Zero inputs    
ventRate = np.zeros(timesteps)

# Constant inputs
alphaRad = np.zeros(timesteps) + 5
t_black_sky = np.zeros(timesteps) + 273.15

# Variable inputs
Q_ig = np.zeros(timesteps_day)
source_igRad = np.zeros(timesteps_day)
for q in range(int(7*timesteps_day/24), int(17*timesteps_day/24)):
    Q_ig[q] = 200 + 80
    source_igRad[q] = 80
Q_ig = np.tile(Q_ig, 60)
source_igRad = np.tile(source_igRad, 60)

q_sol_rad_win_raw = np.loadtxt("inputs/case08_q_sol_win.csv", usecols=(1,2))
solarRad_win = q_sol_rad_win_raw[0:24,:]
solarRad_win[solarRad_win > 100] = solarRad_win[solarRad_win > 100] * 0.15
solarRad_win_adj = np.repeat(solarRad_win, times_per_hour, axis=0)
solarRad_win_in = np.tile(solarRad_win_adj.T, 60).T

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
weatherTemperature = np.tile(t_outside_adj, 60)

equalAirTemp = eq_air_temp.equal_air_temp(solarRad_wall_tiled, 
                                          t_black_sky, 
                                          weatherTemperature, 
                                          sunblind_in_tiled, 
                                          tc.get_eqAirTemp_params(case=8))

# Load constant house parameters
houseData = tc.get_house_data(case=8)

krad = 1

# Define set points (prevent heating or cooling!)
t_set_heating = np.zeros(timesteps)        # in Kelvin
t_set_cooling = np.zeros(timesteps) + 600  # in Kelvin

heater_limit = np.zeros(timesteps) + 1e10
cooler_limit = np.zeros(timesteps) - 1e10

# Calculate indoor air temperature
T_air, Q_hc = low_order_VDI.reducedOrderModelVDI(houseData, weatherTemperature, solarRad_win_in,
                                   equalAirTemp, alphaRad, ventRate, Q_ig, source_igRad, krad,
                                   t_set_heating, t_set_cooling, heater_limit, cooler_limit,
                                   dt=int(3600/times_per_hour))

# Compute averaged results
T_air_c = T_air - 273.15
T_air_mean = np.array([np.mean(T_air_c[i*times_per_hour:(i+1)*times_per_hour]) for i in range(24*60)])

T_air_1 = T_air_mean[0:24]
T_air_10 = T_air_mean[216:240]
T_air_60 = T_air_mean[1416:1440]

# Load reference results    
(T_air_ref_1, T_air_ref_10, T_air_ref_60) = tc.load_res("inputs/case08_res.csv")
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
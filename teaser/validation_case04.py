#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 09 17:40:12 2016

@author: tsz
"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

import low_order_VDI
import validationVDITestcases as tc

# Definition of time horizon
times_per_hour = 60
timesteps = 24 * 60 * times_per_hour # 60 days
timesteps_day = int(24 * times_per_hour)

# Zero inputs    
ventRate = np.zeros(timesteps)
solarRad_in = np.zeros((timesteps,1))
Q_ig = np.zeros(timesteps)

# Constant inputs
alphaRad = np.zeros(timesteps) + 5
equalAirTemp = np.zeros(timesteps) + 295.15 # all temperatures in K
weatherTemperature = np.zeros(timesteps) + 295.15 # in K

# Variable inputs
source_igRad = np.zeros(timesteps_day)
for q in range(int(6*timesteps_day/24), int(18*timesteps_day/24)):
    source_igRad[q] = 1000
source_igRad = np.tile(source_igRad, 60)

# Load constant house parameters
houseData = tc.get_house_data(case=4)

krad = 1

# Define set points (prevent heating or cooling!)
t_set_heating = np.zeros(timesteps)        # in Kelvin
t_set_cooling = np.zeros(timesteps) + 600  # in Kelvin

heater_limit = np.zeros(timesteps) + 1e10
cooler_limit = np.zeros(timesteps) - 1e10

# Calculate indoor air temperature
T_air, Q_hc = low_order_VDI.reducedOrderModelVDI(houseData, weatherTemperature, solarRad_in,
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
(T_air_ref_1, T_air_ref_10, T_air_ref_60) = tc.load_res("inputs/case04_res.csv")
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
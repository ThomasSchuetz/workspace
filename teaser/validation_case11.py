#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 12:17:52 2016

@author: tsz
"""
from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

import low_order_VDI
import validationVDITestcases as tc

import time

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
houseData = tc.get_house_data(case=11)

krad = 1

# Define set points
t_set = np.zeros(timesteps_day) + 273.15 + 22
for q in range(int(6*timesteps_day/24), int(18*timesteps_day/24)):
    t_set[q] = 273.15 + 27
t_set = np.tile(t_set, 60)
t_set_heating = t_set
t_set_cooling = t_set

#heater_limit = np.zeros(timesteps) + 500
#cooler_limit = np.zeros(timesteps) - 500
#
## Calculate indoor air temperature
#T_air, Q_hc = low_order_VDI.reducedOrderModelVDI(houseData, weatherTemperature, solarRad_in,
#                                   equalAirTemp, alphaRad, ventRate, Q_ig, source_igRad, krad,
#                                   t_set_heating, t_set_cooling, heater_limit, cooler_limit,
#                                   dt=int(3600/times_per_hour))

heater_limit = np.zeros((timesteps,3))
cooler_limit = np.zeros((timesteps,3))
heater_limit[:,0] = 500
cooler_limit[:,1] = -500
now = time.time()
# Calculate indoor air temperature
T_air, Q_hc, Q_iw, Q_ow = low_order_VDI.reducedOrderModelVDI(houseData, weatherTemperature, solarRad_in,
                                   equalAirTemp, alphaRad, ventRate, Q_ig, source_igRad, krad,
                                   t_set_heating, t_set_cooling, heater_limit, cooler_limit,
                                   heater_order=np.array([1,2,3]), cooler_order=np.array([2,1,3]),
                                   dt=int(3600/times_per_hour))

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
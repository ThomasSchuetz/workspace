# -*- coding: utf-8 -*-
"""
Created on Tue Dec 08 11:00:25 2015

@author: tsz
"""

from __future__ import division
import numpy as np

def yearly_curves(inputs, nc, nc_cumsum, typicalDays):
    # Construct (yearly) load curves
    # ub = upper bound, lb = lower bound
    clustered = np.zeros_like(inputs)
    for i in range(len(nc_cumsum)):
        if i == 0:
            lb = 0
        else:
            lb = nc_cumsum[i-1]
        ub = nc_cumsum[i]
        
        for j in range(len(inputs)):
            clustered[j, lb:ub] = np.tile(typicalDays[i][j], nc[i])
    
    return clustered

def scaling(inputs, typicalDays, nc):
    sums_inputs = [np.sum(inputs[j,:]) for j in range(inputs.shape[0])]
    scaled = np.array([nc[day] * typicalDays[day,:,:] 
                       for day in range(typicalDays.shape[0])])
    sums_scaled = [np.sum(scaled[:,j,:]) for j in range(inputs.shape[0])]
    scaling_factors = [sums_inputs[j] / sums_scaled[j] 
                       for j in range(inputs.shape[0])]
    scaled_typ_days = [scaling_factors[j] * typicalDays[:,j,:]
                       for j in range(inputs.shape[0])]

    return scaled_typ_days, scaling_factors

def monthly_clustering(inputsTransformed, nc):
    nc_cumsum = np.cumsum(nc)
    
    splits = []
    for demand in inputsTransformed:
        splits.append(np.split(demand, nc_cumsum[0:len(nc_cumsum)-1], axis=1))
    
    typical_days = []
    for cluster in range(len(nc)):
        temp = []
        for demand in splits:
            temp.append(np.mean(demand[cluster], axis=1))
        typical_days.append(temp)

    return np.array(typical_days)

def seasonal_clustering(inputsTransformed, nc, days_month):
    days_month_cumsum = np.cumsum(days_month)
    
    splits = []
    for demand in inputsTransformed:
        temp_split = np.split(demand, [days_month_cumsum[1], 
                                       days_month_cumsum[4], 
                                       days_month_cumsum[7], 
                                       days_month_cumsum[10]], axis=1)
        
        splits.append([np.hstack((temp_split[0], temp_split[4])), # Winter
                       temp_split[2], # Summer
                       np.hstack((temp_split[1], temp_split[3]))]) # Transition

    typical_days = []
    for cluster in range(len(nc)):
        temp = []
        for demand in splits:
            temp.append(np.mean(demand[cluster], axis=1))
        typical_days.append(temp)

    return np.array(typical_days)
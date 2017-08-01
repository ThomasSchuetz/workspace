#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 01 11:28:27 2015

@author: tsz
"""

from __future__ import division
import numpy as np
import math
import k_medoids
import clustering_helpers

def _distances(values, norm=2):
    """
    Compute distance matrix for all data sets (rows of values)
    
    Parameters
    ----------
    values : 2-dimensional array
        Rows represent days and columns values
    norm : integer, optional
        Compute the distance according to this norm. 2 is the standard
        Euklidean-norm.
    
    Return
    ------
    d : 2-dimensional array
        Distances between each data set
    """
    # Initialize distance matrix
    d = np.zeros((values.shape[1], values.shape[1]))

    # Define a function that computes the distance between two days
    dist = (lambda day1, day2, r: 
            math.pow(np.sum(np.power(np.abs(day1 - day2), r)), 1/r))

    # Remember: The d matrix is symmetrical!
    for i in range(values.shape[1]): # loop over first days
        for j in range(i+1, values.shape[1]): # loop second days
            d[i, j] = dist(values[:,i], values[:,j], norm)
    
    # Fill the remaining entries
    d = d + d.T
    
    return d


def cluster(inputs, 
            number_clusters=12, 
            norm=2, 
            time_limit=300,
            mip_gap=0.0):
    """
    Cluster a set of inputs into clusters by solving a k-medoid problem.
    
    Parameters
    ----------
    inputs : 2-dimensional array
        First dimension: Number of different input types.
        Second dimension: Values for each time step of interes.
    number_clusters : integer, optional
        How many clusters shall be computed?
    norm : integer, optional
        Compute the distance according to this norm. 2 is the standard
        Euklidean-norm.
    time_limit : integer, optional
        Time limit for the optimization in seconds
    mip_gap : float, optional
        Optimality tolerance (0: proven global optimum)
    
    Returns
    -------
    scaled_typ_days : 
        Scaled typical demand days. The scaling is based on the annual demands.
    nc : array_like
        Weighting factors of each cluster
    z : 2-dimensional array
        Mapping of each day to the clusters
    """
    # Determine time steps per day
    len_day = 24
    days = int(inputs.shape[1] / len_day)
    
    # Manipulate inputs
    # Initialize arrays
    inputsTransformed = []
    inputsScaled = []
    inputsScaledTransformed = []
    
    # Fill and reshape
    # Scaling to values between 0 and 1, thus all inputs shall have the same
    # weight and will be clustered equally in terms of quality 
    for i in range(inputs.shape[0]):
        vals = inputs[i,:]
        temp = (vals - np.min(vals)) / (np.max(vals) - np.min(vals))
        inputsScaled.append(temp)
        inputsScaledTransformed.append(temp.reshape((len_day, days), order="F"))
        inputsTransformed.append(vals.reshape((len_day, days), order="F"))

    # Put the scaled and reshaped inputs together
    L = np.concatenate(tuple(inputsScaledTransformed))

    # Compute distances
    d = _distances(L, norm)

    # Execute optimization model
    opt_function = k_medoids.k_medoids
    (y, z, obj, times, gap) = opt_function(d, number_clusters, time_limit, mip_gap)
    
    # Section 2.3 and retain typical days
    nc = np.zeros(number_clusters)
    typicalDays = []

    # nc contains how many days are there in each cluster
    nc = []
    for i in xrange(len(y)):
        temp = np.sum(z[i,:])
        if temp > 0:
            nc.append(temp)
            typicalDays.append([ins[:,i] for ins in inputsTransformed])

    typicalDays = np.array(typicalDays)
    nc = np.array(nc, dtype="int")

    # Scaling to preserve original demands
    scaled = clustering_helpers.scaling(inputs=inputs,
                                        typicalDays=typicalDays,
                                        nc=nc)
    (scaled_typ_days, scaling_factors) = scaled
    
    return (scaled_typ_days, nc, z, scaling_factors, times, obj, gap)
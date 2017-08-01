# -*- coding: utf-8 -*-
"""
Created on Tue Dec 08 10:44:08 2015

@author: tsz
"""

from __future__ import division

import clustering_exact
import time
import numpy as np

def cluster(inputs,
            inputs_additional=[],
            method="medoid", 
            number_clusters=12,
            reset_randomness=True,
            norm=2,
            time_limit=300, 
            mip_gap=0.0):
    """
    Parameters
    ----------
    inputs : 2-dimensional array
        First dimension: Number of different input types, e.g. number of features: electricity, heat, solar irradiance.
        Second dimension: Values for each time step of interest.
    inputs_additional : 2-dimensional array
        First dimension: Number of different input types.
        Second dimension: Values for each time step of interest.
        These inputs are not part of the clustering, but their values will be 
        returned in the typical days. These inputs basically have a weighting 
        factor of 0 within the clustering.
    method : string, optional
        - ``reference``: No clustering (reference solution)
        - ``medoid``: k-medoids clustering
        - ``center``: k-center clustering
        - ``median``: k-median clustering
        - ``mean``: k-mean clustering
        - ``monthly``: Typical days are monthly averaged days
    number_clusters : integer, optional
        How many clusters shall be computed?
    reset_randomness : Boolean, optional
        Reset the random seed (True) for reproducibility?
    norm : integer, optional
        Compute the distance according to this norm. 2 is the standard Euklidean-norm.
    time_limit : integer, optional
        Time limit for the optimization in seconds
    mip_gap : float, optional
        Optimality tolerance (0: proven global optimum)
    """
    time_before_clustering = time.time()
    
    # Determine time steps per day
    len_day = 24
    days = int(inputs.shape[1] / len_day)

    typ_inputs_add = []
    inputs_add_trans = []
    
    # Arrange additional inputs in matrix: first dim. - timesteps, second dim. - days
    for i in range(inputs_additional.shape[0]):
        vals = inputs_additional[i,:]
        inputs_add_trans.append(vals.reshape((len_day, days), order="F"))    

###############################################################################    
    medoids = True
    # Compute clusters
    results = clustering_exact.cluster(inputs=inputs,
                                       medoids=medoids,
                                       number_clusters=number_clusters,
                                       norm=norm,
                                       time_limit=time_limit,
                                       mip_gap=mip_gap)
    (scaled_typ_days, nc, z, scaling_factors, times, obj, gap) = results

    # Determine additional input values for the chosen (scaled) typdays
    for ins in inputs_add_trans:
        temp = []
        for i in range(days):
            if np.sum(z[i,:]) > 0:
                temp.append(ins[:,i])
        typ_inputs_add.append(np.array(temp))

###############################################################################        
    time_after_clustering = time.time()
    times["clustering_total"] = time_after_clustering - time_before_clustering
    
    return (scaled_typ_days,
            typ_inputs_add,
            nc, 
            scaling_factors,
            z,
            times, obj, gap)
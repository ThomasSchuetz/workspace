# -*- coding: utf-8 -*-
"""
Created on Tue Aug 01 09:14:12 2017

@author: tsz
"""

from __future__ import division
import gurobipy as gp
import numpy as np
import time

# Implementation of the k-medoids problem, as it is applied in 
# Selection of typical demand days for CHP optimization
# Fernando Domínguez-Muñoz, José M. Cejudo-López, Antonio Carrillo-Andrés and
# Manuel Gallardo-Salazar
# Energy and Buildings. Vol 43, Issue 11 (November 2011), pp. 3036-3043

# Original formulation (hereafter referred to as [1]) can be found in:

# Integer Programming and the Theory of Grouping
# Hrishikesh D. Vinod
# Journal of the American Statistical Association. Vol. 64, No. 326 (June 1969)
# pp. 506-519
# Stable URL: http://www.jstor.org/stable/2283635

def k_medoids(distances, number_clusters, timelimit=100, mipgap=0.0001):
    """
    Parameters
    ----------
    distances : 2d array
        Distances between each pair of node points. `distances` is a 
        symmetrical matrix (dissimmilarity matrix).
    number_clusters : integer
        Given number of clusters.
    timelimit : integer
        Maximum time limit for the optimization.
    """
    
    # Distances is a symmetrical matrix, extract its length
    length = distances.shape[0]
    
    # Create model
    model = gp.Model("k-Medoids-Problem")
    
    time_begin = time.time()
    
    # Create variables
    # Binary variables that are 1 if node i is assigned to cluster j    
    x = model.addVars(length, length, vtype="B", name="x")
    # Binary variables that are 1 if node j is chosen as a cluster
    y = model.addVars(length, vtype="B", name="y")
    
    # Update to introduce the variables to the model
    model.update()
    
    # Set objective - equation 2.1, page 509, [1]
    obj = gp.quicksum(distances[i,j] * x[i,j]
                      for i in xrange(length)
                      for j in xrange(length))
    model.setObjective(obj, gp.GRB.MINIMIZE)
    
    # s.t.
    # Assign all nodes to clusters - equation 2.2, page 509, [1]
    # => x_i cannot be put in more than one group at the same time
    model.addConstrs((sum(x[i,j] for j in xrange(length)) == 1 
                                     for i in xrange(length)), name="assign_x")

    # Maximum number of clusters - equation 2.3, page 509, [1]
    model.addConstr(sum(y[j] for j in xrange(length)) == number_clusters)
    
    # Prevent assigning without opening a cluster - equation 2.4, page 509, [1]
    model.addConstrs((x[i,j] <= y[j] for j in xrange(length) 
                                     for i in xrange(length)), name="x < y")
            
    model.addConstrs((x[j,j] >= y[j] for j in xrange(length)), name="md > y")
            
    # Sum of main diagonal has to be equal to the number of clusters:
    model.addConstr(sum(x[j,j] for j in xrange(length)) == number_clusters)
    
    # Set solver parameters
    model.Params.TimeLimit = timelimit
    model.Params.MIPGap = mipgap    
    
    model.update()
    time_before_opt = time.time()
    
    # Solve the model
    model.optimize()
    
    time_after_opt = time.time()
    
    # Get results
    r_x = np.array([[x[i,j].X for j in range(length)] 
                              for i in range(length)])

    r_y = np.array([y[j].X for j in xrange(length)])

    r_obj = model.ObjVal
    
    time_final = time.time()
    
    times = {"model_gen": time_before_opt - time_begin,
             "opti": time_after_opt - time_before_opt,
             "read_res": time_final - time_after_opt}
    
    return (r_y, r_x.T, r_obj, times, model.MIPGap)
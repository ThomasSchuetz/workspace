#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 07:58:05 2017

@author: Thomas
"""
from __future__ import division
from ortools.linear_solver import pywraplp
import numpy as np

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
    
    # Use CBC
    solver = pywraplp.Solver('k_medoids', 
                             pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)
    
    # Create variables
    x = {} # Binary variables that are 1 if node i is assigned to cluster j
    y = {} # Binary variables that are 1 if node j is chosen as a cluster
    for j in xrange(length):
        y[j] = solver.IntVar(0, 1, "y_"+str(j))
        
        for i in xrange(length):
            x[i,j] = solver.IntVar(0, 1, "x_"+str(i)+"_"+str(j))
    
    # Set objective - equation 2.1, page 509, [1]
    obj = solver.Sum([distances[i,j] * x[i,j]
                          for i in xrange(length)
                          for j in xrange(length)])
    
    # s.t.
    # Assign all nodes to clusters - equation 2.2, page 509, [1]
    # => x_i cannot be put in more than one group at the same time
    for i in xrange(length):
        solver.Add(solver.Sum([x[i,j] for j in xrange(length)]) == 1)
    
    # Maximum number of clusters - equation 2.3, page 509, [1]
    solver.Add(solver.Sum([y[j] for j in xrange(length)]) == number_clusters)
    
    # Prevent assigning without opening a cluster - equation 2.4, page 509, [1]
    for i in xrange(length):
        for j in xrange(length):
            solver.Add(x[i,j] <= y[j])
            
    for j in xrange(length):
        solver.Add(x[j,j] >= y[j])
            
    # Sum of main diagonal has to be equal to the number of clusters:
    solver.Add(solver.Sum([x[j,j] for j in xrange(length)]) == number_clusters)
    
    # Set solver parameters
#    params = pywraplp.MPSolverParameters()
#    params.RELATIVE_MIP_GAP = mipgap
    solver.SetTimeLimit(timelimit*1000)

    # Solve the model
#    objective = solver.Minimize(obj)
    solver.Minimize(obj)
#    solver.Solve(params)
    solver.EnableOutput()
    solver.Solve()
    
    
    # Get results
    r_x = np.array([[x[i,j].SolutionValue() for j in range(length)] 
                                            for i in range(length)])

    r_y = np.array([y[j].SolutionValue() for j in xrange(length)])

    r_obj = solver.Objective().Value()
    
    return (r_y, r_x.T, r_obj, 0)#, model.MIPGap)


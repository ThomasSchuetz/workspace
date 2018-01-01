#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 23 11:13:11 2017

@author: Thomas
"""

from __future__ import division

import gurobipy as gp
import numpy as np
import xlrd
import math

def load_traces(filename, sheetname, sampling=60):
    """
    sampling : secondwise averaging --> 60: minutewise, 3600: hourly
    """
    book = xlrd.open_workbook(filename)
    sheet = book.sheet_by_name(sheetname)
    
    res = {}
    for row in range(1, sheet.nrows):
        key = sheet.cell_value(row, 0)
        res[key] = {"number_cores": int(sheet.cell_value(row,1)),
                    "time_begin": int(math.ceil(sheet.cell_value(row,2)/sampling)),
                    "time_end": int(math.ceil(sheet.cell_value(row,3)/sampling)) + 10, # add flexibility
                    "time": int(math.ceil(sheet.cell_value(row,4)/sampling))}
    
    return res

def extract_relevant_traces(traces, time_end, time_start=0):
    res = {}
    for trace in traces:
        trace_begin = traces[trace]["time_begin"]
        if trace_begin >= time_start and trace_begin <= time_end:
            res[trace] = traces[trace]
    return res

def extract_possible_jobs(jobs, current_time):
    res = {}
    for job in jobs:
        job_begin = jobs[job]["time_begin"]
        job_end = jobs[job]["time_end"]
        if job_begin <= current_time and current_time <= job_end:
            res[job] = jobs[job]
            res[job]["rel_time"] = min(current_time-job_begin, jobs[job]["time"]) + 1
    return res

values = load_traces("traces.xlsx", "DS2")
jobs = extract_relevant_traces(values, time_end=200, time_start=0)

cores_per_node = 12
nodes = 7
power_active_node = 100
power_active_core = 5

time_steps = 200

model = gp.Model()

# Active nodes and used cores per node
node_active = {}
cores_active = {}
for n in range(nodes):
    for t in range(time_steps):
        node_active[n,t] = model.addVar(vtype="B", name="node_active_"+str(n)+"_"+str(t))
        cores_active[n,t] = model.addVar(vtype="I", name="cores_active_"+str(n)+"_"+str(t))

# Assignment of job X to node/cores
assign_jobs_nodes = {}
for n in range(nodes):
    for j in jobs.keys():
        for t in range(time_steps):
            assign_jobs_nodes[n,j,t] = model.addVar(vtype="B", name="assign_"+str(n)+"_"+j+"_"+str(t))

model.update()

obj_power_active = power_active_node * sum(sum(node_active[n,t] for n in range(nodes)) for t in range(time_steps))
obj_power_cores = power_active_core * sum(sum(cores_active[n,t] for n in range(nodes)) for t in range(time_steps))

model.setObjective(obj_power_active + obj_power_cores, gp.GRB.MINIMIZE)

for j in jobs.keys():
    time_begin = jobs[j]["time_begin"]
    time_end = jobs[j]["time_end"]
    time = jobs[j]["time"]
    
    model.addConstr(sum(sum(assign_jobs_nodes[n,j,t] for n in range(nodes)) for t in range(time_begin, time_end-time)) == 1)
    
for t in range(time_steps):
    for n in range(nodes):
        model.addConstr(cores_active[n,t] <= cores_per_node * node_active[n,t])
        
        relevant_jobs = extract_possible_jobs(jobs, t)
        
        jobs_cores = sum(sum(assign_jobs_nodes[n,j,t-tau] for tau in range(relevant_jobs[j]["rel_time"])) * 
                         jobs[j]["number_cores"] for j in relevant_jobs.keys())
        
        model.addConstr(cores_active[n,t] == sum(assign_jobs_nodes[n,j,t] * jobs[j]["number_cores"]))

model.optimize()

res_nodes_active = np.array([[node_active[n,t].X for t in range(time_steps)] for n in range(nodes)])
res_cores_active = np.array([[cores_active[n,t].X for t in range(time_steps)] for n in range(nodes)])
    
res_assign_jobs_nodes = np.array([[[assign_jobs_nodes[n,j,t].X for t in range(time_steps)] 
                                   for j in jobs.keys()] 
                                  for n in range(nodes)])


























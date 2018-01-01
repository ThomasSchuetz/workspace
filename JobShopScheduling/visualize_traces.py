#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 18:54:02 2017

@author: Thomas
"""

from __future__ import division
import xlrd
import numpy as np
import matplotlib.pyplot as plt
import math

def load_traces(filename, sheetname):
    book = xlrd.open_workbook(filename)
    sheet = book.sheet_by_name(sheetname)
    
    res = {}
    for row in range(1, sheet.nrows):
        key = sheet.cell_value(row, 0)
        res[key] = {"number_cores": sheet.cell_value(row,1),
                    "time_begin": math.ceil(sheet.cell_value(row,2)/60),
                    "time_end": math.ceil(sheet.cell_value(row,3)/60),
                    "time": math.ceil(sheet.cell_value(row,4)/60)}
    
    return res

def plot_sum_cores(values, filename=""):
    max_time = int(max([values[key]["time_end"] for key in values.keys()]))

    utilization = np.zeros(max_time)
    
    for key in values.keys():
        t_begin = values[key]["time_begin"]
        t_end = values[key]["time_end"]
        utilization[t_begin:t_end] += values[key]["number_cores"]
    
    plt.figure()
    plt.plot(np.array(range(max_time))/60, utilization)
    plt.xlabel("Time in h")
    plt.savefig("plot_traces/"+filename, bbox_inches="tight", dpi=400)

for i in range(1,7):
    plot_sum_cores(load_traces("traces.xlsx", "DS"+str(i)), "DS"+str(i)+".png")    


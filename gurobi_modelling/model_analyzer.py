# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 13:38:53 2016

@author: tsz
"""

from __future__ import division

import numpy as np
import math

def _is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def find_extreme(filename="model.lp", 
                 direction="min", 
                 section="equations", 
                 reference=None):
    """
    Parse a given lp file to find the location of the min/max values.
    
    See here for more information on the lp file format:
        https://www.gurobi.com/documentation/6.5/refman/lp_format.html
    
    Parameters
    ----------
    filename : string
        Path to the lp file
    direction : string
        'min' for the minimum or any other keyword for the maximum
    section : string
        - 'objective' : First section starting with `Maximize` or `Minimize`
        - 'equations' : Matrix section initiated with `Subject To`
        - 'rhs' : Right hand side of the equations section
        - 'bounds' : Bounds section initiated with `Bounds`
    reference : float
        - For min: find the minimum value that is greater than ``reference``
        - For max: find the minimum value that is less than ``reference``
    """
    
    # Read full text file
    with open(filename, "r") as fin:
        text = fin.readlines()
    
    # Initialize position and row's name
    position = 0
    name = "section does not provide name or does not exist"
    
    # Set dummy values for best_value
    if direction == "min":
        best_value = np.inf
    else:
        best_value = -np.inf
        
    if reference == None:
        reference = - best_value
    
    # Define next sections to look for and set initial position (i)
    if section == "objective":
        next_section = ["Subject To\n", "Bounds\n", "Generals\n", "End\n"]
        i = 0
    elif section in ("equations", "rhs"):
        next_section = ["Bounds\n", "Generals\n", "End\n"]
        i = 1
        while i < len(text) and not text[i] == "Subject To\n":
            i += 1
    else:
        next_section = ["Generals\n", "End\n"]
        i = 1
        while i < len(text) and not text[i] == "Bounds\n":
            i += 1
    
    i += 1
    
    while i < len(text) and not text[i] in next_section:
        line = text[i]
#        print line
#        (line_name, line_eq) = line.split(": ")
#        (eq, rhs) = line_eq.split("=")
        coeffs = line.split(" ")
        
        if section == "equations":
            for coeff in coeffs[:-1]:
                if _is_number(coeff):
                    val = abs(float(coeff))
                    if ((direction == "min" and val < best_value and val > reference) or
                        (direction == "max" and val > best_value and val < reference)):
                        best_value = val
                        position = i+1
#                    name = line_name
        elif section == "rhs":
            if _is_number(coeffs[-1]):
                val = abs(float(coeffs[-1]))
                if ((direction == "min" and val < best_value and val > reference) or
                    (direction == "max" and val > best_value and val < reference)):
                    best_value = val
                    position = i+1
        i += 1
    
    
    return (position, name, best_value)

def find_number(filename="model.lp",
                number=10,
                direction="min",
                section="equations",
                reference=None):
    res = []
    ref = reference
    for i in range(number):
        (pos, name, ref) = find_extreme(filename, direction, section, ref)
        res.append((pos, ref))
    
    return res

#res_min = find_extreme(filename="model.lp", direction="min", section="equations", reference=1e-6)
#res_max = find_extreme(filename="model.lp", direction="max", section="equations")

res_min_rhs = find_number(number=20, section="rhs", direction="min")
res_max_rhs = find_number(number=20, section="rhs", direction="max")
res_min_eq = find_number(number=20, section="equations", direction="min")
res_max_eq = find_number(number=20, section="equations", direction="max")

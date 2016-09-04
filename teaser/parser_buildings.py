#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 11:30:43 2015

@author: tsz
"""

import re

def _parse_line(line):
    """
    """
    # Remove white spaces, quotation marks and line endings
    remove_strings = (" ", '"', ",\n", ");")
    for s in remove_strings:
        line = line.replace(s, "")
    
    # Split line at equal sign
    (key, value) = line.split("=")
    
    # Analyse "value"
    if value[0] == "{":     # Handle arrays
        # Remove brackets and split at comma
        value = value.replace("{","").replace("}","")
        value = value.split(",")
        # Transform to floats instead of strings
        value = [float(v) for v in value]
    elif value == "true":   # Handle "true" and "false"
        value = True
    elif value == "false":
        value = False
    else:                   # Handle other values
        try:    # Transfer to float
            value = float(value)
        except: # Or skip
            pass
    
    
    return (key, value)

def parse_record(filename):
    """
    This function parses house related data into a dictionary
    """

    # Initialize content
    content = {}
    
    # Parse all lines 
    with open(filename) as f:
        # Read all lines at once
        all_lines = f.readlines()
        
        # Skip the first three lines as well as the last one:
        for i in range(3, len(all_lines)-1):
            # Skip empty lines
            if all_lines[i].isspace():
                continue
            # Parse current line
            (key, value) = _parse_line(all_lines[i])
            # Save results
            content[key] = value

    # Return parsed content
    return content


# Example
if __name__ == "__main__":
    filename = "LCS_45_1.mo"

    content = parse_record(filename)
    
    for key in content.keys():
        print(key+ ": " + str(content[key]))
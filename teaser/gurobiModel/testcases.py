#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 14:19:07 2016

@author: Markus
"""
from __future__ import division
import numpy as np

#%% Load a standard result file
def load_res(filename):
    res = np.loadtxt(filename, delimiter=",", skiprows=1) # Skip time step 0
    
    # ignore time
    result = res[:,1:res.shape[1]]
    
    day1 = result[0:24, :]
    day2 = result[24:48, :]
    day3 = result[48:72, :]
    
    return (day1, day2, day3)

#%% Common house inputs
def get_house_data(timesteps, case=1):
    if case in (1,2):
        return {# Thermal zone
                "Vair": 0,
                "alphaRad": np.zeros(timesteps)+5,
                "nOrientations": 1,
        
                # Windows
                "AWin": np.zeros(1),
                "ATransparent": np.zeros(1),
                "alphaWin": 2.7,
                "RWin": 0.00000001,
                "gWin": 1,
                "ratioWinConRad": 0,
                "indoorPortWin": False,
                
                # Exterior Walls
                "AExt": np.array([10.5]),
                "alphaExt": 2.7,
                "nExt": 1,
#                "RExt": np.array([0.00436791293674]),
                "RExt": 0.00436791293674,
                "RExtRem": 0.03895919557,
#                "CExt": np.array([1600848.94]),
                "CExt": 1600848.94,
                "indoorPortExtWalls": False,
                
                # Interior walls
                "AInt": 75.5,
                "alphaInt": 2.24,
                "nInt": 1,
#                "RInt": np.array([0.000595693407511]),
                "RInt": 0.000595693407511,
#                "CInt": np.array([14836354.6282]),
                "CInt": 14836354.6282,
                "indoorPortIntWalls": False,
                
                # theConWall
                "alphaWall": 25*10.5
                }
    elif case in (3,4):
        return {# Thermal zone
                "Vair": 0,
                "alphaRad": np.zeros(timesteps)+5,
                "nOrientations": 1,
        
                # Windows
                "AWin": np.zeros(1),
                "ATransparent": np.zeros(1),
                "alphaWin": 2.7,
                "RWin": 0.00000001,
                "gWin": 1,
                "ratioWinConRad": 0,
                "indoorPortWin": False,
                
                # Exterior Walls
                "AExt": np.array([10.5]),
                "alphaExt": 2.7,
                "nExt": 1,
#                "RExt": np.array([0.00404935160802]),
                "RExt": 0.00404935160802,
                "RExtRem": 0.039330865,
#                "CExt": np.array([47900]),
                "CExt": 47900,
                "indoorPortExtWalls": False,
                
                # Interior walls
                "AInt": 75.5,
                "alphaInt": 2.24,
                "nInt": 1,
#                "RInt": np.array([0.003237138]),
                "RInt": 0.003237138,
#                "CInt": np.array([7297100]),
                "CInt": 7297100,
                "indoorPortIntWalls": False,
                
                # theConWall
                "alphaWall": 25*10.5
                }
    elif case in (5,):
        return {# Thermal zone
                "Vair": 0,
                "alphaRad": np.zeros(timesteps)+5,
                "nOrientations": 1,
        
                # Windows
                "AWin": np.zeros(1),
                "ATransparent": np.array([7]),
                "alphaWin": 2.7,
                "RWin": 0.00000001,
                "gWin": 1,
                "ratioWinConRad": 0.09,
                "indoorPortWin": False,
                
                # Exterior Walls
                "AExt": np.array([10.5]),
                "alphaExt": 2.7,
                "nExt": 1,
#                "RExt": np.array([0.00436791293674]),
                "RExt": 0.00436791293674,
                "RExtRem": 0.03895919557,
#                "CExt": np.array([1600848.94]),
                "CExt": 1600848.94,
                "indoorPortExtWalls": False,
                
                # Interior walls
                "AInt": 75.5,
                "alphaInt": 2.24,
                "nInt": 1,
#                "RInt": np.array([0.000595693407511]),
                "RInt": 0.000595693407511,
#                "CInt": np.array([14836354.6282]),
                "CInt": 14836354.6282,
                "indoorPortIntWalls": False,
                
                # theConWall
                "alphaWall": 25*10.5
                }
    elif case in (6,):
        return {# Thermal zone
                "Vair": 0,
                "alphaRad": np.zeros(timesteps)+5,
                "nOrientations": 1,
        
                # Windows
                "AWin": np.zeros(1),
                "ATransparent": np.zeros(1),
                "alphaWin": 2.7,
                "RWin": 0.00000001,
                "gWin": 1,
                "ratioWinConRad": 0,
                "indoorPortWin": False,
                
                # Exterior Walls
                "AExt": np.array([10.5]),
                "alphaExt": 2.7,
                "nExt": 1,
#                "RExt": np.array([0.004367913]),
                "RExt": 0.004367913,
                "RExtRem": 0.03895917,
#                "CExt": np.array([1600800]),
                "CExt": 1600800,
                "indoorPortExtWalls": False,
                
                # Interior walls
                "AInt": 75.5,
                "alphaInt": 2.24,
                "nInt": 1,
#                "RInt": np.array([0.000595515]),
                "RInt": 0.000595515,
#                "CInt": np.array([14836200]),
                "CInt": 14836200,
                "indoorPortIntWalls": False,
                
                # theConWall
                "alphaWall": 25*10.5
                }
    elif case in (7,):
        return {# Thermal zone
                "Vair": 0,
                "alphaRad": np.zeros(timesteps)+5,
                "nOrientations": 1,
        
                # Windows
                "AWin": np.zeros(1),
                "ATransparent": np.zeros(1),
                "alphaWin": 2.7,
                "RWin": 0.00000001,
                "gWin": 1,
                "ratioWinConRad": 0,
                "indoorPortWin": False,
                
                # Exterior Walls
                "AExt": np.array([10.5]),
                "alphaExt": 2.7,
                "nExt": 1,
#                "RExt": np.array([0.00436791293674]),
                "RExt": 0.00436791293674,
                "RExtRem": 0.03895919557,
#                "CExt": np.array([1600848.94]),
                "CExt": 1600848.94,
                "indoorPortExtWalls": False,
                
                # Interior walls
                "AInt": 75.5,
                "alphaInt": 2.24,
                "nInt": 1,
#                "RInt": np.array([0.000595693407511]),
                "RInt": 0.000595693407511,
#                "CInt": np.array([14836354.6282]),
                "CInt": 14836354.6282,
                "indoorPortIntWalls": False,
                
                # theConWall
                "alphaWall": 25*10.5
                }
    elif case in (8,9):
        return {# Thermal zone
                "Vair": 0,
                "alphaRad": np.zeros(timesteps)+5,
                "nOrientations": 2,
        
                # Windows
                "AWin": np.zeros(2),
                "ATransparent": np.array([7,7]),
                "alphaWin": 2.7,
                "RWin": 0.00000001,
                "gWin": 1,
                "ratioWinConRad": 0.09,
                "indoorPortWin": False,
                
                # Exterior Walls
                "AExt": np.array([10.5,15]),
                "alphaExt": 2.7,
                "nExt": 1,
#                "RExt": np.array([0.0017362530106]),
                "RExt": 0.0017362530106,
                "RExtRem": 0.01913729904,
#                "CExt": np.array([5259932.23]),
                "CExt": 5259932.23,
                "indoorPortExtWalls": False,
                
                # Interior walls
                "AInt": 60.5,
                "alphaInt": 2.12,
                "nInt": 1,
#                "RInt": np.array([0.000668895639141]),
                "RInt": 0.000668895639141,
#                "CInt": np.array([12391363.8631]),
                "CInt": 12391363.8631,
                "indoorPortIntWalls": False,
                
                # theConWall
                "alphaWall": 25 * 25.5
                }
    elif case in (10,):
        return {# Thermal zone
                "Vair": 0,
                "alphaRad": np.zeros(timesteps)+5,
                "nOrientations": 1,
        
                # Windows
                "AWin": np.zeros(1),
                "ATransparent": np.array([7]),
                "alphaWin": 2.7,
                "RWin": 0.00000001,
                "gWin": 1,
                "ratioWinConRad": 0.09,
                "indoorPortWin": False,
                
                # Exterior Walls
                "AExt": np.array([28]),
                "alphaExt": 2.4,
                "nExt": 1,
#                "RExt": np.array([0.00171957697767797]),
                "RExt": 0.00171957697767797,
                "RExtRem": 0.011638548,
#                "CExt": np.array([4338751.41]),
                "CExt": 4338751.41,
                "indoorPortExtWalls": False,
                
                # Interior walls
                "AInt": 58,
                "alphaInt": 2.398,
                "nInt": 1,
#                "RInt": np.array([0.000779671554640369]),
                "RInt": 0.000779671554640369,
#                "CInt": np.array([12333949.4129606]),
                "CInt": 12333949.4129606,
                "indoorPortIntWalls": False,
                
                # theConWall
                "alphaWall": 9.75 * 28
                }
    elif case in (11,):
        return {# Thermal zone
                "Vair": 0,
                "alphaRad": np.zeros(timesteps)+5,
                "nOrientations": 1,
        
                # Windows
                "AWin": np.zeros(1),
                "ATransparent": np.zeros(1),
                "alphaWin": 2.7,
                "RWin": 0.00000001,
                "gWin": 1,
                "ratioWinConRad": 0,
                "indoorPortWin": False,
                
                # Exterior Walls
                "AExt": np.array([10.5]),
                "alphaExt": 2.7,
                "nExt": 1,
#                "RExt": np.array([0.00436791293674]),
                "RExt": 0.00436791293674,
                "RExtRem": 0.03895919557,
#                "CExt": np.array([1600848.94]),
                "CExt": 1600848.94,
                "indoorPortExtWalls": False,
                
                # Interior walls
                "AInt": 75.5,
                "alphaInt": 3,
                "nInt": 1,
#                "RInt": np.array([0.000595693407511]),
                "RInt": 0.000595693407511,
#                "CInt": np.array([14836354.6282]),
                "CInt": 14836354.6282,
                "indoorPortIntWalls": True,
                
                # theConWall
                "alphaWall": 25 * 10.5
                }
    elif case in (12,):
        return {# Thermal zone
                "Vair": 0.1,
                "alphaRad": np.zeros(timesteps)+5,
                "nOrientations": 1,
        
                # Windows
                "AWin": np.zeros(1),
                "ATransparent": np.array([7]),
                "alphaWin": 2.7,
                "RWin": 0.00000001,
                "gWin": 1,
                "ratioWinConRad": 0.09,
                "indoorPortWin": False,
                
                # Exterior Walls
                "AExt": np.array([10.5]),
                "alphaExt": 2.7,
                "nExt": 1,
#                "RExt": np.array([0.00436791293674]),
                "RExt": 0.00436791293674,
                "RExtRem": 0.03895919557,
#                "CExt": np.array([1600848.94]),
                "CExt": 1600848.94,
                "indoorPortExtWalls": False,
                
                # Interior walls
                "AInt": 75.5,
                "alphaInt": 2.24,
                "nInt": 1,
#                "RInt": np.array([0.000595693407511]),
                "RInt": 0.000595693407511,
#                "CInt": np.array([14836354.6282]),
                "CInt": 14836354.6282,
                "indoorPortIntWalls": False,
                
                # theConWall
                "alphaWall": 25*10.5
                }

#%%
def get_eqAirTemp_params(case=8):
    if case in (8,9):
        return {"aExt": 0.7,
                "eExt": 0.9,
                "n": 2,
                "wfWall": np.array([0.05796831135677373, 0.13249899738691134]),
                "wfWin": np.array([0.4047663456281575, 0.4047663456281575]),
                "wfGro": 0,
                "T_Gro": 273.15 + 12,
                "alpha_wall_out": 20,
                "alpha_rad_wall": 5,
                "withLongwave": False}
    elif case in (10,):
        return {"aExt": 0.7,
                "eExt": 0.9,
                "n": 1,
                "wfWall": np.array([0.04646093176283288,]),
                "wfWin": np.array([0.32441554918476245,]),
                "wfGro": 0.6291235190524047,
                "T_Gro": 273.15 + 15,
                "alpha_wall_out": 20,
                "alpha_rad_wall": 5,
                "withLongwave": False}
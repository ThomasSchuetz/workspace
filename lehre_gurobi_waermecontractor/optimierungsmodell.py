# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 18:27:25 2017

@author: tsz
"""

from __future__ import division
import gurobipy as gp
qsum = gp.quicksum

def opt(waermebedarf, par_geraete, par_oekonomie, par_sonst):
    """
    """
    
    # Abkuerzungen einfuehren
    dt = par_sonst["zeitdiskretisierung"]
    anz_zeitschritte = par_sonst["anz_zeitschritte"]
    anz_geraete = {"KWK": len(par_geraete["KWK"].keys()),
                   "Kessel": len(par_geraete["Kessel"].keys())}
    
    
    indices_geraete = [(geraet, nummer) 
                       for nummer in range(1, 1 + anz_geraete[geraet]) 
                       for geraet in ("KWK", "Kessel")]

    indices_geraete_zeit = [(geraet, nummer, zeit)
                            for zeit in range(anz_zeitschritte)
                            for nummer in range(1, 1 + anz_geraete[geraet]) 
                            for geraet in ("KWK", "Kessel")]
    
    indices_kwk_zeit = [("KWK", nummer, zeit)
                        for zeit in range(anz_zeitschritte)
                        for nummer in range(1, 1 + anz_geraete["KWK"])]
    
    model = gp.Model("Auslegung Waermeversorgung")
    
    x = model.addVars(indices_geraete, vtype="B", name="x")
    y = model.addVars(indices_geraete_zeit, vtype="B", name="y")
    
    p = model.addVars(indices_kwk_zeit, vtype="C", name="P")
    q = model.addVars(indices_geraete_zeit, vtype="C", name="Q")
    
    inv = model.addVar(vtype="C", name="Investment")
    jz = model.addVar(vtype="C", name="Jaehrliche_Zahlung")
    
    model.update()
    
    model.setObjective(inv + par_oekonomie["RBF"] * jz, gp.GRB.MINIMIZE)
    
    model.addConstr(inv == qsum(x[geraet,nummer] * par_geraete[geraet][nummer]["c_inv"] 
                                for nummer in range(1, 1 + anz_g[geraet]) 
                                for geraet in ("KWK", "Kessel")), 
                    "Investitionen")
    
    gas_kwk = dt * qsum(qsum((p["KWK", nummer, zeit] + q["KWK", nummer, zeit]) 
                             for zeit in range(anz_zeitschritte)) 
                        / par_geraete["KWK"][nummer]["omega"] 
                        for nummer in range(1, 1 + anz_geraete["KWK"]))
    gas_kessel = dt * qsum(qsum(p["Kessel", nummer, zeit] for zeit in range(anz_zeitschritte)) 
                           / par_geraete["Kessel"][nummer]["eta"] 
                           for nummer in range(1, 1 + anz_geraete["Kessel"]))
    strom_kwk = dt * qsum(qsum(q["KWK", nummer, zeit] for zeit in range(anz_zeitschritte))
                          for nummer in range(1, 1 + anz_geraete["KWK"]))
    
    model.addConstr(jz == par_oekonomie["c_gas"] * (gas_kwk + gas_kessel) - 
                          par_oekonomie["p_el"] * strom_kwk,
                    "JaehrlicheZahlungen")
    
    model.addConstrs( (y[geraet, nummer, zeit] <= x[geraet, nummer] 
                       for zeit in range(anz_zeitschritte)
                       for nummer in range(1, 1 + anz_geraete[geraet]) 
                       for geraet in ("KWK", "Kessel")),
                     "Aktivierung_nur_bei_Kauf")
    
    model.addConstrs( (q[geraet, nummer, zeit] >= 
                       y[geraet, nummer, zeit] * par_geraete[geraet][nummer]["q_min"]
                       for zeit in range(anz_zeitschritte)
                       for nummer in range(1, 1 + anz_geraete[geraet]) 
                       for geraet in ("KWK", "Kessel")),
                     "Thermische_Leistung_Minimum")

    model.addConstrs( (q[geraet, nummer, zeit] <= 
                       y[geraet, nummer, zeit] * par_geraete[geraet][nummer]["q_max"]
                       for zeit in range(anz_zeitschritte)
                       for nummer in range(1, 1 + anz_geraete[geraet]) 
                       for geraet in ("KWK", "Kessel")),
                     "Thermische_Leistung_Maximum")
    
    model.addConstrs( (p["KWK", nummer, zeit] == 
                       q["KWK", nummer, zeit] * par_geraete["KWK"][nummer]["sigma"]
                       for zeit in range(anz_zeitschritte)
                       for nummer in range(1, 1 + anz_geraete["KWK"])),
                     "Elektrische_Leistung_KWK")
    
    model.addConstrs( (qsum(qsum(q[geraet, nummer, zeit] 
                                 for nummer in range(1, 1 + anz_geraete[geraet])) 
                            for geraet in ("KWK", "Kessel")) 
                       == waermebedarf[zeit]
                       for zeit in range(anz_zeitschritte)),
                     "Waermebilanz")
    
    model.optimize()
    
    
    
    
    
    
    
    
    
    
    
    
    
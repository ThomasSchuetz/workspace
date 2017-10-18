# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 12:31:38 2017

@author: tsz
"""

from __future__ import division
import xlrd

def lese_geraete(tabelle, technologie):
    """
    Diese Funktion liest aus dem Reiter "technologie" (also KWK oder Kessel) alle Parameter
    
    Parameter
    ---------
    tabelle : Geoeffnete Tabelle
    """
    
    

def lese_oekonomische_parameter(tabelle):
    """
    Diese Funktion liest aus dem Reiter "OekonomischeParameter" alle Parameter
    
    Parameter
    ---------
    tabelle : Geoeffnete Tabelle
    """
    
    # Erstelle ein neues Dictionary, in welchem die Ergebnisse gespeichert werden
    # Hinweis: cell_value(Zeile, Spalte), wobei die Nummerierung bei 0 beginnt
    ergebnis = {"i": tabelle.cell_value(1,1), 
                "T": tabelle.cell_value(2,1),
                "RBF": tabelle.cell_value(3,1),
                "c_gas": tabelle.cell_value(4,1),
                "p_el": tabelle.cell_value(5,1)}
    
    # Die eingelesenen Parameterwerte werden zurueckgegeben
    return ergebnis

# Der folgende Abschnitt ist lediglich ein Test, um zu ueberpruefen, ob alle 
# Parameter korrekt eingelesen werden.
if __name__ == "__main__":
    excel_datei = "inputs.xlsx"
    # Oeffne eine existierende Excel-Datei
    geoeffnete_datei = xlrd.open_workbook(excel_datei)

    # Finde die korrekte Tabelle innerhalb dieser Excel-Datei
    tabelle_oekonomische_parameter = geoeffnete_datei.sheet_by_name("OekonomischeParameter")
    # Lese die entsprechenden Parameter aus
    parameter_oekonomie = lese_oekonomische_parameter(tabelle_oekonomische_parameter)

    # Lese die Parameter aller verfuegbaren KWK Anlagen und Kessel
    tabelle_kwk = geoeffnete_datei.sheet_by_name("KWK")
    parameter_kwk = lese_geraete(tabelle_kwk, "KWK")
    
    tabelle_kessel = geoeffnete_datei.sheet_by_name("Kessel")
    parameter_kwk = lese_geraete(tabelle_kessel, "Kessel")
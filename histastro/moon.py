#!/bin/env python3

import math as m
import numpy as np
from histastro.constants import pi2,d2r,jd2000,moonRad


def readData(inFile):
    """Read the periodic terms for the ELP82B theory, selected by Meeus"""
    lrTerms = np.genfromtxt(inFile, delimiter=',', skip_header=1,  max_rows=60)  # Longitude and radius (6 columns: 4 args, 2 coefs)
    bTerms  = np.genfromtxt(inFile, delimiter=',', skip_header=61, max_rows=60)  # Latitude (5 columns: 4 args, 1 coef)
    return lrTerms,bTerms


def compute_lbr(jde, lrTerms,bTerms):
    tjc   = (jde-jd2000)/36525  # Julian Centuries after 2000.0 in dynamical time, the T in Meeus, p.163
    tjc2  = tjc**2
    tjc3  = tjc*tjc2
    tjc4  = tjc2**2
    
    # Moon's mean longitude, Meeus p.338:
    lm = (3.8103408236 + 8399.7091116339958*tjc - 2.755176757e-5*tjc2 + 3.239043e-8*tjc3 - 2.6771e-10*tjc4) % pi2
    
    # Delauney arguments [d, ms, mm, f] (Meeus p.338) = [D, l', l, F] in Simon et al. 1994, Sect. 3.5:
    d  = (5.1984665298  + 7771.377144834*tjc - 3.2845e-5*tjc2  + 3.197347e-8*tjc3    - 1.5436512e-10*tjc4) % pi2  # Moon's mean elongation
    ms = (6.240060127   + 628.301955167*tjc  - 2.681e-6*tjc2   + 7.1267017e-10*tjc3)                       % pi2  # Sun's mean anomaly
    mm = (2.355555637   + 8328.691424759*tjc + 1.52566e-4*tjc2 + 2.5041e-7*tjc3      - 1.18633e-9*tjc4)    % pi2  # Moon's mean anomaly
    f  = (1.627905158   + 8433.466158061*tjc - 6.3773e-5*tjc2  - 4.94988e-9*tjc3     + 2.02167e-11*tjc4)   % pi2  # Moon's argument of latitute
    DelArgs = [d,ms,mm,f]  # Delauney arguments
    
    e  = 1 - 0.002516*tjc - 0.0000074*tjc2
    
    # Compute ELP terms:
    lon=0; lat=0; dist=0
    for iLine in range(60):  # Data lines
        argl = 0
        argb = 0
        for iArg in range(4):  # Integer arguments
            argl += lrTerms[iLine,iArg] * DelArgs[iArg]
            argb +=  bTerms[iLine,iArg] * DelArgs[iArg]
        
        lon  += m.sin(argl) * lrTerms[iLine,4] * e**(abs(lrTerms[iLine,1]))
        lat  += m.sin(argb) *  bTerms[iLine,4] * e**(abs(bTerms[iLine,1]))
        dist += m.cos(argl) * lrTerms[iLine,5] * e**(abs(lrTerms[iLine,1]))
    
    
    # Meeus, p.338:
    a1 = (2.090032  +     2.301199 * tjc) % pi2
    a2 = (0.926595  + 8364.7398477 * tjc) % pi2
    a3 = (5.4707345 + 8399.6847253 * tjc) % pi2

    # Meeus, p.342:
    lon += (3958*m.sin(a1) + 1962*m.sin(lm-f) + 318*m.sin(a2)) * d2r * 1e-6
    lat += (2235*m.sin(lm) + 382*m.sin(a3) + 175*m.sin(a1-f) + 175*m.sin(a1+f) + 127*m.sin(lm-mm) - 115*m.sin(lm+mm)) * d2r * 1e-6
    
    # Compute nutation in longitude:
    omg  = 2.18243858558 - 33.7570459367*tjc + 3.6142278e-5*tjc2 + 3.87850944888e-8*tjc3   # Moon's mean lon. of asc.node, Meeus p.144
    ls   = 4.89506386655 + 62.84528862*tjc                                                 # Mean long. Sun, Meeus p.144
    dpsi = -8.338795e-5*m.sin(omg) - 6.39954e-6*m.sin(2*ls) - 1.115e-6*m.sin(2*lm) + 1.018e-6*m.sin(2*omg)
    #dpsi = 0
    
    # Add mean values:
    lon  = (lon + lm + dpsi) % pi2
    dist = dist + 385000.56  # in km
    
    
    diam = 2*m.atan(moonRad/dist)
    return lon,lat,dist,diam


#!/bin/env python3


import math as m
import histastro.vsop as vsop
import histastro.coordinates as coord
import histastro.datetime as dt

# Check obliquity:
JD = dt.julianDay(1000,1,1)
print(JD, m.degrees(coord.obliquity(JD)))
JD = dt.julianDay(-3000,1,1)
print(JD, m.degrees(coord.obliquity(JD)))
    



# Read VSOP data files for Earth and desired planet:
lonTermsE,latTermsE,radTermsE = vsop.readVSOP('data/VSOP87D.ear')
lonTerms,latTerms,radTerms = vsop.readVSOP('data/VSOP87D.jup')

# Compute heliocentric ecliptical coordinates:
for month in range(1,4):
    JD = dt.julianDay(-1000, month, 1)
    eps = coord.obliquity(JD)
    
    # Earth:
    HClonE,HClatE,HCradE = vsop.computeLBR(JD, lonTermsE,latTermsE,radTermsE)
    
    # Planet:
    HClon,HClat,HCrad = vsop.computeLBR(JD, lonTerms,latTerms,radTerms)
    
    # Compute geocentric ecliptical and equatorial coordinates:
    lon,lat,rad = vsop.hc2gc(HClonE,HClatE,HCradE, HClon,HClat,HCrad)
    ra,dec = coord.ecl2eq(lon,lat, eps)
    
    print(month,JD, lon,lat,HCrad,rad, ra,dec)
    
    


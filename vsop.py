#!/bin/env python3

import math as m
import histastro.vsop as vsop
import histastro.coordinates as coord


# Read VSOP data files for Earth and desired planet:
lonTermsE,latTermsE,radTermsE = vsop.readVSOP('VSOP87D.ear')
lonTerms,latTerms,radTerms = vsop.readVSOP('VSOP87D.jup')
#lonTerms,latTerms,radTerms = vsop.readVSOP('VSOP87D.ura')

#xTermsE,yTermsE,zTermsE = vsop.readVSOP('VSOP87C.ear')
#xTerms,yTerms,zTerms = vsop.readVSOP('VSOP87C.jup')


# Compute heliocentric ecliptical coordinates:
for iter in range(1):
    #JDE = 2451545.0 + iter*100
    dYear = iter*50
    JDE = 2451545.0 - dYear*365
    
    # Earth:
    HClonE,HClatE,HCradE = vsop.computeLBR(JDE, lonTermsE,latTermsE,radTermsE)
    #HCxE,HCyE,HCzE = vsop.computeLBR(JDE, xTermsE,yTermsE,zTermsE)
    #print(HClonE,HClatE,HCradE)
    
    # Planet:
    HClon,HClat,HCrad = vsop.computeLBR(JDE, lonTerms,latTerms,radTerms)
    #HCx,HCy,HCz = vsop.computeLBR(JDE, xTerms,yTerms,zTerms)
    #print(HClon,HClat,HCrad)
    
    
    
    # Compute geocentric ecliptical coordinates:
    lon,lat,rad = vsop.hc2gc(HClonE,HClatE,HCradE, HClon,HClat,HCrad)
    print(iter,m.degrees(lon), m.degrees(lat),rad)
    
    # Compute geocentric ecliptical coordinates:
    #lon1,lat1,rad1 = vsop.xyz_hc2lbr_gc(HCxE,HCyE,HCzE, HCx,HCy,HCz)
    #print(iter, m.degrees(lon), m.degrees(lat), m.degrees(rad))
    #print(iter, 2000-dYear, m.degrees(lon-lon1), m.degrees(lat-lat1), rad-rad1)  # ~3000 BCE: dl~3e-3d, db~3e-5d, dr~4e-4 AU
    

# Check obliquity:
JDE = 2448347.5000000000
print(iter, 2000-dYear, JDE, m.degrees(coord.obliquity(JDE)))
JDE = -1931442.5000000000
print(iter, 2000-dYear, JDE, m.degrees(coord.obliquity(JDE)))
    

# Benchmark:
# Starting script:                0.061   +- 0.003 s
# Reading 2 VSOP files:           0.364   +- 0.006 s
# Computing geocentric position:  0.00274 +- 0.00010 s

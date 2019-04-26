#!/bin/env python3

import math as m


def readVSOP(fileName):
    """Read the periodic terms for a heliocentric ecliptical planet position from a VSOP87D file"""

    inFile = open(fileName,'r')
    
    import fortranformat as ff
    formatHeader = ff.FortranRecordReader('(40x,I3, 16x,I1,I8)')         # Block header format
    formatBody   = ff.FortranRecordReader('(79x,F18.11,F14.11,F20.11)')  # Block body format
    
    lonTerms=[]; latTerms=[]; radTerms=[]
    
    for iBlock in range(3*6):  # 3 variables (l,b,r), up to 6 powers (0-5)
        line = inFile.readline()
        var,power,nTerm = formatHeader.read(line)
        #print(var,power,nTerm)
        if line == '': break  # EoF
        
        for iLine in range(nTerm):
            line = inFile.readline()
            a,b,c = formatBody.read(line)
            #print(iLine, var,power, a,b,c)
            
            if var == 1: lonTerms.append([power, a,b,c])  # var=1: ecliptic longitude
            if var == 2: latTerms.append([power, a,b,c])  # var=2: ecliptic latitude
            if var == 3: radTerms.append([power, a,b,c])  # var=3: radial distance

    return lonTerms,latTerms,radTerms


# Compute heliocentric ecliptical coordinates from periodic terms
def computeLBR(JDE, lonTerms,latTerms,radTerms):
    Tjm = (JDE - 2451545.0)/365250.0
    
    lon=0.0; lat=0.0; rad=0.0
    for terms in lonTerms:
        cosTerm = terms[1] * m.cos(terms[2] + terms[3]*Tjm)
        lon += cosTerm * Tjm**terms[0]
    for terms in latTerms:
        cosTerm = terms[1] * m.cos(terms[2] + terms[3]*Tjm)
        lat += cosTerm * Tjm**terms[0]
    for terms in radTerms:
        cosTerm = terms[1] * m.cos(terms[2] + terms[3]*Tjm)
        rad += cosTerm * Tjm**terms[0]
        
    return lon,lat,rad


# Convert from heliocentric spherical to rectangular coordinates, and take the difference (i.e., geocentric
# rectangular coordinates):
def hc2gc(l0,b0,r0, l,b,r):
    x = r * m.cos(b) * m.cos(l)  -  r0 * m.cos(b0) * m.cos(l0)
    y = r * m.cos(b) * m.sin(l)  -  r0 * m.cos(b0) * m.sin(l0)
    z = r * m.sin(b)             -  r0 * m.sin(b0)
    
    # Convert geocentric rectangular to spherical coordinates:
    if x==0 and y==0 and z==0:
        lon = 0.0
        lat = 0.0
        rad = 0.0
    else:
        x2 = x**2
        y2 = y**2
        
        lon = m.atan2(y, x)                # Longitude
        lat = m.atan2(z, m.sqrt(x2 + y2))  # Latitude
        rad = m.sqrt(x2 + y2 + z**2)       # Distance
       
    return lon,lat,rad


# Convert from heliocentric rectangular coordinates to geocentric spherical coordinates):
def xyz_hc2lbr_gc(x0,y0,z0, x,y,z):
    dx = x - x0
    dy = y - y0
    dz = z - z0
    
    # Convert geocentric rectangular to spherical coordinates:
    if dx==0 and dy==0 and dz==0:
        lon = 0.0
        lat = 0.0
        rad = 0.0
    else:
        dx2 = dx**2
        dy2 = dy**2
        
        lon = m.atan2(dy, dx)                 # Longitude
        lat = m.atan2(dz, m.sqrt(dx2 + dy2))  # Latitude
        rad = m.sqrt(dx2 + dy2 + dz**2)       # Distance
       
    return lon,lat,rad



#exit()

# Read VSOP data files for Earth and desired planet:
lonTermsE,latTermsE,radTermsE = readVSOP('VSOP87D.ear')
lonTerms,latTerms,radTerms = readVSOP('VSOP87D.jup')
#lonTerms,latTerms,radTerms = readVSOP('VSOP87D.ura')

#xTermsE,yTermsE,zTermsE = readVSOP('VSOP87C.ear')
#xTerms,yTerms,zTerms = readVSOP('VSOP87C.jup')


# Compute heliocentric ecliptical coordinates:
for iter in range(101):
    #JDE = 2451545.0 + iter*100
    dYear = iter*50
    JDE = 2451545.0 - dYear*365
    
    # Earth:
    HClonE,HClatE,HCradE = computeLBR(JDE, lonTermsE,latTermsE,radTermsE)
    #HCxE,HCyE,HCzE = computeLBR(JDE, xTermsE,yTermsE,zTermsE)
    #print(HClonE,HClatE,HCradE)
    
    # Planet:
    HClon,HClat,HCrad = computeLBR(JDE, lonTerms,latTerms,radTerms)
    #HCx,HCy,HCz = computeLBR(JDE, xTerms,yTerms,zTerms)
    #print(HClon,HClat,HCrad)
    
    
    
    # Compute geocentric ecliptical coordinates:
    lon,lat,rad = hc2gc(HClonE,HClatE,HCradE, HClon,HClat,HCrad)
    print(iter,m.degrees(lon), m.degrees(lat),rad)
    
    # Compute geocentric ecliptical coordinates:
    #lon1,lat1,rad1 = xyz_hc2lbr_gc(HCxE,HCyE,HCzE, HCx,HCy,HCz)
    #print(iter, m.degrees(lon), m.degrees(lat), m.degrees(rad))
    #print(iter, 2000-dYear, m.degrees(lon-lon1), m.degrees(lat-lat1), rad-rad1)  # ~3000 BCE: dl~3e-3°, db~3e-5°, dr~4e-4 AU
    

# Benchmark:
# Starting script:                0.061   +- 0.003 s
# Reading 2 VSOP files:           0.364   +- 0.006 s
# Computing geocentric position:  0.00274 +- 0.00010 s

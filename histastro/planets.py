"""HistAstro planet functions"""

import math as m
import numpy.core as np
import histastro.datetime as dt

pi2 = m.pi*2

def readVSOP(dataDir, pl):
    """Read the periodic terms for a heliocentric ecliptical planet position from a VSOP87D.* file, located in the
    directory specified by dataDir, for planet pl (1-8)"""
    
    exts = ['mer','ven','ear','mar','jup','sat','ura','nep']
    fileName = dataDir+'/VSOP87D.'+exts[pl-1]
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
def computeLBR(jde, lonTerms,latTerms,radTerms):
    Tjm = dt.jd2tjm(jde)  # Time since 2000 in Julian millennia
    
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

    return lon%pi2, lat, rad


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


def plMagn(pl, distPS, distPE, distSE):
    """Compute the magnitude of planet pl (1-2, 4-9)  -  Expl.Suppl.tt.Astr.Almanac 3rd Ed, Table 10.6, p.413 + errata!"""

    #               Mer    Ven1   Ven2   Mars   Jup    Sat    Ur     Nep    Pl
    a0 = np.array([-0.60, -4.47,  0.98, -1.52, -9.40, -8.88, -7.19, -6.87, -1.01])
    a1 = np.array([ 4.98,  1.03, -1.02,   1.6,   0.5,   4.4,   0.2,     0,     0]) * 1e-2
    a2 = np.array([-4.88,  0.57,     0,     0,     0,     0,     0,     0,     0]) * 1e-4
    a3 = np.array([ 3.02,  0.13,     0,     0,     0,     0,     0,     0,     0]) * 1e-6
    
    phAng = np.degrees( np.arccos( (distPS**2 + distPE**2 - distSE**2) / (2*distPS*distPE) ) )  # Phase angle (deg!)
    
    if(pl==2 and phAng>163.6): pl = 3  # Venus 2
    pl = pl-1  # 1-9 -> 0-8
    
    mag = 5*np.log10(distPS*distPE) + a0[pl] + a1[pl]*phAng + a2[pl]*phAng**2 + a3[pl]*phAng**3
    return mag
    

def satRingMagn(JD, lon,lat):
    """Compute the magnitude of Saturn's rings from the JD and Saturns geocentric, ecliptical coordinates (in rad)"""
    tJC = dt.jd2tjc(JD)  # Time since 2000 in Julian centuries
    
    incl = 0.49                # Inclination of Saturn's rotation axis (rad)
    ascNod = 2.96 + 0.024*tJC  # Ascending node of Saturn's orbit (rad)
    
    sinB = np.sin(incl) * np.cos(lat) * np.sin(lon-ascNod)  -  np.cos(incl) * np.sin(lat)
    satRingMagn = -2.60*abs(sinB) + 1.25*(sinB)**2
    
    return satRingMagn
    #  As is: mean abs. dev. from full expression: 0.014m, max: 0.041m (10^5 trials, last 5ka)  (using phi iso DeltaU)
    

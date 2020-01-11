#  Copyright (c) 2019-2020  Marc van der Sluys - marc.vandersluys.nl
#   
#  This file is part of the HistAstro Python package,
#  see: http://astro.ru.nl/~sluys/HistAstro/
#   
#  This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#  
#  This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
#  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License along with this code.  If not, see
#  <http://www.gnu.org/licenses/>.


"""Planet functions for HistAstro."""

import math as m
import numpy.core as np
import histastro.datetime as dt

pi2 = m.pi*2


def readVSOP(dataDir, pl):
    """
    Read the periodic terms for a heliocentric ecliptical planet position from a VSOP87D.* file.
    
    Args:
      dataDir (str):  Directory where the VSOP87D.* files are located.
      pl (int):       Planet ID: 1-8 = Mercury - Neptune.
    
    Returns:
      tuple (double,double,double):  Tuple containing (lonTerms, latTerms, radTerms):
    
      - lonTerms (double):  Numpy array containing VSOP87D periodic terms for heliocentric ecliptical
                            longitude.
      - latTerms (double):  Numpy array containing VSOP87D periodic terms for heliocentric ecliptical latitude.
      - radTerms (double):  Numpy array containing VSOP87D periodic terms for heliocentric distance.
    
    References:
      - [Bretagnon & Francou, A&A 202, 309 (1988)](https://ui.adsabs.harvard.edu/abs/1988A%26A...202..309B).
      - [VSOP87D.* data files](http://cdsarc.u-strasbg.fr/viz-bin/Cat?cat=VI/81): click on FTP, download VSOP87D.*.
    
    """
    
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
        # print(var,power,nTerm)
        if line == '': break  # EoF
        
        for iLine in range(nTerm):
            line = inFile.readline()
            a,b,c = formatBody.read(line)
            # print(iLine, var,power, a,b,c)
            
            if var == 1: lonTerms.append([power, a,b,c])  # var=1: ecliptic longitude
            if var == 2: latTerms.append([power, a,b,c])  # var=2: ecliptic latitude
            if var == 3: radTerms.append([power, a,b,c])  # var=3: radial distance
            
    return lonTerms,latTerms,radTerms


def computeLBR(jde, lonTerms,latTerms,radTerms):
    """
    Compute the heliocentric ecliptical coordinates for a planet from its VSOP87D periodic terms.
    
    Args:
      jde (double):       Julian Day in dynamical time, i.e., corrected for Delta T (days).
      lonTerms (double):  Numpy array containing VSOP87D periodic terms for heliocentric ecliptical
                          longitude.
      latTerms (double):  Numpy array containing VSOP87D periodic terms for heliocentric ecliptical latitude.
      radTerms (double):  Numpy array containing VSOP87D periodic terms for heliocentric distance.
    
    Returns:
      tuple (double,double,double):  Tuple containing (lon, lat, rad):
      
        - lon (double):  Heliocentric ecliptical longitude (rad).
        - lat (double):  Heliocentric ecliptical latitude (rad).
        - rad (double):  Heliocentric distance (km).
    
    References:
      - [Bretagnon & Francou, A&A 202, 309 (1988)](https://ui.adsabs.harvard.edu/abs/1988A%26A...202..309B).
    
    """
    
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
        
    return lon % pi2, lat, rad



def hc2gc(l0,b0,r0, l,b,r):
    """
    Convert the heliocentric spherical coordinates of an object to geocentric spherical coordinates.
    
    Args:
      l0 (double):  Heliocentric ecliptical longitude of the Earth (rad).
      b0 (double):  Heliocentric ecliptical latitude of the Earth (rad).
      r0 (double):  Heliocentric ecliptical distance of the Earth (km).
    
      l (double):   Heliocentric ecliptical longitude of the other object (rad).
      b (double):   Heliocentric ecliptical latitude of the other object (rad).
      r (double):   Heliocentric ecliptical distance of the other object (km).
    
    Returns:
      tuple (double,double,double):  Tuple containing (lon, lat, rad):
      
        - lon (double):  Geocentric ecliptical longitude (rad).
        - lat (double):  Geocentric ecliptical latitude (rad).
        - rad (double):  Geocentric distance (km).
    
    Note:
      - The heliocentric spherical coordinates of the Earth and the other object are first converted to
        rectangular coordinates.  Then the difference between the two sets is taken, yielding geocentric
        rectangular coordinates.  Finally, this difference is converted back to geocentric spherical
        coordinates.
      - When the coordinates of the Earth are replaced by those of a different object, the coordinates of the
        second object are computed, as seen from the centre of the first object.
      - This function is useful when using the VSOP87D files.
    
    """
    
    x = r * m.cos(b) * m.cos(l)  -  r0 * m.cos(b0) * m.cos(l0)
    y = r * m.cos(b) * m.sin(l)  -  r0 * m.cos(b0) * m.sin(l0)
    z = r * m.sin(b)             -  r0 * m.sin(b0)
    
    # Convert geocentric rectangular to geocentric spherical coordinates:
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



def xyz_hc2lbr_gc(x0,y0,z0, x,y,z):
    """
    Convert the heliocentric rectangular coordinates of an object to geocentric spherical coordinates.
    
    Args:
      x0 (double):  Heliocentric ecliptical x-coordinate of the Earth.
      y0 (double):  Heliocentric ecliptical y-coordinate of the Earth.
      z0 (double):  Heliocentric ecliptical z-coordinate of the Earth.
    
      x (double):   Heliocentric ecliptical x-coordinate of the other object.
      y (double):   Heliocentric ecliptical y-coordinate of the other object.
      z (double):   Heliocentric ecliptical z-coordinate of the other object.
    
    Returns:
      tuple (double,double,double):  Tuple containing (lon, lat, rad):
      
        - lon (double):  Geocentric ecliptical longitude (rad).
        - lat (double):  Geocentric ecliptical latitude (rad).
        - rad (double):  Geocentric distance (km).
    
    Note:
      - The distance units of the arguments must all be the same.
      - When the coordinates of the Earth are replaced by those of a different object, the coordinates of the
        second object are computed, as seen from the centre of the first object.
      - This function is useful when using the VSOP87C files.
    
    """
    
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



def magnPlanet(pl, distPS, distPE, distSE):
    """Compute the apparent visual magnitude of a planet.
    
    Args:
      pl (int):   Planet ID (1-2, 4-9 for Mercury-Venus, Mars-Pluto).)
      distPS:     Heliocentric distance of the planet (AU).
      distPE:     Geocentric distance of the planet (AU).
      distSE:     Heliocentric distance of the Earth (AU).
    
    Returns:
      double:  Apparent visual magnitude of the planet.
    
    References:
      - [Explanatory Supplement to the Astronomical Almanac 3rd Ed, Table 10.6, p.413 +
        errata](https://aa.usno.navy.mil/publications/docs/exp_supp.php)

    """
    
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
    


def magnSatRing(JD, lon,lat):
    """Compute the apparent visual magnitude of Saturn's rings.
    
    Args:
      JD (double):   Julian Day (days).
      lon (double):  Geocentric ecliptical longitude of Saturn (rad).
      lat (double):  Geocentric ecliptical latitude of Saturn (rad).
    
    Returns:
      double:  magnSatRing: Apparent visual magnitude of the planet.
    
    References:
      - [Explanatory Supplement to the Astronomical Almanac 3rd Ed, Table 10.6, p.413 +
        errata](https://aa.usno.navy.mil/publications/docs/exp_supp.php)
    
    Note:
      - This is a simplified expression, which uses phi instead of Delta U.  The mean absolute deviation from
        the full expression: is 0.014m, the maximum deviation found: 0.041m in 10^5 trials, over the last 5000
        years.

    """
    
    tJC = dt.jd2tjc(JD)  # Time since 2000 in Julian centuries
    
    incl = 0.49                # Inclination of Saturn's rotation axis (rad)
    ascNod = 2.96 + 0.024*tJC  # Ascending node of Saturn's orbit (rad)
    
    sinB = np.sin(incl) * np.cos(lat) * np.sin(lon-ascNod)  -  np.cos(incl) * np.sin(lat)
    magnSatRing = -2.60*abs(sinB) + 1.25*(sinB)**2
    
    return magnSatRing
    

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


"""Coordinate transformations and related functions for HistAstro."""


# Modules:
import math as m
import numpy.core as np
from histastro.constants import pi2, r2d,as2r, earthRad,AU
import histastro.datetime as dt



def obliquity(jd):
    """
    Compute the obliquity of the ecliptic in radians from the JD(E).
    
    Arguments:
      jd (double):  Julian day (days).
    
    Returns:
      double:  eps: Obliquity of the ecliptic (rad).
    
    References:
      - Seidelman 1992, Eq. 3.222-1.
    
    """
    
    tJC = dt.jd2tjc(jd)  # Time in Julian centuries since J2000.0
    eps = 0.409092804 - 2.269655e-4*tJC - 2.86e-9*tJC**2 + 8.78967e-9*tJC**3  # Obliquity of the ecliptic (rad)
    
    return eps


def eq2ecl(ra,dec, eps):
    """
    Convert equatorial coordinates to ecliptical.

    Arguments:
      ra (double):   Right ascension (rad).
      dec (double):  Declination (rad).
      eps (double):  Obliquity of the ecliptic (rad).
    
    Returns:
      tuple (double,double): tuple containing (lon, lat):
    
        - lon (double):  Ecliptical longitude (rad).
        - lat (double):  Ecliptical latitude (rad).
    
    """
    
    lon = np.arctan2( np.sin(ra)  * m.cos(eps) + np.tan(dec) * m.sin(eps),  np.cos(ra) ) % pi2
    lat =  np.arcsin( np.sin(dec) * m.cos(eps) - np.cos(dec) * m.sin(eps) * np.sin(ra) )
    
    return lon,lat


def ecl2eq(lon,lat, eps):
    """Convert (geocentric) spherical ecliptical coordinates to spherical equatorial coordinates.
    
    Arguments:
      lon (double):  Ecliptical longitude (rad).
      lat (double):  Ecliptical latitude (rad).
      eps (double):  Obliquity of the ecliptic (rad).
    
    Returns:
      tuple (double,double): tuple containing (ra, dec):
    
        - ra (double):   Right ascension (rad).
        - dec (double):  Declination (rad).
    
    References:
      - [Explanatory Supplement to the Astronomical Almanac 3rd Ed,
        Eq.14.43](https://aa.usno.navy.mil/publications/docs/exp_supp.php)

    """
    
    ra  = np.arctan2( np.sin(lon) * np.cos(eps)  -  np.tan(lat) * np.sin(eps),  np.cos(lon) ) % pi2
    dec =  np.arcsin( np.sin(lat) * np.cos(eps)  +  np.cos(lat) * np.sin(eps) * np.sin(lon) )
    
    return ra,dec


def par2horiz(ha,dec, phi):
    """Convert parallactic coordinates to horizontal.

    Arguments:
      ha (double):   Hour angle (rad).
      dec (double):  Declination (rad).
      phi (double):  Geographical latitude (rad, N>0).
    
    Returns:
      tuple (double,double): tuple containing (az, alt):
    
        - az (double):   Azimuth (rad, S=0).
        - alt (double):  Altitude (rad).
    
    """
    
    az  = np.arctan2( np.sin(ha),   np.cos(ha) * m.sin(phi) - np.tan(dec) * m.cos(phi) ) % pi2
    alt = np.arcsin(  np.sin(dec) * m.sin(phi) + np.cos(ha) * np.cos(dec) * m.cos(phi) )
    
    return az,alt


def properMotion(startJD,targetJD, ra,dec, pma,pmd):
    """Compute the proper motion from startJD to targetJD for the given positions and proper motions.
    
    Arguments:
      startJD (double):   Julian day of the initial epoch (days).
      targetJD (double):  Julian day of the target epoch (days).
    
      ra (double):        Right ascension (numpy array, rad).
      dec (double):       Declination (numpy array, rad).
    
      pma (double):       Proper motion in right ascension (numpy array, rad/yr).
      pmd (double):       Proper motion in declination (numpy array, rad/yr).

    Returns:
      tuple (double,double): tuple containing (raTarget, decTarget):
    
        - raTarget (double):   Right ascension for the target epoch (rad).
        - decTarget (double):  Declination for the target epoch (rad).

    """
    
    dtYr   = (targetJD - startJD)/365.25
    raTarget  = ra  + pma*dtYr / np.cos(dec)
    decTarget = dec + pmd*dtYr
    
    return raTarget,decTarget


def precessHip(jd, ra,dec):
    """Compute precession in equatorial coordinates from the Hipparcos equinox (J2000) to that of the specified JD.
    
    Arguments:
      jd (double):   Julian day (days).
      ra (double):   Right ascension (rad).
      dec (double):  Declination (rad).
    
    Returns:
      tuple (double,double): tuple containing (raTarget, decTarget):
    
        - raNew (double):   Right ascension for the target equinox (rad).
        - decNew (double):  Declination for the target equinox (rad).
    
    """
    
    tJC  = dt.jd2tjc(jd)  # Time in Julian centuries since J2000.0
    tJC2 = tJC**2
    tJC3 = tJC*tJC2
    
    zeta  = (2306.2181*tJC + 0.30188*tJC2 + 0.017998*tJC3)*as2r
    z     = (2306.2181*tJC + 1.09468*tJC2 + 0.018203*tJC3)*as2r
    theta = (2004.3109*tJC - 0.42665*tJC2 - 0.041833*tJC3)*as2r
    
    raNew  = (np.arctan2( np.sin(ra + zeta) * np.cos(dec),  np.cos(ra + zeta) * m.cos(theta) * np.cos(dec) - m.sin(theta) * np.sin(dec) ) + z) % pi2
    decNew = np.arcsin( np.cos(ra + zeta) * m.sin(theta) * np.cos(dec)  +  m.cos(theta) * np.sin(dec) )
    
    return raNew,decNew


def geoc2topoc_ecl(gcLon,gcLat, gcDist,gcRad, eps,lst, obsLat,obsEle=0, debug=False):
    """Convert spherical ecliptical coordinates from the geocentric to the topocentric system.
    
    Arguments:
      gcLon (double):   Geocentric ecliptic longitude (rad).
      gcLat (double):   Geocentric ecliptic latitude (rad).
      gcDist (double):  Geocentric distance (AU).
      gcRad (double):   Geocentric semi-diameter (rad).
      
      eps (double):     Obliquity of the ecliptic (rad).
      lst (double):     Local sidereal time (rad).
      
      obsLat (double):  Geographical latitude of the observer (rad).
      obsEle (double):  Altitude/elevation of the observer above sea level (metres, optional, default value = 0).
      
      debug (double):   Print debug output (True/False, optional, default value = True).
    
    Returns:
      tuple (double,double,double): tuple containing (tcLon, tcLat, tcRad):
    
        - tcLon (double):  Topocentric ecliptic longitude (rad).
        - tcLat (double):  Topocentric ecliptic latitude (rad).
        - tcRad (double):  Topocentric semi-diameter (rad).
    
    """
    
    # Meeus, Ch.11, p.82:
    Req = earthRad*1000      # Equatorial radius of the Earth in metres (same units as the elevation)
    #                        (http://earth-info.nga.mil/GandG/publications/tr8350.2/wgs84fin.pdf)
    RpolEq = 0.996647189335  # Rpol/Req = 1-f: flattening of the Earth - WGS84 ellipsoid 
    
    u  = m.atan(RpolEq*m.tan(obsLat))
    RsinPhi = RpolEq*m.sin(u) + obsEle/Req * m.sin(obsLat)
    RcosPhi = m.cos(u)        + obsEle/Req * m.cos(obsLat)
    
    sinHp = m.sin(earthRad/AU)/(gcDist/AU)  # Sine of the horizontal parallax, Meeus, Eq. 40.1
    
    # Meeus, Ch.40, p.282:
    N  = m.cos(gcLon)*m.cos(gcLat) - RcosPhi*sinHp*m.cos(lst)
    
    tcLon = m.atan2( m.sin(gcLon)*m.cos(gcLat) - sinHp*(RsinPhi*m.sin(eps) + RcosPhi*m.cos(eps)*m.sin(lst)), N ) % pi2  # Topocentric longitude
    tcLat = m.atan((m.cos(tcLon)*(m.sin(gcLat) - sinHp*(RsinPhi*m.cos(eps) - RcosPhi*m.sin(eps)*m.sin(lst))))/N)         # Topocentric latitude
    tcRad = m.asin(m.cos(tcLon)*m.cos(tcLat)*m.sin(gcRad)/N)                                                   # Topocentric semi-diameter
    
    # print(gcDist, gcDist*gcRad/tcRad)
    
    if debug:
        print()
        print('geoc2topoc_ecl():')
        print('%10s  %25s  %25s' % ('', 'rad/km/...','deg'))
        print()
        print('%10s  %25.15f' % ('Req: ', Req) )
        print('%10s  %25.15f' % ('RpolEq: ', RpolEq) )
        print()
        print('%10s  %25.15f  %25.15f' % ('u: ', u, u*r2d) )
        print('%10s  %25.15f' % ('RsinPhi: ', RsinPhi) )
        print('%10s  %25.15f' % ('RcosPhi: ', RcosPhi) )
        print()
        print('%10s  %25.15f' % ('sinHp: ', sinHp) )
        print('%10s  %25.15f  %25.15f' % ('N: ', N, N*r2d) )
        print()
        print('%10s  %25.15f  %25.15f' % ('tcLon: ', tcLon, tcLon*r2d) )
        print('%10s  %25.15f  %25.15f' % ('tcLat: ', tcLat, tcLat*r2d) )
        print('%10s  %25.15f  %25.15f' % ('tcRad: ', tcRad, tcRad*r2d) )
        
    return tcLon,tcLat,tcRad




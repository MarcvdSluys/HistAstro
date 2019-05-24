"""HistAstro coordinate functions"""

# Modules:
import math as m
import numpy.core as np
from histastro.constants import pi,pi2, r2d,as2r, earthRad,AU
import histastro.datetime as dt



def obliquity(jd):
    """Compute the obliquity of the ecliptic in radians from the JD(E)  (Seidelman 1992, Eq. 3.222-1)"""
    
    tJC = dt.jd2tjc(jd)  # Time in Julian centuries since J2000.0
    eps = 0.409092804 - 2.269655e-4*tJC - 2.86e-9*tJC**2 + 8.78967e-9*tJC**3  # Obliquity of the ecliptic (rad)
    return eps


def eq2ecl(ra,dec, eps):
    """Convert equatorial coordinates to ecliptical"""
    
    lon = np.arctan2( np.sin(ra)  * m.cos(eps) + np.tan(dec) * m.sin(eps),  np.cos(ra) ) % pi2
    lat =  np.arcsin( np.sin(dec) * m.cos(eps) - np.cos(dec) * m.sin(eps) * np.sin(ra) )
    return lon,lat


def ecl2eq(lon,lat, eps):
    """Convert (geocentric) spherical ecliptical coordinates lon, lat (and eps) to spherical equatorial coordinates RA, Dec.
    See Expl. Supl. t.t. Astronimical Almanac 3rd Ed, Eq.14.43"""
    
    ra  = np.arctan2( np.sin(lon) * np.cos(eps)  -  np.tan(lat) * np.sin(eps),  np.cos(lon) ) % pi2
    dec =  np.arcsin( np.sin(lat) * np.cos(eps)  +  np.cos(lat) * np.sin(eps) * np.sin(lon) )
    return ra,dec


def par2horiz(ha,dec, phi):
    """Convert parallactic coordinates to horizontal"""
    
    az  = np.arctan2( np.sin(ha),   np.cos(ha) * m.sin(phi) - np.tan(dec) * m.cos(phi) ) % pi2
    alt = np.arcsin(  np.sin(dec) * m.sin(phi) + np.cos(ha) * np.cos(dec) * m.cos(phi) )
    return az,alt


def properMotion(startJD,targetJD, ra,dec, pma,pmd):
    """Compute the proper motion from startJD to targetJD for the positions given in (numpy arrays) ra and dec
    (in rad) and proper motions in pma,pmd (in rad/yr)"""
    
    dtYr   = (targetJD - startJD)/365.25
    raTarget  = ra  + pma*dtYr / np.cos(dec)
    decTarget = dec + pmd*dtYr
    return raTarget,decTarget


def precessHip(jd, ra,dec):
    """Compute precession in equatorial coordinates from the Hipparcos epoch (J2000) to the specified JD"""
    
    tJC  = dt.jd2tjc(jd)  # Time in Julian centuries since J2000.0
    tJC2 = tJC**2
    tJC3 = tJC*tJC2
    
    zeta  = (2306.2181*tJC + 0.30188*tJC2 + 0.017998*tJC3)*as2r
    z     = (2306.2181*tJC + 1.09468*tJC2 + 0.018203*tJC3)*as2r
    theta = (2004.3109*tJC - 0.42665*tJC2 - 0.041833*tJC3)*as2r
    
    raNew  = (np.arctan2( np.sin(ra + zeta) * np.cos(dec),  np.cos(ra + zeta) * m.cos(theta) * np.cos(dec) - m.sin(theta) * np.sin(dec) ) + z) % pi2
    decNew = np.arcsin( np.cos(ra + zeta) * m.sin(theta) * np.cos(dec)  +  m.cos(theta) * np.sin(dec) )
    return raNew,decNew


def geoc2topoc_ecl(gcLon,gcLat, gcDist,gcRad, eps,lst, obsLat,obsEle=0):
    """ Convert spherical ecliptical coordinates from the geocentric to the topocentric system
    
    Input parameters:
    - gcLon   Geocentric ecliptic longitude (rad)
    - gcLat   Geocentric ecliptic latitude (rad)
    - gcDist  Geocentric distance
    - gcRad   Geocentric semi-diameter (rad)
    
    - eps   Obliquity of the ecliptic (rad)
    - lst   Local sidereal time (rad)
    
    - obsLat   Latitude of the observer (rad)
    - obsEle   Altitude/elevation of the observer above sea level (metres, optional)
    
    Return values:
    - tcLon   Topocentric ecliptic longitude (rad)
    - tcLat   Topocentric ecliptic latitude (rad)
    - tcRad   Topocentric semi-diameter (rad)
    
    See:
    - Meeus, Astronomical Algorithms, 1998, Ch. 11 and 40
    """
    
    # Meeus, Ch.11, p.82:
    ba = 0.996647189335    # b/a = 1-f: flattening of the Earth - WGS84 ellipsoid 
    re = earthRad*1000     # Equatorial radius of the Earth in metres (same units as the elevation)
    #                      (http://earth-info.nga.mil/GandG/publications/tr8350.2/wgs84fin.pdf)
    
    u  = m.atan(ba*m.tan(obsLat))
    rs = ba*m.sin(u) + obsEle/re * m.sin(obsLat)
    rc = m.cos(u)    + obsEle/re * m.cos(obsLat)
    
    sinHp = m.sin(earthRad/AU)/(gcDist/AU)  # Sine of the horizontal parallax, Meeus, Eq. 40.1
    
    # Meeus, Ch.40, p.282:
    N  = m.cos(gcLon)*m.cos(gcLat) - rc*sinHp*m.cos(lst)
    
    tcLon = m.atan2( m.sin(gcLon)*m.cos(gcLat) - sinHp*(rs*m.sin(eps) + rc*m.cos(eps)*m.sin(lst)) , N ) % pi2  # Topocentric longitude
    tcLat = m.atan((m.cos(tcLon)*(m.sin(gcLat) - sinHp*(rs*m.cos(eps) - rc*m.sin(eps)*m.sin(lst))))/N)         # Topocentric latitude
    tcRad = m.asin(m.cos(tcLon)*m.cos(tcLat)*m.sin(gcRad)/N)                                                   # Topocentric semi-diameter
    
    # print(gcDist, gcDist*gcRad/tcRad)
    
    return tcLon,tcLat,tcRad




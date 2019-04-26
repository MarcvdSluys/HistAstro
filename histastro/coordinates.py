"""HistAstro coordinate functions"""

# Modules:
import math as m
import numpy as np

# Constants:
pi    = m.pi
pi2   = 2*pi
r2d   = m.degrees(1)  # Radians to degrees
d2r   = 1.0/r2d       # Degrees to radians
#r2h  = r2d/15
h2r   = d2r*15        # Hours to radians
as2r  = d2r/3.6e3     # Arcseconds to radians
mas2r = as2r/1000.0   # Milliarcseconds to radians


def obliquity(JD):
    """Compute the obliquity of the ecliptic in radians from the JD(E)  (Seidelman 1992, Eq. 3.222-1)"""
    tJC = (JD - 2451545.0)/36525
    #eps = 23.4392911*d2r
    #eps += (-46.815*tJC - 0.00059*tJC**2 + 0.001813*tJC**3)*as2r
    eps = 0.409092804 - 2.269655e-4*tJC - 2.86e-9*tJC**2 + 8.78967e-9*tJC**3  # Checked
    return eps

def eq2ecl(ra,dec, eps):
    """Convert equatorial coordinates to ecliptical"""
    lon = np.arctan2( np.sin(ra)  * m.cos(eps) + np.tan(dec) * m.sin(eps),  np.cos(ra) ) % pi2
    lat =  np.arcsin( np.sin(dec) * m.cos(eps) - np.cos(dec) * m.sin(eps) * np.sin(ra) )
    return lon,lat


def par2horiz(ha,dec, phi):
    """Convert parallactic coordinates to horizontal"""
    az  = np.arctan2( np.sin(ha),   np.cos(ha) * m.sin(phi) - np.tan(dec) * m.cos(phi) ) % pi2
    alt = np.arcsin(  np.sin(dec) * m.sin(phi) + np.cos(ha) * np.cos(dec) * m.cos(phi) )
    return az,alt

def properMotion(startJD,targetJD, ra,dec, pma,pmd):
    """Compute the proper motion from startJD to targetJD for the positions given in (numpy arrays) ra and dec
    (in rad) and proper motions in pma,pmd (in rad/yr)

    """
    dtYr = (targetJD - startJD)/365.25
    raOld  = ra  + pma*dtYr / np.cos(dec)
    decOld = dec + pmd*dtYr
    return raOld,decOld


def precessHip(jd, ra,dec):
    """Compute precession in equatorial coordinates from the Hipparcos epoch (J2000) to the specified JD"""
    tJC = (jd - 2451545.0)/36525  # Time in Julian centuries since J2000.0
    tJC2 = tJC**2
    tJC3 = tJC*tJC2
    
    zeta  = (2306.2181*tJC + 0.30188*tJC2 + 0.017998*tJC3)*as2r
    z     = (2306.2181*tJC + 1.09468*tJC2 + 0.018203*tJC3)*as2r
    theta = (2004.3109*tJC - 0.42665*tJC2 - 0.041833*tJC3)*as2r
    
    raNew  = np.arctan2( np.sin(ra + zeta) * np.cos(dec),  np.cos(ra + zeta) * m.cos(theta) * np.cos(dec) - m.sin(theta) * np.sin(dec) ) + z
    decNew = np.arcsin( np.cos(ra + zeta) * m.sin(theta) * np.cos(dec)  +  m.cos(theta) * np.sin(dec) )
    return raNew,decNew





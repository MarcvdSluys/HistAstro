"""HistAstro date and time functions"""


# Modules:
import math as m
from histastro.constants import pi2, jd1900,jd2000


def julianDay(year,month,day):
    """Compute the Julian Day for a given year, month and (decimal) day UT"""
    year0 = year
    if month <= 2:  # Jan and Feb are month 13 and 14 of the previous year
       year -= 1
       month += 12
       
    b = 0; a=0
    if year0 > 1582:     # Assume a Gregorian date
       a = m.floor(year/100.0)
       b = 2 - a + m.floor(a/4.0)
    
    jd = m.floor(365.25*(year+4716)) + m.floor(30.6001*(month+1)) + day + b - 1524.5
    return jd


def jd2tjc(jd):
    """Return the time since 2000 expressed in Julian centuries"""
    
    return (jd - 2451545.0)/36525


def jd2tjm(jd):
    """Return the time since 2000 expressed in Julian millennia"""
    
    return (jd - 2451545.0)/365250


def gmst(jd):
    """Calculate Greenwich Mean Sidereal Time for any instant, in radians.
    Explanatory Supplement to the Astronomical Almanac, 3rd ed, Eq. 6.66 (2012)"""
    
    tjd  = jd - jd2000                      # Julian Days after 2000.0 UT
    tjd2 = tjd**2
    tjd4 = tjd2**2
           
    gmst = 4.89496121088131 + 6.30038809894828323*tjd + 5.05711849e-15*tjd2 - 4.378e-28*tjd2*tjd - 8.1601415e-29*tjd4 \
        - 2.7445e-36*tjd4*tjd  # Eq. 6.66, removed Delta-T term, hence replaced the first term
    
    return gmst % pi2


def DeltaT(jd):
    """Return a rough estimation for the value of DeltaT (s), using a lenghtening of the day of 1.8 ms/century and
    that the minimum of the parabola is DeltaT=0 in 1900."""
    #return 0.5 * 1.8e-3/86400/(36525*86400) * ((jd-jd1900)*86400)**2
    return 0.5 * 1.8e-3 / 36525 * (jd-jd1900)**2


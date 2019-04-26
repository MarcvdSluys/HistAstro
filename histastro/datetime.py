"""HistAstro date and time functions"""


# Modules:
import math as m


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


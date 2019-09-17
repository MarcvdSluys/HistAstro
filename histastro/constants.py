#!/bin/env python3

import math as m

pi   = m.pi           # pi
pi2  = 2*pi           # 2 pi
pio2 = pi/2           # pi/2

r2d  = m.degrees(1)   # Radians to degrees
d2r  = m.radians(1)   # Degrees to radians

h2r   = d2r*15        # Hours to radians
as2r  = d2r/3.6e3     # Arcseconds to radians
mas2r = as2r/1000.0   # Milliarcseconds to radians

jd1820 = 2385801      # JD in 1820  (for Delta-T fit)
jd1900 = 2415021      # JD in 1900
jd2000 = 2451545      # JD in 2000.0

earthRad = 6378.1366  # Earth radius in km
moonRad = 1737.5      # Moon radius in km
AU = 1.495978707e8    # Astronomical unit in km


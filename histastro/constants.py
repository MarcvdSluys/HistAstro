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


"""Definition of constants for HistAstro."""

# Modules:
import math as m

pi   = m.pi;           """pi"""
pi2  = 2*pi;           """2 pi"""
pio2 = pi/2;           """pi/2"""

r2d  = m.degrees(1);   """Radians to degrees"""
d2r  = m.radians(1);   """Degrees to radians"""

h2r   = d2r*15;        """Hours to radians"""
as2r  = d2r/3.6e3;     """Arcseconds to radians"""
mas2r = as2r/1000.0;   """Milliarcseconds to radians"""

jd1820 = 2385801;      """JD in 1820  (for Delta-T fit)"""
jd1900 = 2415021;      """JD in 1900"""
jd2000 = 2451545;      """JD in 2000.0"""

earthRad = 6378.1366;  """Earth radius in km"""
moonRad = 1737.5;      """Moon radius in km"""

AU = 1.495978707e8;    """Astronomical unit in km"""


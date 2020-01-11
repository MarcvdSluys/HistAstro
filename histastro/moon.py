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


"""Moon functions for HistAstro."""


import math as m
import numpy as np
from histastro.constants import pi2,r2d,jd2000,moonRad  # ,d2r


def readData(inFile='data/moonposMeeus.csv'):
    """
    Return the periodic terms for the ELP82B theory from moonposMeeus.csv.

      - Two arrays aand return them in two arrays: one for longitude and distance, and one for latitude.

    Args: 
      inFile (str): Name of the input file, including (relative or absolute) path.  Optional, default
      value = 'data/moonposMeeus.csv'.

    Returns:
      tuple (double,double):  Tuple containing (lrTerms, bTerms):

        - lrTerms (double): Numpy array containing six columns with periodic terms for longitude and distance.
          The first four columns contain arguments for both variables, the last two the coefficients for
          longitude and distance:
    
          - 1 (int): Multiplication factor for the mean elongation of the Moon.
          - 2 (int): Multiplication factor for the mean anomaly of the Sun.
          - 3 (int): Multiplication factor for the mean anomaly of the Moon.
          - 4 (int): Multiplication factor for the Moon's argument of latitude.
          - 5 Coefficient of the sine of the argument, for the longitude (rad).
          - 6 Coefficient of the cosine of the argument, for the distance (km).
    
        - bTerms (double): Numpy array containing five columns periodic terms for latitude.  The first four
          columns contain arguments, the last contains the coefficients for latitude:
    
          - 1 (int): Multiplication factor for the mean elongation of the Moon.
          - 2 (int): Multiplication factor for the mean anomaly of the Sun.
          - 3 (int): Multiplication factor for the mean anomaly of the Moon.
          - 4 (int): Multiplication factor for the Moon's argument of latitude.
          - 5 Coefficient of the sine of the argument, for the latitude (rad).
    
    References:
      - [Jean Meeus: Astronomical Algorithms, 2nd Ed. (1998)](https://www.willbell.com/math/MC1.HTM), Ch.47.

    """
    
    lrTerms = np.genfromtxt(inFile, delimiter=',', skip_header=1,  max_rows=60)  # Longitude and radius (6 columns: 4 args, 2 coefs)
    bTerms  = np.genfromtxt(inFile, delimiter=',', skip_header=61, max_rows=60)  # Latitude (5 columns: 4 args, 1 coef)
    
    return lrTerms,bTerms


def compute_lbr(jde, lrTerms,bTerms, debug=False):
    """
    Compute the geocentric ecliptic coordinates of the Moon for the equinox of date for the given JDE.
    
    Args:
      jde (double):      Julian Day in dynamical time, i.e., corrected for Delta T (days).
      lrTerms (double):  Numpy array with ELP82 periodic terms for longitude and distance.
      bTerms (double):   Numpy array with ELP82 periodic terms for latitude.
      debug (bool):      Produce debug output.  Optional, default value = False.
    
    Returns:
      tuple (double,double,double,double):  Tuple containing (lon, lat, dist, diam):
    
        - lon (double):   Geocentric ecliptic longitude of the Moon (rad).
        - lat (double):   Geocentric ecliptic longitude of the Moon (rad).
        - dist (double):  Geocentric distance of the Moon (km).
        - diam (double):  Geocentric apparent diameter of the Moon (rad).
    
    Notes:
      - A reduced version of the ELP82 theory is used.
      - The necessary ELP82 terms can be read from file using the function readData().
    
    References:
      - [Chapront-Touze & Chapront, A&A 190, 342 (1988)](https://ui.adsabs.harvard.edu/abs/1988A%26A...190..342C).
    
    """
    
    tjc   = (jde-jd2000)/36525  # Julian Centuries after 2000.0 in dynamical time
    tjc2  = tjc**2
    tjc3  = tjc*tjc2
    tjc4  = tjc2**2
    
    # Moon's mean longitude, Meeus p.338:
    lm = (3.8103408236 + 8399.7091116339958*tjc - 2.755176757e-5*tjc2 + 3.239043e-8*tjc3 - 2.6771e-10*tjc4) % pi2
    
    # Delauney arguments [d, ms, mm, f] (Meeus p.338) = [D, l', l, F] in Simon et al. 1994, Sect. 3.5:
    d  = (5.1984665298  + 7771.377144834*tjc - 3.2845e-5*tjc2  + 3.197347e-8*tjc3    - 1.5436512e-10*tjc4) % pi2  # Moon's mean elongation
    ms = (6.240060127   + 628.301955167*tjc  - 2.681e-6*tjc2   + 7.1267017e-10*tjc3)                       % pi2  # Sun's mean anomaly
    mm = (2.355555637   + 8328.691424759*tjc + 1.52566e-4*tjc2 + 2.5041e-7*tjc3      - 1.18633e-9*tjc4)    % pi2  # Moon's mean anomaly
    f  = (1.627905158   + 8433.466158061*tjc - 6.3773e-5*tjc2  - 4.94988e-9*tjc3     + 2.02167e-11*tjc4)   % pi2  # Moon's argument of latitute
    DelArgs = [d,ms,mm,f]  # Delauney arguments
    
    e  = 1 - 0.002516*tjc - 0.0000074*tjc2
    
    # Compute ELP82 terms:
    lon=0; lat=0; dist=0
    for iLine in range(60):  # Data lines
        argl = 0
        argb = 0
        for iArg in range(4):  # Integer arguments
            argl += lrTerms[iLine,iArg] * DelArgs[iArg]
            argb +=  bTerms[iLine,iArg] * DelArgs[iArg]
        
        lon  += m.sin(argl) * lrTerms[iLine,4] * e**(abs(lrTerms[iLine,1]))
        lat  += m.sin(argb) *  bTerms[iLine,4] * e**(abs(bTerms[iLine,1]))
        dist += m.cos(argl) * lrTerms[iLine,5] * e**(abs(lrTerms[iLine,1]))
        
    # Save to print intermediate results:
    lon1  = lon
    lat1  = lat
    dist1 = dist
    
    
    # Perturbations by other planets, and flattening of the Earth:
    # Meeus, p.338:
    a1 = (2.090032  +     2.301199 * tjc) % pi2  # Influence of Venus
    a2 = (0.926595  + 8364.7398477 * tjc) % pi2  # Influence of Jupiter
    a3 = (5.4707345 + 8399.6847253 * tjc) % pi2
    
    # Meeus, p.342:
    # dlon = (3958*m.sin(a1) + 1962*m.sin(lm-f) + 318*m.sin(a2)) * d2r * 1e-6
    # dlat = (-2235*m.sin(lm) + 382*m.sin(a3) + 175*m.sin(a1-f) + 175*m.sin(a1+f) + 127*m.sin(lm-mm) - 115*m.sin(lm+mm)) * d2r * 1e-6
    dlon =  6.908e-5*m.sin(a1) + 3.4243e-5*m.sin(lm-f) + 5.55e-6*m.sin(a2)
    dlat = -3.9008e-5*m.sin(lm) + 6.667e-6*m.sin(a3) + 3.0543e-6*m.sin(a1-f) + 3.0543e-6*m.sin(a1+f) + 2.2166e-6*m.sin(lm-mm) - 2.007e-6*m.sin(lm+mm)
    
    
    # Compute nutation in longitude:
    # omg  = 2.18243858558 - 33.7570459367*tjc + 3.6142278e-5*tjc2 + 3.87850944888e-8*tjc3   # Moon's mean lon. of asc.node, Meeus p.144
    # ls   = 4.89506386655 + 62.84528862*tjc                                                 # Mean long. Sun, Meeus p.144
    # dpsi = -8.338795e-5*m.sin(omg) - 6.39954e-6*m.sin(2*ls) - 1.115e-6*m.sin(2*lm) + 1.018e-6*m.sin(2*omg)
    dpsi = 0
    
    # Add mean values:
    lon  = (lon + dlon + lm + dpsi) % pi2
    lat  = lat + dlat
    dist = dist + 385000.56  # in km
    
    diam = 2*m.atan(moonRad/dist)
    
    if debug:
        print()
        print('compute_lbr():')
        print('%10s  %25s  %25s' % ('', 'rad/km/...','deg'))
        print()
        print('%10s  %25.15f' % ('jde:  ', jde) )
        print('%10s  %25.15f' % ('tjc:  ', tjc) )
        print()
        print('%10s  %25.15f  %25.15f' % ('lm:  ', lm, lm*r2d) )
        print()
        print('%10s  %25.15f  %25.15f' % ('d:   ', d, d*r2d) )
        print('%10s  %25.15f  %25.15f' % ('ms:  ', ms, ms*r2d) )
        print('%10s  %25.15f  %25.15f' % ('mm:  ', mm, mm*r2d) )
        print('%10s  %25.15f  %25.15f' % ('f:   ', f, f*r2d) )
        print()
        print('%10s  %25.15f' % ('e:   ', e) )
        print()
        print('%10s  %25.15f  %25.15f' % ('a1:  ', a1, a1*r2d) )
        print('%10s  %25.15f  %25.15f' % ('a2:  ', a2, a2*r2d) )
        print('%10s  %25.15f  %25.15f' % ('a3:  ', a3, a3*r2d) )
        print()
        print("ELP82 sums:")
        print('%10s  %25.15f  %25.15f' % ('lon: ', lon1, lon1*r2d) )
        print('%10s  %25.15f  %25.15f' % ('lat: ', lat1, lat1*r2d) )
        print('%10s  %25.15f' % ('dist:', dist1) )
        print()
        print('%10s  %25.15f  %25.15f' % ('dlon: ', dlon, dlon*r2d) )
        print('%10s  %25.15f  %25.15f' % ('dlat: ', dlat, dlat*r2d) )
        print()
        print('%10s  %25.15f  %25.15f' % ('lon:  ', lon, lon*r2d) )
        print('%10s  %25.15f  %25.15f' % ('lat:  ', lat, lat*r2d) )
        print('%10s  %25.15f' % ('dist:  ', dist) )
        print('%10s  %25.15f  %25.15f' % ('diam:  ', diam, diam*r2d) )
        print()
        #    print('%10s  %25.15f  %25.15f' % ('x:  ', x*r2d) )
    
    return lon,lat,dist,diam



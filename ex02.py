#!/bin/env python3

# Modules:
import math as m
import numpy as np
import matplotlib.pyplot as plt

import histastro.coordinates as coord
import histastro.datetime as dt



# Constants:
pi    = m.pi
pi2   = 2*pi
r2d   = m.degrees(1)  # Radians to degrees
d2r   = 1.0/r2d       # Degrees to radians
h2r   = d2r*15        # Hours to radians
as2r  = d2r/3.6e3     # Arcseconds to radians
mas2r = as2r/1000.0   # Milliarcseconds to radians


    
# Read the Hipparcos catalogue
# Columns: 0: hip#, 1: vmag, 2: ra (rad), 3: dec (rad), 4: pmRA (mas/yr), 5: pmDec (mas/yr), 6: ErRA (?), 7:
# ErDec (?), 8: ErPa (mas/yr), 9: ErPd (mas/yr)

hip = np.loadtxt('data/BrightStars.csv', skiprows=1, delimiter=',', usecols=(1,2,3))    # Read the numbers in columns 1-3: V, ra, dec - amazingly, this can also read BrightStars.csv.gz!

mag = hip[:,0]
ra  = hip[:,1]                           # Right ascension (rad)
dec = hip[:,2]                           # Declination (rad)

# Set magnitude and map limits:
Mlim = 7.0  # Magnitude limit
sizes = 30*(0.5 + (Mlim-mag)/3.0)**2     # Scale inversely with magnitude.  Square, since scatter() uses surface area


# Plot horizontal map:
phi = 51.178*d2r  # Nijmegen

ha = -6*h2r - ra  # At the spring equinox and sunrise ra_sun = 0, ha_sun=-6h
az,alt = coord.par2horiz(ha,dec, phi)


plt.figure(figsize=(10,5.5))                   # Set png size to 1000x550 (dpi=100)
#plt.figure(figsize=(20,11))                   # Set png size to 2000x1100 (dpi=100)

azMin  = 225*d2r
azMax  = 305*d2r
altMin = 0*d2r
altMax = 40*d2r

Mlim = 4.5  # Magnitude limit
sizes = 20*(0.5 + (Mlim-mag)/3.0)**2     # Scale inversely with magnitude.  Square, since scatter() uses surface area

sel = np.logical_and(az > azMin, az < azMax)
sel = np.logical_and(sel, alt > altMin)
sel = np.logical_and(sel, alt < altMax)
sel = np.logical_and(sel, mag < Mlim)
#sel = az < 1e6  # Select all stars for plotting

# Create a scatter plot:
plt.scatter(az[sel]*r2d, alt[sel]*r2d, s=sizes[sel])

plt.axis('scaled')                                      # Set axes to a 'square grid' by moving the plot box inside the figure
plt.axis([azMin*r2d,azMax*r2d, altMin*r2d,altMax*r2d])  # Select Aries (RA=26-50 deg, dec=10-30 deg)
plt.xlabel(r'A ($^\circ$)')                             # Label the horizontal axis
plt.ylabel(r'h ($^\circ$)')                             # Label the vertical axis - use LaTeX for symbols

plt.tight_layout()                    # Use narrow margins
plt.savefig("ex02_horizontal.png")    # Save the plot as png
#plt.show()                           # Show the plot to screen
plt.close()                           # Close the plot



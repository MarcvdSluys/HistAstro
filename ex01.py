#!/bin/env python3

# Modules:
import math as m
import numpy as np
import matplotlib.pyplot as plt

import histastro.coordinates as coord



# Constants:
r2d   = m.degrees(1)  # Radians to degrees
d2r   = 1.0/r2d       # Degrees to radians


    


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

raMin  = 26.0*d2r
raMax  = 50.0*d2r
decMin = 10.0*d2r
decMax = 30.0*d2r


# Convert equatorial to ecliptic coordinates:
eps = 0.40909280  # For 2000 in rad
lon,lat = coord.eq2ecl(ra,dec, eps)


# Use five coordinates for the corners of the map, in order to plot four lines below
lonMinMax,latMinMax = coord.eq2ecl([raMin,raMax,raMax,raMin,raMin],[decMin,decMin,decMax,decMax,decMin], eps)
lonMin = min(lonMinMax)
lonMax = max(lonMinMax)
latMin = min(latMinMax)
latMax = max(latMinMax)

# Select stars to plot:
sel = np.logical_and(lon > lonMin, lon < lonMax)
sel = np.logical_and(sel, lat > latMin)
sel = np.logical_and(sel, lat < latMax)
#sel = lon < 1e6  # Select all stars for plotting



# Plot ecliptic map:
plt.figure(figsize=(10,7))                   # Set png size to 1000x700 (dpi=100)

# Create a scatter plot:
plt.scatter(lon[sel]*r2d, lat[sel]*r2d, s=sizes[sel])

# Plot outline of equatorial map:
plt.plot(lonMinMax*r2d, latMinMax*r2d, ':')


plt.axis('scaled')                                        # Set axes to a 'square grid' by moving the plot box inside the figure
plt.axis([lonMax*r2d,lonMin*r2d, latMin*r2d,latMax*r2d])  # Select Aries (RA=26-50 deg, dec=10-30 deg)
plt.xlabel(r'$\lambda_{2000}$ ($^\circ$)')                # Label the horizontal axis
plt.ylabel(r'$\beta_{2000}$ ($^\circ$)')                  # Label the vertical axis - use LaTeX for symbols

plt.tight_layout()                       # Use narrow margins
plt.savefig("ex01_ecliptical.png")       # Save the plot as png
#plt.show()                              # Show the plot to screen
plt.close()                              # Close the plot


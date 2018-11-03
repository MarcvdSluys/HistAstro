#!/bin/env python3

import numpy as np
import math

r2d = math.degrees(1)  # Radians to degrees

# Read the input file, skipping the first two lines:
hip = np.loadtxt('combihip.csv', skiprows=2, delimiter=',')

# Columns: 0: hip#, 1: vmag, 2: ra (rad), 3: dec (rad), 4: pmRA (mas/yr), 5: pmDec (mas/yr), 6: ErRA (?), 7:
# ErDec (?), 8: ErPa (mas/yr), 9: ErPd (mas/yr)


import matplotlib
matplotlib.use('Agg')  # Agg backend doesn't need an X server and is ~5x faster

import matplotlib.pyplot as plt

#plt.figure(figsize=(7,7))                   # Set png size to 1000x700 (dpi=100)

sizes = 30*(0.5 + (7.0-hip[:,1])/3.0)**2     # Scale inversely with magnitude.  Square, since .scatter() uses surface area.
ra  = hip[:,2]*r2d                           # Right ascension (deg!)
dec = hip[:,3]*r2d                           # Declination

# Create a scatter plot:
plt.scatter(ra, dec, s=sizes)

plt.axis('equal')                            # Set axes to a 'square grid' by changing the x,y limits to match image size - should go before .axis([])
plt.axis([50.0,26.0, 10.0,30.0])             # Select Aries (RA=26-50 deg, dec=10-30 deg)
plt.xlabel(r'$\alpha_{2000}$ ($^\circ$)')    # Label the horizontal axis
plt.ylabel(r'$\delta_{2000}$ ($^\circ$)')    # Label the vertical axis - use LaTeX for symbols

plt.tight_layout()                           # Use narrow margins
plt.savefig("hipparcos.png")                 # Save the plot as png
plt.close()                                  # Close the plot


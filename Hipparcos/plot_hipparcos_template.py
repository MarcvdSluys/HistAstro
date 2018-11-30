#!/bin/env python3

import math
import numpy as np
import matplotlib.pyplot as plt

r2d = math.degrees(1)  # Radians to degrees
d2r = 1.0/r2d          # Degrees to radians

# Read the input file, skipping the first two lines:
hip = np.loadtxt('newcombi.dat', skiprows=2, delimiter=',', usecols=(0,1,2,3,4,5))   # Read the numbers from columns 1-6

# Columns: 0: hip#, 1: vmag, 2: ra (rad), 3: dec (rad), 4: pmRA (mas/yr), 5: pmDec (mas/yr), 6: ErRA (?), 7:
# ErDec (?), 8: ErPa (mas/yr), 9: ErPd (mas/yr)

# Use comprehensible array names:
mag   = hip[:,1]                             # Magnitude
ra    = hip[:,2]                             # Right ascension (rad)
dec   = hip[:,3]                             # Declination (rad)
sizes = 30*(0.5 + (7.0-mag)/3.0)**2          # Scale stars inversely with magnitude.  Square, since .scatter() uses surface area.

# Set the plot boundaries:
raMin  = 26.0*d2r
raMax  = 50.0*d2r
decMin = 10.0*d2r
decMax = 30.0*d2r

# Select the stars within the boundaries:
sel = np.logical_and(ra > raMin, ra < raMax)
sel = np.logical_and(sel, dec > decMin)
sel = np.logical_and(sel, dec < decMax)
sel = ra < 1e6  # Select all stars for plotting

plt.figure(figsize=(9,7))                                 # Set png size to 900x700 (dpi=100)

# Create a scatter plot:
plt.scatter(ra[sel]*r2d, dec[sel]*r2d, s=sizes[sel])

#plt.axis('equal')                                         # Set axes to a 'square grid' by releasing the x,y limits to match image size - should come before .axis([])
plt.axis('scaled')                                        # Set axes to a 'square grid' by making the plot box smaller than the figure
plt.axis([raMax*r2d,raMin*r2d, decMin*r2d,decMax*r2d])    # Select Aries (note that RA is flipped!)
plt.xlabel(r'$\alpha_{2000}$ ($^\circ$)')                 # Label the horizontal axis
plt.ylabel(r'$\delta_{2000}$ ($^\circ$)')                 # Label the vertical axis - use LaTeX for symbols

plt.tight_layout()                                        # Use narrow margins
#plt.savefig("hipparcos.png")                              # Save the plot as png
plt.show()                                                # Plot to screen
plt.close()                                               # Close the plot

print(bla)

#!/bin/env python3

import math as m
import numpy as np
import matplotlib.pyplot as plt

r2d = m.degrees(1)  # Radians to degrees
d2r = 1.0/r2d       # Degrees to radians

# Read columns 2-4 from the input file, skipping the first line:
hip = np.loadtxt('hipparcos.csv', skiprows=1, delimiter=',', usecols=(1,2,3))

# Use comprehensible array names:
mag   = hip[:,0]                             # Magnitude
ra    = hip[:,1]                             # Right ascension (rad)
dec   = hip[:,2]                             # Declination (rad)
sizes = 30*(0.5 + (7.0-mag)/3.0)**2          # Scale stars inversely with mag.

# Set the plot boundaries:
raMin  = 26.0*d2r
raMax  = 50.0*d2r
decMin = 10.0*d2r
decMax = 30.0*d2r

# Select the stars within the boundaries:
sel = np.logical_and(ra > raMin, ra < raMax)
sel = np.logical_and(sel, dec > decMin)
sel = np.logical_and(sel, dec < decMax)
#sel = ra < 1e6  # Select all stars

plt.figure(figsize=(9,7))                  # Set png size to 900x700 (dpi=100)

# Create a scatter plot:
plt.scatter(ra[sel]*r2d, dec[sel]*r2d, s=sizes[sel])

plt.axis('scaled')                                        # Use a 'square grid'
plt.axis([raMax*r2d,raMin*r2d, decMin*r2d,decMax*r2d])    # Apply plot bounds
plt.xlabel(r'$\alpha_{2000}$ ($^\circ$)')                 # Label horiz. axis
plt.ylabel(r'$\delta_{2000}$ ($^\circ$)')                 # Label vert. axis

plt.tight_layout()               # Use narrow margins
plt.show()                       # Plot to screen
#plt.savefig("hipparcos.pdf")     # Save the plot as pdf
plt.close()                      # Close the plot

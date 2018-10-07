#!/bin/env python3

import numpy as np
import math

r2d = math.degrees(1)
r2h = r2d/15

# Skip the first two lines of the input file:
hip = np.loadtxt('combihip.csv', skiprows=2, delimiter=',')
# print(hip)


import matplotlib
matplotlib.use('Agg')  # Agg backend doesn't need an X server and is ~5x faster

import matplotlib.pyplot as plt

plt.xkcd()  # Plot everything that follows in XKCD style (needed 2x somehow)
plt.xkcd()  # Plot everything that follows in XKCD style
plt.figure(figsize=(10, 7))                    # Set png size to 1000x700 (dpi=100)

ra  = hip[:,2]*r2h                             # Right ascension
dec = hip[:,3]*r2d                             # Declination
#sizes = 1.0/(2+hip[:,1])                       # Symbol size scales inversely with magnitude - ok for whole sky
sizes = 10000.0/(2+hip[:,1])**4                # Symbol size scales inversely with magnitude for Aries

# Create a scatter plot:
plt.scatter(ra, dec, s=sizes)

# plt.xlim(24,0)                                 # Flip the x-axis range when plotting the whole sky
plt.axis([3.34,1.73, 10.0,30.0])               # Select Aries (RA=1.73-3.34h, dec=10-30 deg)
plt.xlabel(r'$\alpha$ (h)')                    # Label the horizontal axis
plt.ylabel(r'$\delta$ ($^\circ$)')             # Label the vertical axis - use LaTeX for symbols

plt.tight_layout()                             # Use narrow margins
# plt.show()                                     # Show the plot to screen
plt.savefig("hipparcos.png")                   # Save the plot as png
plt.close()                                    # Close the plot


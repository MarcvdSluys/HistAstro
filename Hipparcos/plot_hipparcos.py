#!/bin/env python3

def eq2ecl(ra,dec, eps):
    import numpy as np
    lon = np.arctan2( np.sin(ra)  * np.cos(eps) + np.tan(dec) * np.sin(eps),  np.cos(ra) )
    lat =  np.arcsin( np.sin(dec) * np.cos(eps) - np.cos(dec) * np.sin(eps) * np.sin(ra) )
    return lon,lat

import numpy as np
import math

import time

t0 = time.perf_counter() 


r2d = math.degrees(1)
#r2h = r2d/15
d2r = 1.0/r2d

# Read the input file, skipping the first two lines:
#hip = np.loadtxt('combihip.csv', skiprows=2, delimiter=',')  # Works (old file, no text)
#hip = np.loadtxt('combihip.csv', skiprows=2, delimiter=',', usecols=(0,1,2,3,4,5))  # Works (old file, no text)
#hip = np.genfromtxt('newcombi.dat', skip_header=1, delimiter=',')  # WORKS, but text fields become nan
#hip = np.genfromtxt('newcombi.dat', skip_header=1, delimiter=',', dtype=None)  # WORKS, but get hip[15544][13] iso hip[15544,13] 

t1 = time.perf_counter() 
hip    = np.loadtxt('newcombi.dat', skiprows=2, delimiter=',', usecols=(0,1,2,3,4,5))             # Read the numbers
hiptxt = np.loadtxt('newcombi.dat', skiprows=2, delimiter=',', usecols=(10,11,12), dtype=np.str)  # Read the text columns
t2 = time.perf_counter() 

# Columns: 0: hip#, 1: vmag, 2: ra (rad), 3: dec (rad), 4: pmRA (mas/yr), 5: pmDec (mas/yr), 6: ErRA (?), 7:
# ErDec (?), 8: ErPa (mas/yr), 9: ErPd (mas/yr)

#hip = hip.reshape((15544,13))
#print(type(hip))
#print(hip.shape)
#print(hip[1])
#print(hiptxt[1])
#print(hip[1][11])

sizes = 30*(0.5 + (7.0-hip[:,1])/3.0)**2     # Scale inversely with magnitude.  Square, since scatter() uses surface area
#ra  = hip[:,2]*r2h                          # Right ascension (h)
ra  = hip[:,2]*r2d                           # Right ascension (deg)
dec = hip[:,3]*r2d                           # Declination
pma = hip[:,4]/3.6e6                         # pmRA, mas/yr -> deg/yr
pmd = hip[:,5]/3.6e6                         # pmDec, mas/yr -> deg/yr


# Correct for proper motion:
startEpoche = 1992.25
targetEpoche = -1500
dt = targetEpoche - startEpoche

raOld  = ra  + pma*dt / np.cos(dec*d2r)
decOld = dec + pmd*dt


# Convert to ecliptical coordinates:
#import astropy.coordinates as coord
#import astropy.units as u
#hip1 = coord.ICRS(ra=ra*u.degree, dec=dec*u.degree, pm_ra_cosdec=pma*u.degree/u.yr, pm_dec=pmd*u.degree/u.yr)
#print(hip1)


eps = 0.40931975  # For 2000 in rad
lon,lat = eq2ecl(ra,dec, eps)

raMin  = 26.0
raMax  = 50.0
decMin = 10.0
decMax = 30.0

sel = np.logical_and(ra > raMin, ra < raMax)
sel = np.logical_and(sel, dec > decMin)
sel = np.logical_and(sel, dec < decMax)
#sel = ra < 1e6  # Select all stars for plotting

t3 = time.perf_counter() 



import matplotlib
matplotlib.use('Agg')  # Agg backend doesn't need an X server and is ~5x faster

import matplotlib.pyplot as plt


# Plot equatorial:
plt.xkcd()  # Plot everything that follows in XKCD style (needed 2x somehow)
plt.xkcd()  # Plot everything that follows in XKCD style
plt.figure(figsize=(10,7))                   # Set png size to 1000x700 (dpi=100)

# Create a scatter plot:
t4 = time.perf_counter() 
plt.scatter(ra[sel], dec[sel], s=sizes[sel])

#plt.scatter(raOld[sel], decOld[sel], s=sizes[sel])
t5 = time.perf_counter() 

#plt.xlim(24,0)                              # Flip the x-axis range when plotting the whole sky
#plt.axis('equal')                            # Set axes to a 'square grid' by changing the x,y limits to match image size - should go before .axis([])
plt.axis('scaled')                          # Set axes to a 'square grid' by moving the plot box inside the figure
#plt.axis('square')                          # Set axes to a 'square grid' by moving the plot box inside the figure and setting xmax-xmin = ymax-ymin
plt.axis([raMax,raMin, decMin,decMax])             # Select Aries (RA=26-50 deg, dec=10-30 deg)
plt.xlabel(r'$\alpha_{2000}$ ($^\circ$)')           # Label the horizontal axis
plt.ylabel(r'$\delta_{2000}$ ($^\circ$)')           # Label the vertical axis - use LaTeX for symbols

plt.tight_layout()                           # Use narrow margins
plt.savefig("hipparcos_equatorial.pdf")                 # Save the plot as png
plt.close()                                  # Close the plot









t6 = time.perf_counter()
t9 = time.perf_counter()

print("Total run time:  %0.2f s" % (t9-t0))
print("  read file:     %0.2f s" % (t2-t1))
print("  total plot:    %0.2f s" % (t6-t3))
print("    init:        %0.2f s" % (t4-t3))
print("    scatter:     %0.2f s" % (t5-t4))
print("    finish:      %0.2f s" % (t6-t5))



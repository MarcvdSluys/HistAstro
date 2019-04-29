#!/bin/env python3

# Modules:
import math as m
import numpy as np
import matplotlib.pyplot as plt

import histastro.coordinates as coord
import histastro.datetime as dt



# Constants:
r2d   = m.degrees(1)  # Radians to degrees
d2r   = 1.0/r2d       # Degrees to radians
h2r   = d2r*15        # Hours to radians
as2r  = d2r/3.6e3     # Arcseconds to radians
mas2r = as2r/1000.0   # Milliarcseconds to radians


    
# Read the Hipparcos catalogue
# Columns: 0: hip#, 1: vmag, 2: ra (rad), 3: dec (rad), 4: pmRA (mas/yr), 5: pmDec (mas/yr), 6: ErRA (?), 7:
# ErDec (?), 8: ErPa (mas/yr), 9: ErPd (mas/yr)

hip = np.loadtxt('data/BrightStars.csv', skiprows=1, delimiter=',', usecols=(1,2,3,4,5))    # Read columns 1-5: V, ra, dec, pma, pmd

mag = hip[:,0]
ra  = hip[:,1]                           # Right ascension (rad)
dec = hip[:,2]                           # Declination (rad)
pma = hip[:,3]*mas2r                     # pmRA, mas/yr -> rad/yr
pmd = hip[:,4]*mas2r                     # pmDec, mas/yr -> rad/yr



# Verify Julian days (ex. 3.1):
print("Julian days:")
print(-3000, dt.julianDay(-3000, 1, 1.5))
print( 1000, dt.julianDay(1000,  1, 1.5))
print( 2000, dt.julianDay(2000,  1, 1.5))



# Prepare horizontal maps:
phi = 51.178*d2r

ha = -6*h2r - ra  # At the spring equinox and sunrise ra_sun = 0, ha_sun=-6h
az,alt = coord.par2horiz(ha,dec, phi)

jdHip = dt.julianDay(1991, 4, 1)

azMin  = 225*d2r
azMax  = 305*d2r
altMin = 0*d2r
altMax = 40*d2r

Mlim = 4.5  # Magnitude limit
sizes = 20*(0.5 + (Mlim-mag)/3.0)**2     # Scale inversely with magnitude.  Square, since scatter() uses surface area




# See the effect of proper motion (ex. 3.2):
year = -10000
jd = dt.julianDay(year, 1, 1)
raOld,decOld = coord.properMotion(jdHip,jd, ra,dec, pma,pmd)

# Convert equatorial -> parallactic -> horizontal coordinates:
haOld = -6*h2r - raOld  # At the spring equinox and sunrise ra_sun = 0, ha_sun=-6h
azOld,altOld = coord.par2horiz(haOld,decOld, phi)

sel = np.logical_and(az > azMin, az < azMax)
sel = np.logical_and(sel, alt > altMin)
sel = np.logical_and(sel, alt < altMax)
sel = np.logical_and(sel, mag < Mlim)
#sel = az < 1e6  # Select all stars for plotting


# Create scatter plots:
plt.figure(figsize=(10,5.5))                   # Set png size to 1000x550 (dpi=100)
plt.scatter(az[sel]*r2d, alt[sel]*r2d, s=sizes[sel])
plt.scatter(azOld[sel]*r2d, altOld[sel]*r2d, s=sizes[sel], alpha=1.)

plt.axis('scaled')                          # Set axes to a 'square grid' by moving the plot box inside the figure
plt.axis([azMin*r2d,azMax*r2d, altMin*r2d,altMax*r2d])             # Select Aries (RA=26-50 deg, dec=10-30 deg)
plt.xlabel(r'A ($^\circ$)')           # Label the horizontal axis
plt.ylabel(r'h ($^\circ$)')           # Label the vertical axis - use LaTeX for symbols

plt.tight_layout()                           # Use narrow margins
plt.savefig("ex03_properMotion_-10000.png")  # Save the plot as png
plt.close()                                  # Close the plot










# See the effect of precession (ex. 3.3):
year = 1000
jd = dt.julianDay(year, 1, 1)
raOld,decOld = coord.precessHip(jd, ra,dec)

# Convert equatorial -> parallactic -> horizontal coordinates:
haOld = -6*h2r - raOld  # At the spring equinox and sunrise ra_sun = 0, ha_sun=-6h
azOld,altOld = coord.par2horiz(haOld,decOld, phi)

sel = np.logical_and(az > azMin, az < azMax)
sel = np.logical_and(sel, alt > altMin)
sel = np.logical_and(sel, alt < altMax)
sel = np.logical_and(sel, mag < Mlim)
#sel = az < 1e6  # Select all stars for plotting

plt.figure(figsize=(10,5.5))                   # Set png size to 1000x550 (dpi=100)

# Create scatter plots:
plt.scatter(az[sel]*r2d, alt[sel]*r2d, s=sizes[sel])
plt.scatter(azOld[sel]*r2d, altOld[sel]*r2d, s=sizes[sel], alpha=1.)


plt.axis('scaled')                          # Set axes to a 'square grid' by moving the plot box inside the figure
plt.axis([azMin*r2d,azMax*r2d, altMin*r2d,altMax*r2d])             # Select Aries (RA=26-50 deg, dec=10-30 deg)
plt.xlabel(r'A ($^\circ$)')           # Label the horizontal axis
plt.ylabel(r'h ($^\circ$)')           # Label the vertical axis - use LaTeX for symbols

plt.tight_layout()                           # Use narrow margins
plt.savefig("ex03_precession_1000.png")      # Save the plot as png
plt.close()                                  # Close the plot








# Does the bear dip in the ocean in Athens in 800BC? - correct for proper motion and precession (ex 3.4)
year = -800  
jd = dt.julianDay(year, 1, 1)
raOld,decOld = coord.properMotion(jdHip,jd, ra,dec, pma,pmd)
raOld,decOld = coord.precessHip(jd, raOld,decOld)

rMax  = 60*d2r  # Plot limit (radius)
Mlim  = 5.0     # Magnitude limit
sizes = 20*(0.5 + (Mlim-mag)/3.0)**2     # Scale inversely with magnitude.  Square, since scatter() uses surface area


# Plot polar equatorial map:
plt.figure(figsize=(10,10))                   # Set png size to 1000x1000 (dpi=100)
ax = plt.subplot(111, projection='polar')

# Compute r and theta from ra and dec:
r = m.pi/2 - dec
theta = -ra

# Compute old r and theta from old ra and dec:
rOld = m.pi/2 - decOld  # Ensure raOld, decOld are for 800 BCE
thetaOld = -raOld

# Select bright stars close to the NP:
sel    = (r < rMax) & (mag < Mlim)
selOld = (rOld < rMax) & (mag < Mlim)

# Make scatter plots:
ax.scatter(theta[sel],        r[sel]*r2d,        s=sizes[sel])
ax.scatter(thetaOld[selOld],  rOld[selOld]*r2d,  s=sizes[selOld])

# Draw a circle at 38deg around NP to indicate the horizon in Athens:
rCirc = np.ones(101)*38
thCirc = np.arange(101)/100*m.pi*2
ax.plot(thCirc, rCirc, 'r')

ax.set_ylim(0, rMax*r2d)

plt.tight_layout()                             # Use narrow margins
plt.savefig("ex03_Bear.png")                   # Save the plot as png
plt.close()                                    # Close the plot


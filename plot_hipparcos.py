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
#r2h  = r2d/15
h2r   = d2r*15        # Hours to radians
as2r  = d2r/3.6e3     # Arcseconds to radians
mas2r = as2r/1000.0   # Milliarcseconds to radians


    


# Time the execution:
import time
t0 = time.perf_counter() 

# Read the input file, skipping the first two lines:
#hip = np.loadtxt('combihip.csv', skiprows=2, delimiter=',')  # Works (old file, no text)
#hip = np.loadtxt('combihip.csv', skiprows=2, delimiter=',', usecols=(0,1,2,3,4,5))  # Works (old file, no text)
#hip = np.genfromtxt('BrightStars.csv', skip_header=1, delimiter=',')  # WORKS, but text fields become nan
#hip = np.genfromtxt('BrightStars.csv', skip_header=1, delimiter=',', dtype=None)  # WORKS, but get hip[15544][13] iso hip[15544,13] 

t1 = time.perf_counter() 
hip    = np.loadtxt('BrightStars.csv', skiprows=1, delimiter=',', usecols=(0,1,2,3,4,5))             # Read the numbers in columns 0-5: HIP, V, ra, dec, pma, pmd - amazingly, this can also read BrightStars.csv.gz!
#hiptxt = np.loadtxt('BrightStars.csv', skiprows=1, delimiter=',', usecols=(7,10,11), dtype=np.str)  # Read the text columns
t2 = time.perf_counter() 

# Columns: 0: hip#, 1: vmag, 2: ra (rad), 3: dec (rad), 4: pmRA (mas/yr), 5: pmDec (mas/yr), 6: ErRA (?), 7:
# ErDec (?), 8: ErPa (mas/yr), 9: ErPd (mas/yr)

#hip = hip.reshape((15544,13))
#print(type(hip))
#print(hip.shape)
#print(hip[1])
#print(hiptxt[1])
#print(hip[1][11])

Mlim = 7.0  # Magnitude limit
mag = hip[:,1]
sizes = 30*(0.5 + (Mlim-mag)/3.0)**2     # Scale inversely with magnitude.  Square, since scatter() uses surface area
ra  = hip[:,2]                           # Right ascension (rad)
dec = hip[:,3]                           # Declination (rad)
pma = hip[:,4]*mas2r                     # pmRA, mas/yr -> rad/yr
pmd = hip[:,5]*mas2r                     # pmDec, mas/yr -> rad/yr


# Correct for proper motion:
t3 = time.perf_counter()
jd1 = dt.julianDay( 1991, 4, 1)
jd2 = dt.julianDay(-1500, 1, 1)
raOld,decOld = coord.properMotion(jd1,jd2, ra,dec, pma,pmd)

#print(-3000, dt.julianDay(-3000, 1, 1.5))
#print(1000, dt.julianDay(1000, 1, 1.5))
#print(2000, dt.julianDay(2000, 1, 1.5))

t4 = time.perf_counter() 
raMin  = 26.0*d2r
raMax  = 50.0*d2r
decMin = 10.0*d2r
decMax = 30.0*d2r

sel = np.logical_and(ra > raMin, ra < raMax)
sel = np.logical_and(sel, dec > decMin)
sel = np.logical_and(sel, dec < decMax)
#sel = ra < 1e6  # Select all stars for plotting




#import matplotlib
#matplotlib.use('Agg')  # Agg backend doesn't need an X server and is ~5x faster



# Plot equatorial map:
#plt.style.use('dark_background')
plt.figure(figsize=(10,7))                   # Set png size to 1000x700 (dpi=100)

# Make a scatter plot.  s contains the *surface areas* of the circles:
plt.scatter(ra[sel]*r2d, dec[sel]*r2d, s=sizes[sel])

plt.scatter(raOld[sel]*r2d, decOld[sel]*r2d, s=sizes[sel])

#plt.xlim(24,0)                              # Flip the x-axis range when plotting the whole sky
#plt.axis('equal')                            # Set axes to a 'square grid' by changing the x,y limits to match image size - should go before .axis([])
plt.axis('scaled')                          # Set axes to a 'square grid' by moving the plot box inside the figure
#plt.axis('square')                          # Set axes to a 'square grid' by moving the plot box inside the figure and setting xmax-xmin = ymax-ymin
plt.axis([raMax*r2d,raMin*r2d, decMin*r2d,decMax*r2d])             # Select Aries (RA=26-50 deg, dec=10-30 deg)
plt.xlabel(r'$\alpha_{2000}$ ($^\circ$)')           # Label the horizontal axis
plt.ylabel(r'$\delta_{2000}$ ($^\circ$)')           # Label the vertical axis - use LaTeX for symbols

plt.tight_layout()                           # Use narrow margins
plt.savefig("hipparcos_equatorial.png")                 # Save the plot as png
#plt.show()                              # Show the plot to screen
plt.close()                                  # Close the plot



# Convert to ecliptic coordinates:
t5 = time.perf_counter() 
eps = 0.40909280  # For 2000 in rad
lon,lat = coord.eq2ecl(ra,dec, eps)


# I use five coordinates for the corners of the map, in order to plot four lines below.
lonMinMax,latMinMax = coord.eq2ecl([raMin,raMax,raMax,raMin,raMin],[decMin,decMin,decMax,decMax,decMin], eps)
lonMin = min(lonMinMax)
lonMax = max(lonMinMax)
latMin = min(latMinMax)
latMax = max(latMinMax)

sel = np.logical_and(lon > lonMin, lon < lonMax)
sel = np.logical_and(sel, lat > latMin)
sel = np.logical_and(sel, lat < latMax)
#sel = lon < 1e6  # Select all stars for plotting






# Plot ecliptic map:
plt.figure(figsize=(10,7))                   # Set png size to 1000x700 (dpi=100)

# Create a scatter plot:
plt.scatter(lon[sel]*r2d, lat[sel]*r2d, s=sizes[sel])

plt.plot(lonMinMax*r2d, latMinMax*r2d, ':')


plt.axis('scaled')                          # Set axes to a 'square grid' by moving the plot box inside the figure
plt.axis([lonMax*r2d,lonMin*r2d, latMin*r2d,latMax*r2d])             # Select Aries (RA=26-50 deg, dec=10-30 deg)
plt.xlabel(r'$\lambda_{2000}$ ($^\circ$)')           # Label the horizontal axis
plt.ylabel(r'$\beta_{2000}$ ($^\circ$)')           # Label the vertical axis - use LaTeX for symbols

plt.tight_layout()                           # Use narrow margins
plt.savefig("hipparcos_ecliptic.png")                 # Save the plot as png
#plt.show()                              # Show the plot to screen
plt.close()                                  # Close the plot

t6 = time.perf_counter()




# Plot horizontal map:
phi = 51.178*d2r

ha = -6*h2r - ra  # At the spring equinox and sunrise ra_sun = 0, ha_sun=-6h
az,alt = coord.par2horiz(ha,dec, phi)

# Correct for proper motion and precession:
#year = -10000
#year = 1000
year = -800  # Does the bear dip in the ocean in Athens in 800BC?
jd1 = dt.julianDay(1991, 4, 1)
jd2 = dt.julianDay(year, 1, 1)
raOld,decOld = coord.properMotion(jd1,jd2, ra,dec, pma,pmd)
raOld,decOld = coord.precessHip(jd2, raOld,decOld)

haOld = -6*h2r - raOld  # At the spring equinox and sunrise ra_sun = 0, ha_sun=-6h
azOld,altOld = coord.par2horiz(haOld,decOld, phi)


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
sel = mag < Mlim  # Magnitude limit only
#plt.scatter(azOld[sel]*r2d, altOld[sel]*r2d, s=sizes[sel], alpha=1.)


plt.axis('scaled')                          # Set axes to a 'square grid' by moving the plot box inside the figure
plt.axis([azMin*r2d,azMax*r2d, altMin*r2d,altMax*r2d])             # Select Aries (RA=26-50 deg, dec=10-30 deg)
plt.xlabel(r'A ($^\circ$)')           # Label the horizontal axis
plt.ylabel(r'h ($^\circ$)')           # Label the vertical axis - use LaTeX for symbols

plt.tight_layout()                           # Use narrow margins
plt.savefig("hipparcos_horizontal.png")                 # Save the plot as png
#plt.show()                              # Show the plot to screen
plt.close()                                  # Close the plot

t7 = time.perf_counter()










# Plot polar equatorial map:
plt.figure(figsize=(7,7))                   # Set png size to 1000x700 (dpi=100)
ax = plt.subplot(111, projection='polar')

# Compute r and theta from ra and dec:
r = m.pi/2 - dec
theta = -ra

# Compute r and theta from precessed ra and dec:
r = m.pi/2 - decOld  # Ensure raOld, decOld are for 800 BCE
theta = -raOld

rMax = 60*d2r  # Plot limit
Mlim = 5.0  # Magnitude limit
sizes = 20*(0.5 + (Mlim-mag)/3.0)**2     # Scale inversely with magnitude.  Square, since scatter() uses surface area

# Select bright stars close to the NP:
sel = (r < rMax) & (mag < Mlim)

# Make a scatter plot.  s contains the *surface areas* of the circles:
ax.scatter(theta[sel],    r[sel]*r2d,    s=sizes[sel])

# Draw a circle at 38deg around NP to indicate the horizon in Athens:
rCirc = np.ones(101)*38
thCirc = np.arange(101)/100*m.pi*2
ax.plot(thCirc, rCirc, 'r')

ax.set_ylim(0,rMax*r2d)

plt.tight_layout()                             # Use narrow margins
plt.savefig("hipparcos_equatorial_polar.png")  # Save the plot as png
#plt.show()                                     # Show the plot to screen
plt.close()                                    # Close the plot










t9 = time.perf_counter()

print("Total run time:  %0.2f s" % (t9-t0))
print("  read file:     %0.2f s" % (t2-t1))
print("  pm:            %0.2f s" % (t4-t3))
print("  plot eq:       %0.2f s" % (t5-t4))
print("  plot ecl:      %0.2f s" % (t6-t5))
print("  plot az:       %0.2f s" % (t7-t6))



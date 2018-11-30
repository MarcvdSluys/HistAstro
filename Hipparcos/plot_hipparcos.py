#!/bin/env python3

# Modules:
import math
import numpy as np
import matplotlib.pyplot as plt

# Constants:
pi    = math.pi
pi2   = 2*pi
r2d   = math.degrees(1)  # Radians to degrees
#r2h  = r2d/15
d2r   = 1.0/r2d          # Degrees to radians
mas2r = d2r/3.6e6      # Milliarcseconds to radians
h2r   = d2r*15          # Hours to radians

def eq2ecl(ra,dec, eps):
    lon = np.arctan2( np.sin(ra)  * np.cos(eps) + np.tan(dec) * np.sin(eps),  np.cos(ra) )
    lat =  np.arcsin( np.sin(dec) * np.cos(eps) - np.cos(dec) * np.sin(eps) * np.sin(ra) )
    return lon,lat

def par2horiz(ha,dec, phi):
    az  = np.arctan2( np.sin(ha),   np.cos(ha) * np.sin(phi) - np.tan(dec) * np.cos(phi) ) % pi2
    alt = np.arcsin(  np.sin(dec) * np.sin(phi) + np.cos(ha) * np.cos(dec) * np.cos(phi) )
    return az,alt

    
# Time the execution:
import time
t0 = time.perf_counter() 

xkcd = False

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

Mlim = 7.0  # Magnitude limit
mag = hip[:,1]
sizes = 30*(0.5 + (Mlim-mag)/3.0)**2     # Scale inversely with magnitude.  Square, since scatter() uses surface area
ra  = hip[:,2]                           # Right ascension (rad)
dec = hip[:,3]                           # Declination (rad)
pma = hip[:,4]*mas2r                     # pmRA, mas/yr -> rad/yr
pmd = hip[:,5]*mas2r                     # pmDec, mas/yr -> rad/yr


# Correct for proper motion:
t3 = time.perf_counter() 
startEpoche  = 1992.25
targetEpoche = -1500
dt = targetEpoche - startEpoche

raOld  = ra  + pma*dt / np.cos(dec)
decOld = dec + pmd*dt


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



# Plot equatorial:
if xkcd:
    plt.xkcd()  # Plot everything that follows in XKCD style (needed 2x somehow)
    plt.xkcd()  # Plot everything that follows in XKCD style

#plt.style.use('dark_background')
plt.figure(figsize=(10,7))                   # Set png size to 1000x700 (dpi=100)

# Create a scatter plot:
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




t5 = time.perf_counter() 
eps = 0.40931975  # For 2000 in rad
lon,lat = eq2ecl(ra,dec, eps)


lonMinMax,latMinMax = eq2ecl([raMin,raMax,raMax,raMin,raMin],[decMin,decMin,decMax,decMax,decMin], eps)
lonMin = min(lonMinMax)
lonMax = max(lonMinMax)
latMin = min(latMinMax)
latMax = max(latMinMax)

sel = np.logical_and(lon > lonMin, lon < lonMax)
sel = np.logical_and(sel, lat > latMin)
sel = np.logical_and(sel, lat < latMax)
#sel = lon < 1e6  # Select all stars for plotting






# Plot equatorial:
if xkcd:
    plt.xkcd()  # Plot everything that follows in XKCD style (needed 2x somehow)
    plt.xkcd()  # Plot everything that follows in XKCD style
    
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




# Horizontal map:
phi = 51.178*d2r

ha = -6*h2r - ra  # At the spring equinox and sunrise ra_sun = 0, ha_sun=-6h
az,alt = par2horiz(ha,dec, phi)

if xkcd:
    plt.xkcd()  # Plot everything that follows in XKCD style (needed 2x somehow)
    plt.xkcd()  # Plot everything that follows in XKCD style
    
plt.figure(figsize=(10,3))                   # Set png size to 1000x700 (dpi=100)

azMin  = 225*d2r
azMax  = 305*d2r
altMin = 0*d2r
altMax = 20*d2r

Mlim = 4.5  # Magnitude limit
sizes = 20*(0.5 + (Mlim-mag)/3.0)**2     # Scale inversely with magnitude.  Square, since scatter() uses surface area

sel = np.logical_and(az > azMin, az < azMax)
sel = np.logical_and(sel, alt > altMin)
sel = np.logical_and(sel, alt < altMax)
sel = np.logical_and(sel, mag < Mlim)
#sel = az < 1e6  # Select all stars for plotting

# Create a scatter plot:
plt.scatter(az[sel]*r2d, alt[sel]*r2d, s=sizes[sel])


plt.axis('scaled')                          # Set axes to a 'square grid' by moving the plot box inside the figure
plt.axis([azMin*r2d,azMax*r2d, altMin*r2d,altMax*r2d])             # Select Aries (RA=26-50 deg, dec=10-30 deg)
plt.xlabel(r'A ($^\circ$)')           # Label the horizontal axis
plt.ylabel(r'h ($^\circ$)')           # Label the vertical axis - use LaTeX for symbols

plt.tight_layout()                           # Use narrow margins
plt.savefig("hipparcos_horizontal.png")                 # Save the plot as png
#plt.show()                              # Show the plot to screen
plt.close()                                  # Close the plot

t7 = time.perf_counter()














t9 = time.perf_counter()

print("Total run time:  %0.2f s" % (t9-t0))
print("  read file:     %0.2f s" % (t2-t1))
print("  pm:            %0.2f s" % (t4-t3))
print("  plot eq:       %0.2f s" % (t5-t4))
print("  plot ecl:      %0.2f s" % (t6-t5))
print("  plot az:       %0.2f s" % (t7-t6))



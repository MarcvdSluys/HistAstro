#!/bin/env python3

import numpy as np
import fortranformat as ff

def readVSOP(fileName):
    inFile = open(fileName,'r')
    
    formatHeader = ff.FortranRecordReader('(40x,I3, 16x,I1,I8)')         # Block header format
    formatBody   = ff.FortranRecordReader('(79x,F18.11,F14.11,F20.11)')  # Block body format
    
    lon=[]; lat=[]; rad=[]
    
    for iBlock in range(3*6):  # 3 variables (l,b,r), up to 6 powers (0-5)
        line = inFile.readline()
        var,power,nTerm = formatHeader.read(line)
        #print(var,power,nTerm)
        if line == '': break  # EoF
        
        for iLine in range(nTerm):
            line = inFile.readline()
            a,b,c = formatBody.read(line)
            #print(iLine, var,power, a,b,c)
            
            if var == 1: lon.append([power, a,b,c])  # var=1: ecliptic longitude
            if var == 2: lat.append([power, a,b,c])  # var=2: ecliptic latitude
            if var == 3: rad.append([power, a,b,c])  # var=3: radial distance

    return lon,lat,rad


# Read a VSOP data file:
lon,lat,rad = readVSOP('VSOP87D.jup')
#lon,lat,rad = readVSOP('VSOP87D.ura')

print(np.shape(lon))
print(lon[0])
print(lon[-1])

print(np.shape(lat))
print(lat[0])
print(lat[-1])

print(np.shape(rad))
print(rad[0])
print(rad[-1])

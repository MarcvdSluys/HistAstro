#!/bin/env python3

import math as m
import numpy as np
import fortranformat as ff

def readVSOP(fileName):
    inFile = open(fileName,'r')
    
    formatHeader = ff.FortranRecordReader('(40x,I3, 16x,I1,I8)')         # Block header format
    formatBody   = ff.FortranRecordReader('(79x,F18.11,F14.11,F20.11)')  # Block body format
    
    lonTerms=[]; latTerms=[]; radTerms=[]
    
    for iBlock in range(3*6):  # 3 variables (l,b,r), up to 6 powers (0-5)
        line = inFile.readline()
        var,power,nTerm = formatHeader.read(line)
        #print(var,power,nTerm)
        if line == '': break  # EoF
        
        for iLine in range(nTerm):
            line = inFile.readline()
            a,b,c = formatBody.read(line)
            #print(iLine, var,power, a,b,c)
            
            if var == 1: lonTerms.append([power, a,b,c])  # var=1: ecliptic longitude
            if var == 2: latTerms.append([power, a,b,c])  # var=2: ecliptic latitude
            if var == 3: radTerms.append([power, a,b,c])  # var=3: radial distance

    return lonTerms,latTerms,radTerms


def computeLBR(JDE, lonTerms,latTerms,radTerms):
    tau = (JDE - 2451545.0)/365250.0
    
    lon=0.0; lat=0.0; rad=0.0
    for terms in lonTerms:
        cosTerm = terms[1] * m.cos(terms[2] + terms[3]*tau)
        lon += cosTerm * tau**terms[0]
    for terms in latTerms:
        cosTerm = terms[1] * m.cos(terms[2] + terms[3]*tau)
        lat += cosTerm * tau**terms[0]
    for terms in radTerms:
        cosTerm = terms[1] * m.cos(terms[2] + terms[3]*tau)
        rad += cosTerm * tau**terms[0]
    return lon,lat,rad


# Read a VSOP data files for Earth and desired planet:
lonTermsE,latTermsE,radTermsE = readVSOP('VSOP87D.ear')
lonTerms,latTerms,radTerms = readVSOP('VSOP87D.jup')
#lonTerms,latTerms,radTerms = readVSOP('VSOP87D.ura')

print(np.shape(lonTerms))
print(lonTerms[0])
print(lonTerms[-1])

print(np.shape(latTerms))
print(latTerms[0])
print(latTerms[-1])

print(np.shape(radTerms))
print(radTerms[0])
print(radTerms[-1])


# Compute heliocentric position:
JDE = 2451545.0

# Earth:
HClonE,HClatE,HCradE = computeLBR(JDE, lonTermsE,latTermsE,radTermsE)
print(HClonE,HClatE,HCradE)

# Planet:
HClon,HClat,HCrad = computeLBR(JDE, lonTerms,latTerms,radTerms)
print(HClon,HClat,HCrad)


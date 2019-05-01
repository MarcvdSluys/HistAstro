"""HistAstro Moon functions"""

import math as m
import numpy.core as np
import sys

# Constants:
d2r = m.radians(1)  # Degrees to radians
r2d = m.degrees(1)  # Radians to degrees
r2as = r2d*3600     # Radians to arcseconds
pi = m.pi
pi2 = pi*2
pio2 = pi/2
jd2000 = 2451545

# Global variable
modeInit = 999  # := uninitialised

#***************************************************************************************************
def elp_mpp02_initialise_and_read_files(mode):
    global modeInit
    print("Initialise and read files:", modeInit)
    #integer, save :: modeInit = 999       # Test whether data have been initialised
    
    # Initializing of constants and reading the files:
    ierr = 0
    #print*,mode,modeInit
    if(mode!=modeInit):
        w,eart,peri, zeta,dela,  p,delnu,dele,delg,delnp,delep,dtasm,am,  p1,p2,p3,p4,p5, q1,q2,q3,q4,q5 = elp_mpp02_initialise(mode)
        ierr = elp_mpp02_read_files(w,eart,peri, zeta,dela,  p,delnu,dele,delg,delnp,delep,dtasm,am,  p1,p2,p3,p4,p5, q1,q2,q3,q4,q5)
        print("init_and_read ierr:", ierr)
        if(ierr!=0): return ierr
        
        modeInit = mode
        print("modeInit:", modeInit)
        return ierr,  w,eart,peri, zeta,dela,  p,delnu,dele,delg,delnp,delep,dtasm,am,  p1,p2,p3,p4,p5, q1,q2,q3,q4,q5
        
    return ierr
#***************************************************************************************************



#***************************************************************************************************
#> \brief  Initialization of the constants and parameters used for the evaluation of the ELP/MPP02 series
##
## \param mode  Index of the corrections to the constants: 0: LLR observations for 1970-2000, 1: DE405 ephemeris for 1950-2060
##
## \retval elp_mpp02_constants  Set of the constants of ELPMPP02 solution (module)
##
## \note
##
## - Remarks:
##    The nominal values of some constants have to be corrected.  There are two sets of corrections, one of which can be chosen
##    using the parameter 'mode' (used in elp_mpp02_initialise()):
##    - mode=0, the constants are fitted to LLR observations provided from 1970 to 2001; it is the default value;
##    - mode=1, the constants are fitted to DE405 ephemeris over one century (1950-2060); the lunar angles W1, W2, W3 receive also additive corrections to the secular coefficients.
## 
## - Moon constants:
##    nu        : mean motion of the Moon (W1(1,1))                 (Nu)
##    g         : half coefficient of sin(F) in latitude         (Gamma)
##    e         : half coefficient of sin(l) in longitude            (E)
##    np        : mean motion of EMB (eart(1))                      (n')
##    ep        : eccentricity of EMB                               (e')
##
##    p is the precession rate and t is the time
##
## \see
##   - ELPdoc: Lunar solution ELP, version ELP/MPP02,  Jean Chapront and Gerard Francou, October 2002

def elp_mpp02_initialise(mode):
    global modeInit
    print("Initialise:", modeInit)
    #use TheSky_elp_mpp02_constants, only: w,eart,peri, zeta,dela,   p,delnu,dele,delg,delnp,delep,dtasm,am,   p1,p2,p3,p4,p5,q1,q2,q3,q4,q5
    
    #real(double) :: bp(5,2)
    
    # Constant for the correction to the constant of precession:
    Dprec = -0.29965                                  # Source: IAU 2000A
    
    
    #bp = np.array([[0.311079095,-0.4482398e-2], [-0.110248500e-2,0.1056062e-2], [0.50928e-4,-0.103837907],
    #               [0.6682870e-3,-0.129807200e-2], [-0.1780280e-3,-0.37342e-4]]) # (5,2)
    bp = np.array([[0,0,0], [0,0.311079095,-0.4482398e-2], [0,-0.110248500e-2,0.1056062e-2], [0,0.50928e-4,-0.103837907],
                   [0,0.6682870e-3,-0.129807200e-2], [0,-0.1780280e-3,-0.37342e-4]])  # (5,2) -> (6,3) with first row/column 0
    
    if(mode<0 or mode>1): sys.exit('elp_mpp02_initialise(): mode must have value 0 or 1, not %i' % mode)
    
    # Constants for the evaluation of the partial derivatives:
    am     =  0.074801329           # Ratio of the mean motions (EMB / Moon)
    alpha  =  0.002571881           # Ratio of the semi-major axis (Moon / EMB)
    dtasm  =  (2*alpha)/(3*am)  # (2*alpha) / (3*am)
    xa     =  (2*alpha)/3
    
    
    # Corrections to constants:
    if(mode==0):  # Default - LLR
       # Values of the corrections to the constants fitted to LLR.  Fit 13-05-02 (2 iterations) except Phi and eps w2_1 and w3_1
       # See ELPdoc, Table 3 and paper, Table 1
       Dw1_0   = -0.10525
       Dw2_0   =  0.16826
       Dw3_0   = -0.10760
       Deart_0 = -0.04012
       Dperi   = -0.04854
       Dw1_1   = -0.32311
       Dgam    =  0.00069
       De      =  0.00005
       Deart_1 =  0.01442
       Dep     =  0.00226
       Dw2_1   =  0.08017
       Dw3_1   = -0.04317
       Dw1_2   = -0.03794
    else:  # DE 405
       # Values of the corrections to the constants fitted to DE405 over the time interval (1950-2060)
       Dw1_0   = -0.07008
       Dw2_0   =  0.20794
       Dw3_0   = -0.07215
       Deart_0 = -0.00033
       Dperi   = -0.00749
       Dw1_1   = -0.35106
       Dgam    =  0.00085
       De      = -0.00006
       Deart_1 =  0.00732
       Dep     =  0.00224
       Dw2_1   =  0.08017
       Dw3_1   = -0.04317
       Dw1_2   = -0.03743

    # Fundamental arguments (Moon and EMB - ELPdoc, Table 1):
    w = np.zeros((4,5))  # Actually, need [1:3,0:4]
    
    # W1: mean longitude of the Moon:
    w[1,0]  = elp_dms2rad(218,18,59.95571+Dw1_0)      # Source: ELP
    w[1,1]  = (1732559343.73604+Dw1_1)/r2as           # Source: ELP
    w[1,2]  = (        -6.8084 +Dw1_2)/r2as           # Source: DE405
    w[1,3]  =          0.66040e-2/r2as                  # Source: ELP
    w[1,4]  =         -0.31690e-4/r2as                  # Source: ELP
    
    # W2: mean longitude of the lunar perigee:
    w[2,0]  = elp_dms2rad( 83,21,11.67475+Dw2_0)      # Source: ELP
    w[2,1]  = (  14643420.3171 +Dw2_1)/r2as           # Source: DE405
    w[2,2]  = (       -38.2631)/r2as                  # Source: DE405
    w[2,3]  =         -0.45047e-1/r2as                  # Source: ELP
    w[2,4]  =          0.21301e-3/r2as                  # Source: ELP
    
    # W3: mean longitude of the lunar ascending node:
    w[3,0]  = elp_dms2rad(125, 2,40.39816+Dw3_0)      # Source: ELP
    w[3,1]  = (  -6967919.5383 +Dw3_1)/r2as           # Source: DE405
    w[3,2]  =          6.3590/r2as                    # Source: DE405
    w[3,3]  =          0.76250e-2/r2as                  # Source: ELP
    w[3,4]  =         -0.35860e-4/r2as                  # Source: ELP
    
    # Earth-Moon (EMB) elements:
    # Te: mean longitude of EMB:
    eart = np.zeros(5)
    eart[0] = elp_dms2rad(100,27,59.13885+Deart_0)    # Source: VSOP2000
    eart[1] = (129597742.29300 +Deart_1)/r2as         # Source: VSOP2000
    eart[2] =         -0.020200/r2as                  # Source: ELP
    eart[3] =          0.90000e-5/r2as                  # Source: ELP
    eart[4] =          0.15000e-6/r2as                  # Source: ELP
    
    # Pip: mean longitude of the perihelion of EMB:
    peri = np.zeros(5)
    peri[0] = elp_dms2rad(102,56,14.45766+Dperi)      # Source: VSOP2000
    peri[1] =       1161.24342/r2as                   # Source: VSOP2000
    peri[2] =          0.529265/r2as                  # Source: VSOP2000
    peri[3] =         -0.11814e-3/r2as                  # Source: VSOP2000
    peri[4] =          0.11379e-4/r2as                  # Source: VSOP2000
    
    # Corrections to the secular terms of Moon angles.  This gives a better (long-term?) fit
    #   to DE 406.  See ELPdoc, Table 6/paper, Table 4, line 2:
    if(mode==1):  # DE 405 / DE 406
       w[1,3] -= 0.00018865/r2as
       w[1,4] -= 0.00001024/r2as
       
       w[2,2] -= 0.00470602/r2as
       w[2,3] -= 0.00025213/r2as
       
       w[3,2] -= 0.00261070/r2as
       w[3,3] -= 0.00010712/r2as

    
    # Corrections to the mean motions of the Moon angles W2 and W3, infered from the modifications of the constants:
    x2     =   w[2,1] / w[1,1]
    x3     =   w[3,1] / w[1,1]
    y2     =   am*bp[1,1] + xa*bp[5,1]
    y3     =   am*bp[1,2] + xa*bp[5,2]
    
    d21    =   x2 - y2
    d22    =   w[1,1] * bp[2,1]
    d23    =   w[1,1] * bp[3,1]
    d24    =   w[1,1] * bp[4,1]
    d25    =   y2/am
    
    d31    =   x3 - y3
    d32    =   w[1,1] * bp[2,2]
    d33    =   w[1,1] * bp[3,2]
    d34    =   w[1,1] * bp[4,2]
    d35    =   y3/am
    
    Cw2_1  =  d21*Dw1_1+d25*Deart_1+d22*Dgam+d23*De+d24*Dep
    Cw3_1  =  d31*Dw1_1+d35*Deart_1+d32*Dgam+d33*De+d34*Dep
    
    w[2,1] +=  Cw2_1/r2as
    w[3,1] +=  Cw3_1/r2as
    
    # Arguments of Delaunay:
    dela = np.zeros((5,5))  # Need [1:4,0:4] - del is a keyword in Python
    for iD in range(5):     # do iD=0,4
       dela[1,iD] = w[1,iD]  - eart[iD]                 # D   =  W1 - Te + 180 degrees
       dela[2,iD] = w[1,iD]  - w[3,iD]                  # F   =  W1 - W3
       dela[3,iD] = w[1,iD]  - w[2,iD]                  # l   =  W1 - W2   mean anomaly of the Moon
       dela[4,iD] = eart[iD] - peri[iD]                 # l'  =  Te - Pip  mean anomaly of EMB
    
    dela[1,0] = dela[1,0] + pi
    
    # Planetary arguments: mean longitudes for J2000 (from VSOP2000):
    p = np.zeros((9,5))  # p(8,0:4)
    p[1,0] = elp_dms2rad(252, 15,  3.216919)         # Mercury
    p[2,0] = elp_dms2rad(181, 58, 44.758419)         # Venus
    p[3,0] = elp_dms2rad(100, 27, 59.138850)         # EMB (eart(0))
    p[4,0] = elp_dms2rad(355, 26,  3.642778)         # Mars
    p[5,0] = elp_dms2rad( 34, 21,  5.379392)         # Jupiter
    p[6,0] = elp_dms2rad( 50,  4, 38.902495)         # Saturn
    p[7,0] = elp_dms2rad(314,  3,  4.354234)         # Uranus
    p[8,0] = elp_dms2rad(304, 20, 56.808371)         # Neptune
    
    # Planetary arguments: mean motions (from VSOP2000):
    p[1,1] = 538101628.66888/r2as                    # Mercury
    p[2,1] = 210664136.45777/r2as                    # Venus
    p[3,1] = 129597742.29300/r2as                    # EMB (eart(1))
    p[4,1] =  68905077.65936/r2as                    # Mars
    p[5,1] =  10925660.57335/r2as                    # Jupiter
    p[6,1] =   4399609.33632/r2as                    # Saturn
    p[7,1] =   1542482.57845/r2as                    # Uranus
    p[8,1] =    786547.89700/r2as                    # Neptune
    
    p[1:9,2:5] = 0  # p(1:8,2:4) = 0
    
    
    # Zeta: Mean longitude of the Moon W1 + Rate of precession (pt):
    zeta = np.zeros(5)
    zeta[0] = w[1,0]
    zeta[1] = w[1,1] + (5029.0966+Dprec)/r2as
    zeta[2] = w[1,2]
    zeta[3] = w[1,3]
    zeta[4] = w[1,4]
    
    # Corrections to the parameters: Nu, E, Gamma, n' et e' (Source: ELP):
    delnu  = (+0.55604+Dw1_1)/r2as/w[1,1]                 # Correction to the mean motion of the Moon
    dele   = (+0.01789+De)/r2as                           # Correction to the half coefficient of sin(l) in longitude
    delg   = (-0.08066+Dgam)/r2as                         # Correction to the half coefficient of sin(F) in latitude
    delnp  = (-0.06424+Deart_1)/r2as/w[1,1]               # Correction to the mean motion of EMB
    delep  = (-0.12879+Dep)/r2as                          # Correction to the eccentricity of EMB
    
    # Precession of the longitude of the ascending node of the mean ecliptic of date on fixed ecliptic J2000 (Laskar, 1986):
    # P: sine coefficients:
    p1 =  0.10180391e-4
    p2 =  0.47020439e-6
    p3 = -0.5417367e-9
    p4 = -0.2507948e-11
    p5 =  0.463486e-14
    
    # Q: cosine coefficients:
    q1 = -0.113469002e-3
    q2 =  0.12372674e-6
    q3 =  0.1265417e-8
    q4 = -0.1371808e-11
    q5 = -0.320334e-14
    
    return w,eart,peri, zeta,dela,  p,delnu,dele,delg,delnp,delep,dtasm,am,  p1,p2,p3,p4,p5, q1,q2,q3,q4,q5
#***************************************************************************************************
  
  

#***************************************************************************************************
#> \brief  Read the six data files containing the ELP/MPP02 series
##
## \retval ierr      File rrror index: ierr=0: no error, ierr=1: file error
##
## \note
## - module elp_mpp02_constants:  Set of the constants of ELP/MPP02 solution (input)
## - module elp_mpp02_series:  Series of the ELP/MPP02 solution (output)
  
def elp_mpp02_read_files(w,eart,peri, zeta,dela,  p,delnu,dele,delg,delnp,delep,dtasm,am,  p1,p2,p3,p4,p5, q1,q2,q3,q4,q5):
    print("Read files:")
    #use TheSky_elp_mpp02_series, only: cmpb,fmpb,nmpb,   cper,fper,nper
    #use TheSky_elp_mpp02_constants, only: zeta,dela,   p,delnu,dele,delg,delnp,delep,dtasm,am
    
    #integer :: ilu(4),ifi(16)

    # Global variables?:
    nmpb = np.zeros((4,4))  # nmpb(3,3)
    cmpb = np.zeros(2646)   # cmpb(max1), max1 = 2645
    fmpb = np.zeros((5,2646))  # fmpb(0:4,max1)
    nper = np.zeros((4,4,4))  # nper(3,0:3,3)
    cper = np.zeros(33257)    #  cper(max2), max2=33256
    fper = np.zeros((5,33257))  # fper(0:4,max2)
    
    
    # Read the Main Problem series:
    ir=0
    
    ilu = np.zeros(5)  # int(!) ilu(4)
    a = 0.
    b = np.zeros(6)  # double b(5)
    #ierr=1
    nerr=0
    
    # Name of the (here single) ELPMPP02 file:
    fileName = 'data/elp_mpp02.dat'
    inFile = open(fileName,'r')
    #inquire(file=trim(fileName), exist=fexist)
    #if(not fexist) file_open_error_quit(trim(fileName), 1, 1)  # 1: input file
    
    import fortranformat as ff
    formatMainHeader = ff.FortranRecordReader('(25x,I10)')              # Block header format
    formatMainBody   = ff.FortranRecordReader('(4I3,2x,F13.5,5F12.2)')  # Block body format
    
    
    for iFile in range(1,4):  # do iFile=1,3  # These used to be three files
        line = inFile.readline()
        nmpb[iFile,1] = formatMainHeader.read(line)[0]
        #print(nmpb[iFile,1])
        #if(nerr!=0): return 3
        
        nmpb[iFile,2] = ir+1
        nmpb[iFile,3] = nmpb[iFile,1] + nmpb[iFile,2] - 1
        
        for iLine in range(1, int(round(nmpb[iFile,1]))+1):  # do iLine=1,nmpb(iFile,1)
            line = inFile.readline()
            ilu[1],ilu[2],ilu[3],ilu[4], a, b[1],b[2],b[3],b[4],b[5] = formatMainBody.read(line)
            #print(iFile,iLine, "   ", ilu, a, b)
            #if(nerr!=0): return 4
            
            ir=ir+1
            tgv = b[1] + dtasm*b[5]
            if(iFile==3):  a -= 2*a*delnu/3
            cmpb[ir] = a + tgv*(delnp-am*delnu) + b[2]*delg + b[3]*dele + b[4]*delep
            
            for k in range(5):  # do k=0,4
                fmpb[k,ir] = 0
                for i in range(1,5):  # do i=1,4
                    fmpb[k,ir] += ilu[i] * dela[i,k]
                    
                    if(iFile==3): fmpb[0,ir] += pio2

    
    # Read the Perturbations series:
    ir=0
    ipt = 0
    icount = 0
    s = 0.0
    c = 0.0
    ifi = np.zeros(17)  # int ifi(16)
    nper = np.zeros((4,4,4))  # nper(3,0:3,3)
    
    formatPertHeader = ff.FortranRecordReader('(25x,2I10)')         # Perturbation header format
    formatPertBody   = ff.FortranRecordReader('(I5,2D20.13,16I3)')  # Perturbation body format
    
    for iFile in range(1,4):  # do iFile=1,3  # These used to be three files
        for it in range(4):   # do it=0,3
            #read(ip,'(25x,2I10)', iostat=nerr,end=100) nper(iFile,it,1),ipt
            #if(nerr!=0): return 6
            line = inFile.readline()
            #print(line)
            nper[iFile,it,1],ipt = formatPertHeader.read(line)
            #print(iFile,iLine, "   ", nper[iFile,it,1],ipt)
            
            nper[iFile,it,2] = ir+1
            nper[iFile,it,3] = nper[iFile,it,1] + nper[iFile,it,2] - 1
            if(nper[iFile,it,1]==0): continue  # cycle
            
            nLines = int(round(nper[iFile,it,1]))
            for iLine in range(1, nLines+1):  # do iLine=1,nper(iFile,it,1)
                #read(ip,'(I5,2D20.13,16I3)', iostat=nerr,end=100) icount,s,c,ifi
                #if(nerr!=0): return 7
                line = inFile.readline()
                icount,s,c,ifi[1],ifi[2],ifi[3],ifi[4],ifi[5],ifi[6],ifi[7],ifi[8],ifi[9],ifi[10],ifi[11],ifi[12],ifi[13],ifi[14],ifi[15],ifi[16] = formatPertBody.read(line)
                #if(iLine==1 or iLine==nLines): print(iFile,iLine, "     ", icount,s,c,ifi)
                
                ir = ir+1
                cper[ir] = m.sqrt(c**2+s**2)
                pha = m.atan2(c,s)
                if(pha<0): pha = pha+pi2
                
                for k in range(5):  # do k=0,4
                    fper[k,ir] = 0
                    if(k==0): fper[k,ir] = pha
                    for i in range(1,5):  # do i=1,4
                        fper[k,ir] += ifi[i] * dela[i,k]
                        
                    for i in range(5,13):  # do i=5,12
                        fper[k,ir] += ifi[i] * p[i-4,k]
                    
                    fper[k,ir] += ifi[13] * zeta[k]
                
    inFile.close()
    
    # Exit:
    ierr=0
    return ierr
    
    # End of file error:
# 100 continue
    ierr=9
    
    icount = icount  # Suppress 'variable set but not used' compiler warnings
    ipt = ipt        # Suppress 'variable set but not used' compiler warnings
    
    return ierr
#***************************************************************************************************
  
  
#***************************************************************************************************
#> \brief Function for the conversion: sexagesimal degrees -> radians
def elp_dms2rad(deg,min,sec):
    return (deg+min/60+sec/3600)*d2r
#***************************************************************************************************


#*********************************************************************************************************************************
#> \brief  Compute the spherical lunar coordinates using the ELP2000/MPP02 lunar theory in the dynamical mean ecliptic and
##           equinox of J2000.
##
## \param jd    Julian day to compute Moon position for
## \param mode  Index of the corrections to the constants: 0-Fit to LLR observations, 1-Fit to DE405 1950-2060 (historical)
##
## \retval  lon  Ecliptic longitude (rad)
## \retval  lat  Ecliptic latitude (rad)
## \retval  rad  Distance (AU)

def elp_mpp02_lbr(jd, mode):
    print("Compute lbr:")
    #real(double) :: xyz(3),vxyz(3)
    
    xyz,vxyz, ierr = elp_mpp02_xyz(jd, mode)
    
    # Compute ecliptic l,b,r:
    rad = m.sqrt(sum(xyz**2))
    lon = m.atan2(xyz[2], xyz[1])
    lat = m.asin(xyz[3]/rad)
    
    # precess_ecl(jd2000,jd, lon,lat)
    
    rad = rad/1.49597870700e8  # km -> AU
    
    #write(*,'(/, F20.5, 2F14.7,F14.5)') jd, rev(lon)*r2d, lat*r2d, rad
    return lon,lat,rad
#*********************************************************************************************************************************


#***************************************************************************************************
#> \brief  Compute the rectangular lunar coordinates using the ELP/MPP02 lunar theory in the dynamical mean ecliptic and equinox of J2000.
##
## \param jd    Julian day to compute Moon position for
## \param mode  Index of the corrections to the constants: 0-Fit to LLR observations, 1-Fit to DE405 1950-2060 (historical)
## 
## \retval xyz   Geocentric rectangular coordinates:
##               - xyz(1) : Position X (km)
##               - xyz(2) : Position Y (km)
##               - xyz(3) : Position Z (km)
## \retval vxyz  Geocentric rectangular velocities:
##               - vxyz(1) : Velocity X' (km/day)
##               - vxyz(2) : Velocity Y' (km/day)
##               - vxyz(3) : Velocity Z' (km/day)
## \retval ierr  File error index - ierr=0: no error, ierr=1: file error
##
## \note
##  - The subroutine elp_mpp02() uses two modules:
##    - elp_mpp02_constants:  Constants of the solution ELP/MPP02 (input),
##    - elp_mpp02_series:     Series of the solution ELP/MPP02 (input).
##
##  - The nominal values of some constants have to be corrected.  There are two sets of corrections, which can be selected
##    using the parameter 'mode' (used in elp_mpp02_initialise()).
##    - mode=0, the constants are fitted to LLR observations provided from 1970 to 2001; it is the default value;
##    - mode=1, the constants are fitted to DE405 ephemeris over one century (1950-2060); the lunar angles W1, W2, W3
##              receive also additive corrections to the secular coefficients ('historical mode').
##    When the mode is changed, the data must be reinitialised and the data file reread.
##
##  - Solutions (discussed) in the paper:
##    - ELP (original):
##      - ELP2000-82: using VSOP82 (1983)
##      - ELP2000-85: new mean lunar arguments, higher truncation level, longer time range (1988)
##      - ELP2000-82B, here called "ELP": ELP2000-82, using mean lunar arguments from ELP2000-85 (19??)
##    - ELP/MPP01:  using latest planetary perturbations from MPP01 and VSOP2000, but simpler than MPP01
##    - ELP/MPP02:  ELP/MPP01, but for some arguments back to ELP + different selection of perturbations + lower truncation.  Good fit with DE 405 in [1950,2060]
##    - ELP/MPP02*: improved secular arguments, better long-term comparison to DE 405/406 [-3000,2500]
##    - ELP/MPP02(LLR): ELP/MPP02(*?), optimised for lunar ranging since 1970
##    - ELPa: ELP + few Poisson terms (tested in the current study only?)
##    - ELPa*: ELPa + better secular arguments (as in ELP/MPP02*)
##  - It is not entirely clear which version is given below, but we can hope it is ELP/MPP02*.  However, the subroutine
##      elp82b_lbr() above is known to underperform (by a factor of 10) in accuracy.
  
def elp_mpp02_xyz(jd, mode):
    print("Compute xyz:")
    #use TheSky_elp_mpp02_series, only: cmpb,fmpb,nmpb,   cper,fper,nper
    #use TheSky_elp_mpp02_constants, only: w, p1,p2,p3,p4,p5, q1,q2,q3,q4,q5
    
    # Constants:
    a405=384747.9613701725
    aelp=384747.980674318
    sc=36525  # Moon mean distance for DE405 und ELP; Julian century in days
    
    # Global variables?
    nmpb = np.zeros((4,4))  # nmpb(3,3)
    cmpb = np.zeros(2646)   # cmpb(max1), max1 = 2645
    fmpb = np.zeros((5,2646))  # fmpb(0:4,max1)
    nper = np.zeros((4,4,4))  # nper(3,0:3,3)
    cper = np.zeros(33257)    #  cper(max2), max2=33256
    fper = np.zeros((5,33257))  # fper(0:4,max2)
    w = np.zeros((4,5))  # w(3,0:4)
    p1=0; p2=0; p3=0; p4=0; p5=0;  q1=0; q2=0; q3=0; q4=0; q5=0;
    #real(double) :: rjd, t(-1:4),v(6)
    
    # Initialise data and read files if needed:
    ierr = elp_mpp02_initialise_and_read_files(mode)
    if(ierr!=0): sys.exit('Could not read ELP-MPP02 files')
    
    
    # Initialization of time powers:
    rjd  = jd - jd2000  # Reduced JD - JD since 2000
    t = np.zeros(5)  # t(-1:4) - -1 used to prevent compiler warnings?
    t[0] = 1
    t[1] = rjd/sc       # t: time since 2000 in Julian centuries
    t[2] = t[1]**2      # t^2
    t[3] = t[2]*t[1]    # t^3
    t[4] = t[2]**2      # t^4
    
    # Evaluation of the series: substitution of time in the series
    v = np.zeros(7)  # v(6)
    for iVar in range(1,4):  # do iVar=1,3  # iVar=1,2,3: Longitude, Latitude, Distance
        v[iVar] = 0
        v[iVar+3] = 0
        
        # Main Problem series:
        for iLine in range(int(round(nmpb[iVar,2])), int(round(nmpb[iVar,3]))+1):  # do iLine=nmpb(iVar,2),nmpb(iVar,3)
            x = cmpb[iLine]
            y = fmpb[0,iLine]
            yp = 0
            
            for k in range(1,5):  # do k=1,4
                y  = y  +   fmpb[k,iLine] * t[k]
                yp = yp + k*fmpb[k,iLine] * t[k-1]
          
            v[iVar]   = v[iVar]    +  x * m.sin(y)
            v[iVar+3] = v[iVar+3]  +  x *yp * m.cos(y)
            
            
        # Perturbations series:
        for it in range(4):  # do it=0,3
            for iLine in range(int(round(nper[iVar,it,2])), int(round(nper[iVar,it,3]))+1):  # do iLine=nper(iVar,it,2),nper(iVar,it,3)
                x = cper[iLine]
                y = fper[0,iLine]
                xp = 0
                yp = 0
                if(it!=0): xp = it * x * t[it-1]
                
                for k in range(1,5):  # do k=1,4
                    y +=  fper[k,iLine] * t[k]
                    yp += k * fper[k,iLine] * t[k-1]
             
            v[iVar] += x * t[it] * m.sin(y)
            v[iVar+3] += xp * m.sin(y) + x * t[it] * yp * m.cos(y)
    
    # Compute the spherical coordinates for the mean inertial ecliptic and equinox of date:
    v[1]   = v[1]/r2as + w[1,0] + w[1,1]*t[1] + w[1,2]*t[2] + w[1,3]*t[3] + w[1,4]*t[4]  # Longitude + mean longitude (rad)
    v[2]   = v[2]/r2as                                                                   # Latitude (rad)
    v[3]   = v[3] * a405 / aelp                                                          # Distance (km)
    
    #lambda = v[1]
    #beta = v[2]
    #rad = v[3]
    #lambda = lambda + [5029.0966*t[1] + 1.1120*t[2] + 0.000077*t[3] - 0.00002353*t[4]  -  0.29965*t(1]) * as2r  # Precession from J2000 to EoD(?), but only in longitude#
    #write(*,'(2F14.7,F14.5)') rev(lambda)*r2d, beta*r2d, rad
    
    
    # Compute the rectangular coordinates (for the EoD?):
    clamb  = m.cos(v[1])
    slamb  = m.sin(v[1])
    cbeta  = m.cos(v[2])
    sbeta  = m.sin(v[2])
    cw     = v[3]*cbeta
    sw     = v[3]*sbeta
    
    x1     = cw*clamb
    x2     = cw*slamb
    x3     = sw
    
    # Is this simply precession in rectangular coordinates from EoD to J2000?
    pw     = (p1 + p2*t[1] + p3*t[2] + p4*t[3] + p5*t[4]) * t[1]
    qw     = (q1 + q2*t[1] + q3*t[2] + q4*t[3] + q5*t[4]) * t[1]
    
    ra     = 2*m.sqrt(1 - pw**2 - qw**2)
    pwqw   = 2*pw*qw
    pw2    = 1 - 2*pw**2
    qw2    = 1 - 2*qw**2
    pwra   = pw*ra
    qwra   = qw*ra
    
    xyz = np.zeros(4)  # xyz(3)
    xyz[1] =  pw2*x1  + pwqw*x2 + pwra*x3
    xyz[2] =  pwqw*x1 + qw2*x2  - qwra*x3
    xyz[3] = -pwra*x1 + qwra*x2 + (pw2+qw2-1)*x3
    
    #xyz[1] = x1
    #xyz[2] = x2
    #xyz[3] = x3
    
    
    # Compute the rectangular velocities for the equinox J2000:
    v[4]   = v[4]/r2as + w[1,1] + 2*w[1,2]*t[1] + 3*w[1,3]*t[2] + 4*w[1,4]*t[3]
    v[5]   = v[5]/r2as
    
    xp1    = (v[6]*cbeta - v[5]*sw)*clamb - v[4]*x2
    xp2    = (v[6]*cbeta - v[5]*sw)*slamb + v[4]*x1
    xp3    = v[6]*sbeta  + v[5]*cw
    
    ppw    = p1 + (2*p2 + 3*p3*t[1] + 4*p4*t[2] + 5*p5*t[3]) * t[1]
    qpw    = q1 + (2*q2 + 3*q3*t[1] + 4*q4*t[2] + 5*q5*t[3]) * t[1]
    ppw2   = -4*pw*ppw
    qpw2   = -4*qw*qpw
    ppwqpw = 2*(ppw*qw + pw*qpw)
    rap    = (ppw2+qpw2)/ra
    ppwra  = ppw*ra + pw*rap
    qpwra  = qpw*ra + qw*rap
    
    vxyz = np.zeros(4)  # vxyz(3)
    vxyz[1] = (pw2*xp1 + pwqw*xp2 + pwra*xp3  +  ppw2*x1 + ppwqpw*x2 + ppwra*x3) / sc
    vxyz[2] = (pwqw*xp1 + qw2*xp2 - qwra*xp3  +  ppwqpw*x1 + qpw2*x2 - qpwra*x3) / sc
    vxyz[3] = (-pwra*xp1 + qwra*xp2 + (pw2+qw2-1)*xp3  -  ppwra*x1 + qpwra*x2 + (ppw2+qpw2)*x3) / sc
    
    return xyz,vxyz, ierr
#***************************************************************************************************
  
  
  

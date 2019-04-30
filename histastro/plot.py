"""HistAstro plot-related functions"""

def mag2size(Mlim, mag, scale=1):
    """Convert magnitudes (mag) to disc area ('size') for plotting with pyplot.scatter(), given a magnitude limit Mlim"""
    
    size = scale * 30*(0.5 + (Mlim-mag)/3.0)**2     # Disc area as a function of magnitude
    return size


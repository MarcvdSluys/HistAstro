"""Plot-related functions for HistAstro."""

def mag2size(Mlim, mag, scale=1):
    """
    Convert magnitudes (mag) to disc area ('size') for plotting with pyplot.scatter(), given a magnitude limit Mlim.
    
    Args:
      Mlim (double):   Magnitude limit.
      mag (double):    Magnitude (Numpy array).
      scale (double):  Scale factor (default value = 1).
    
    Returns:
      double:  Area sizes for the stellar discs (Numpy array)
    
    """
    
    size = scale * 30*(0.5 + (Mlim-mag)/3.0)**2     # Disc area as a function of magnitude
    
    return size


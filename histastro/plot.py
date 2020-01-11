#  Copyright (c) 2019-2020  Marc van der Sluys - marc.vandersluys.nl
#   
#  This file is part of the HistAstro Python package,
#  see: http://astro.ru.nl/~sluys/HistAstro/
#   
#  This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#  
#  This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
#  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License along with this code.  If not, see
#  <http://www.gnu.org/licenses/>.


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


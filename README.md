# HistAstro #

A Python package for historical-astronomy calculations of Sun, Moon and planets.


## Installation ##

This package can be installed using `pip install histastro`.  This should automatically install the
dependency packages `numpy` and `fortranformat` if they haven't been installed already.


### Data files ###

You will need to manually download the eight
data files `VSOP87D.xxx` from ftp://ftp.imcce.fr/pub/ephem/planets/vsop87/, one for each of the planets, and
save them in a directory of your choice (e.g. the subdirectory `data/` of your working directory).  To compute
the Moon position, the file `moonposMeeus.csv` is needed, which can be downloaded from this page:
http://astro.ru.nl/~sluys/HistAstro/  



## Author and licence ##

* Author: Marc van der Sluys
* Contact: http://astro.ru.nl/~sluys/
* Website: [Github](https://github.com/MarcvdSluys/HistAstro), [PyPI](https://pypi.org/project/histastro/)
* Licence: [GPLv3+](https://www.gnu.org/licenses/gpl.html)


## References ##

* [Chapront & Francou (2003)](https://ui.adsabs.harvard.edu/abs/2003A%26A...404..735C/abstract)
* [FTP data files](ftp://cyrano-se.obspm.fr/pub/2_lunar_solutions/2_elpmpp02) &mdash; in case [FTP urls don't work in Markdown](https://github.com/gollum/gollum/issues/759): ftp://cyrano-se.obspm.fr/pub/2_lunar_solutions/2_elpmpp02
* This Python code is adapted from the Fortran implementation in [libTheSky](http://libthesky.sourceforge.net/)

# HistAstro #

A Python package for historical-astronomy calculations of Sun, Moon and planets.


## Installation ##

This package can be installed using `pip install histastro`.  This should automatically install the
dependency packages `numpy` and `fortranformat` if they haven't been installed already.  If you are installing
by hand, you have to ensure these packages are installed as well.


### Data files ###

You will need to manually download the eight data files
[`VSOP87D.xxx`](http://cdsarc.u-strasbg.fr/viz-bin/Cat?cat=VI/81) (click on FTP), one for each of the planets,
and save them in a directory of your choice (e.g. the subdirectory `data/` of your working directory).  To
compute the Moon position, the file `moonposMeeus.csv` is needed, which can be downloaded from this page:
http://astro.ru.nl/~sluys/HistAstro/


## HistAstro pages ##

* [Radboud University](http://astro.ru.nl/~sluys/HistAstro/): HistAstro homepage: data files and code documentation
* [Pypi](https://pypi.org/project/histastro/): HistAstro Python package
* [GitHub](https://github.com/MarcvdSluys/HistAstro): HistAstro source code


## Author and licence ##

* Author: Marc van der Sluys
* Contact: http://astro.ru.nl/~sluys/
* Licence: [GPLv3+](https://www.gnu.org/licenses/gpl.html)


## References ##

* Meeus, [Astronomical algorithms](https://www.willbell.com/math/MC1.HTM), 2nd Ed.
* VSOP87 planet theory: [Bretagnon & Francou, 1988A+A...202..309B](https://ui.adsabs.harvard.edu/abs/1988A%26A...202..309B/)
* [VSOP87 data files](http://cdsarc.u-strasbg.fr/viz-bin/Cat?cat=VI/81)
* This Python code is adapted from the Fortran implementation in [libTheSky](http://libthesky.sourceforge.net/)
* [Celestial mechanics in a nutshell (CMiaNS)](https://cmians.sourceforge.io/)

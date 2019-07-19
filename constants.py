#!/usr/bin/env python

__author__ = "Pablo Marchant"
__credits__ = ["Pablo Marchant"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Pablo Marchant"
__email__ = "pamarca@gmail.com"

"""Various constants, all in CGS
these are the MESA defined values
astro variables are from Bahcall et al, ApJ 618 (2005) 1049-1056
"""

clight = 2.99792458e10
cgrav = 6.67428e-8
Msun = 1.9892e33
Rsun = 6.9598e10
Lsun = 3.8418e33
au = 1.495978921e13
secyear = 24*3600*365.25

def merger_time(a0,m1,m2):
	Beta = (64/5) * cgrav**3 * m1 * m2 * (m1 + m2)*Msun**3/ clight**5
	Beta = Beta/(3.1536*1e16)
	merger_time = (a0**4)/(4*Beta)
	return merger_time




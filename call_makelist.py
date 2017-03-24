from numpy import *
import numpy as np
import os
import csv

from astropy import wcs
import astropy.units as u
from astropy import coordinates as coord

from astroquery.sdss import SDSS
from astroquery.vizier import Vizier

import matplotlib.pyplot as plt

from input_skymaker import makelist

dec_list = np.linspace(0,65,25)
ra_list = np.linspace (160,200,25)

for i in range(1,26):
    fname="GOTO_01_20170323_{0:05}".format(i)
    makelist(ra_list[i-1],dec_list[i-1],fname)

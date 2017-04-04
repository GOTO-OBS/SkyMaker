import pyfits
import os
from datetime import datetime
from astropy.io.fits import getheader
import numpy as np
from astropy.io import fits

date=datetime.strftime(datetime.now(), '%Y%m%d')

for file in os.listdir(date):
    if file.endswith(".fits"):
        data, header = fits.getdata(os.path.join(date,file), header=True)
        #print header
        header.append('IMGTYPE','INSTRUME','RUN') 
        header['IMGTYPE']= 'OBJECT'
        print header['IMGTYPE']
        print header
        fits.writeto(str(os.path.join(date,file)), data, header, clobber=True)

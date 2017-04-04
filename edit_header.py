import pyfits
import os
from datetime import datetime
from astropy.io.fits import getheader
import numpy as np
from astropy.io import fits

def edit_header(date,ccd,visit,\
                ccd_ra,ccd_dec,\
                mount_ra, mount_dec,\
                fwhm):

    #Construct filename and read file:
    fname="GOTO_"+ccd+"_"+date+"_"+visit+".fits"                
    data, header = fits.getdata(os.path.join(date,fname), header=True)

    #Add the necessary information:
    header.append('IMGTYPE','INSTRUME','RUN') 
    header['IMGTYPE']= 'OBJECT'

    #You'll need to add:
    #WCS, FWHM, mount (RA,Dec), CCD (RA,Dec),
    #Date, CCD, visit number.
    #Remember: (ccd_ra,cdd_dec) refers to the central\
    #pixel of the CCD.

    #Save file with new header information:
    fits.writeto(str(os.path.join(date,fname)), data, header, clobber=True)

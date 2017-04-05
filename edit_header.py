import pyfits
import os
from datetime import datetime
from astropy.io.fits import getheader
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import Angle


def edit_header(date,ccd,visit,\
                ccd_ra,ccd_dec,\
                mount_ra, mount_dec,\
                fwhm):

    #Construct filename and read file:
    fname="GOTO_"+ccd+"_"+date+"_"+visit+".fits"                
    data, header = fits.getdata(os.path.join(date,fname), header=True)

    #Add the necessary information:
    header.append('IMGTYPE','INSTRUME','RUN','FILTER')
    header.append('TEL-RA','TEL-DEC','CRVAL1','CRVAL2') 
    header.append('EPOCH','EQUINOX','CTYPE1','CTYPE2')

    header['DATE-OBS'] = str(date)
    header['IMGTYPE'] = 'OBJECT'

    header['INSTRUME'] =str(ccd)
    header['RUN'] = str(visit)
    header['FILTER'] = 'Clear'

    header['TEL-RA'] = Angle(mount_ra, u.deg).to_string(unit=u.hour, sep=':')
    header['TEL-DEC'] = Angle(mount_dec, u.deg).to_string(unit=u.degree, sep=':')
    header['CRVAL1'] = ccd_ra
    header['CRVAL2'] = ccd_dec

    header['EPOCH'] = 2000
    header['EQUINOX'] =2000
    header['CTYPE1'] = 'RA---TAN'
    header['CTYPE2'] = 'DEC---TAN'
    

    #You'll need to add:
    #WCS, FWHM, mount (RA,Dec), CCD (RA,Dec),
    #Date, CCD, visit number.
    #Remember: (ccd_ra,cdd_dec) refers to the central\
    #pixel of the CCD.

    #Save file with new header information:
    fits.writeto(str(os.path.join(date,fname)), data, header, clobber=True)

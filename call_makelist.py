import numpy as np
from astropy import wcs

#from input_skymaker import makelist

def mount_pointing(x, y, \
                   width=4., height=4., \
                   w_olap=0.1, h_olap=0.1):

    #Calculate declination of the pointing:
    d_dec = y*(height-h_olap)
    dec = d_dec

    #Use the Haversine formula to calculate
    #where the RA pointing is.
    d_ra = x*(width - w_olap)
    theta = np.pi*d_ra/180.
    ra = 2.*np.arcsin(np.sin(theta/2.)/\
                      np.cos(np.pi*dec/180.))
    ra = 180.*ra/np.pi
    
    return ra, dec

def ccd_pointing(ra, dec, \
                 width=2, height=2,
                 w_olap=0.1, h_olap=0.1):

    #Calculate the declination of the ccd pointing:
    ccd_dec = dec + np.array([-1,-1,1,1])*(height-h_olap)/2.
    
    #Use the Haversine formula to calculate
    #the RA of the centre of each ccd.
    d_ra = np.array([-1,1,-1,1])*(width-w_olap)/2.
    theta = np.pi*d_ra/180.
    ccd_ra = 2.*np.arcsin(np.sin(theta/2.)/\
                          np.cos(np.pi*ccd_dec/180.))  
    ccd_ra = 180.*ccd_ra/np.pi + ra
    
    return ccd_ra, ccd_dec
    
ra, dec = mount_pointing(5,10)
ccd_ra, ccd_dec = ccd_pointing(ra, dec)
for ra in np.ndenumerate(ccd_ra):
    print ra[1], ccd_dec[ra[0]]

#    fname="GOTO_01_20170323_{0:05}".format(i)
#    makelist(ra_list[i-1],dec_list[i-1],fname)

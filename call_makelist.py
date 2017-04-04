import numpy as np
from astropy import wcs

from input_skymaker import makelist

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
    
    sin_theta = np.sin(theta/2.)/\
                       np.cos(np.pi*dec/180.)
    opp = False
    if sin_theta>1:
        sin_theta = 2.-sin_theta
        opp = True
        
    ra = 2.*np.arcsin(sin_theta)
    ra = 180.*ra/np.pi

    if opp:
        ra = 360-ra
    print ra, dec
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

xs = np.arange(40, 45, 1)
ys = np.arange(0, 5, 1)

i = 0
for y in ys:
    for x in xs:
        i = i+1
        visit = "{0:04}".format(i)
        
        mount_ra, mount_dec = mount_pointing(x, y)
        ccd_ras, ccd_decs = ccd_pointing(mount_ra, mount_dec)

        #for j in np.arange(ccd_ras.size):
        #    ccd = "{0:02}".format(j+1)
        #    fname="GOTO_"+ccd+"_20170331_"+visit
        #    makelist(ccd_ras[j], ccd_decs[j], fname)

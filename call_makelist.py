import numpy as np
from astropy import wcs
import os
from datetime import datetime

from input_skymaker import makelist

#This is a trivial change.

def mount_pointing(x, y, \
                   width=4., height=4., \
                   w_olap=0.1, h_olap=0.1):

    #Calculate declination of the pointing:
    d_dec = y*(height-h_olap)
    dec = d_dec

    #Use the Haversine formula to calculate
    #where the RA pointing is.

    #Width of one pointing in radians:
    theta = np.pi*(width - w_olap)/180.
    
    #Delta RA for one pointing:
    d_ra1 = 2.*np.arcsin(np.sin(theta/2.)/\
                         np.cos(np.pi*dec/180.))

    #RA of x:
    ra = 180.*(x*d_ra1)/np.pi

    if ra >= 360:
        print "x out of range. RA > 360."
        print "Maximum x for y={} is:".format(y), np.int(np.floor(2.*np.pi/d_ra1))
    
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


#Date of obs:
date=datetime.strftime(datetime.now(), '%Y%m%d')

#Width and height of each chip in degs:
xpix=8176
ypix=6132
xsi=1.24*xpix/3600.
ysi=1.24*ypix/3600.

#Width of each pointing (10' overlap between CCDs):
xpsi = 2.*xsi-(10./60.)
ypsi = 2.*ysi-(10./60.)

xs = np.array([1])
ys = np.array([1])
i = 0
for y in ys:
    for x in xs:
        i = i+1
        visit = "{0:04}".format(i)
        
        #Get (RA,Dec) of mount pointing:
        #Overlap between pointings is 15':
        mount_ra, mount_dec = mount_pointing(x, y,\
                                             width=xpsi, height=ypsi,\
                                             w_olap=15./60., h_olap=15./60.)
                                             
        #Get (RA,Dec) of CCDs:
        #Overlap between CCDs is 10':
        ccd_ras, ccd_decs = ccd_pointing(mount_ra, mount_dec,\
                                         width=xsi, height=ysi,\
                                         w_olap=10./60., h_olap=10./60.)
        print mount_ra, mount_dec
        for j in np.arange(ccd_ras.size):
            ccd = "{0:02}".format(j+1)
            fname="GOTO_"+ccd+"_"+date+"_"+visit
            makelist(ccd_ras[j],ccd_decs[j],"templist")
            if not os.path.exists(date):
                 os.makedirs(date)
            with open('goto.conf') as infile, open(date+"/"+str(fname)+'.conf', 'w') as outfile:
                for line in infile:
                    line = line.replace('goto.fits', date+"/"+str(fname)+'.fits')
                    line = line.replace('2.96065834', '1')
                    outfile.write(line)

            os.system('./bin/sky templist.list -c '+ date+"/"+str(fname)+'.conf')   

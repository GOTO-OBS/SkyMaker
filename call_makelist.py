import numpy as np
from astropy import wcs
import os
from datetime import datetime
from input_skymaker import makelist
from edit_header import edit_header

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


#Date of obs and generate directory
#for that date:
date=datetime.strftime(datetime.now(), '%Y%m%d')
if not os.path.exists(date):
    os.makedirs(date)
    os.makedirs(date+'/fits')
    os.makedirs(date+'/conf')
    os.makedirs(date+'/list')
    os.makedirs(date+'/reg')
    
#Define the "mean" FWHM for that night:
fwhm_n = 10.**np.random.normal(loc=np.log10(0.9), scale=0.15)

#Width and height of each chip in degs:
xpix=8176
ypix=6132
xsi=1.24*xpix/3600.
ysi=1.24*ypix/3600.

#Width of each pointing (10' overlap between CCDs):
xpsi = 2.*xsi-(10./60.)
ypsi = 2.*ysi-(10./60.)

xs = np.array([32,33,34,35,31,32,33,34,31,32,33,34,30,31,32,33])
ys = np.array([3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6])

#xs = np.array([32])
#ys = np.array([3])

i = 0
for m in np.arange(xs.size):
    i = i+1
    visit = "{0:04}".format(i)
        
        #Get (RA,Dec) of mount pointing:
        #Overlap between pointings is 15':
    mount_ra, mount_dec = mount_pointing(xs[m], ys[m],\
                                         width=xpsi, height=ypsi,\
                                         w_olap=15./60., h_olap=15./60.)
                                             
        #Get (RA,Dec) of CCDs:
        #Overlap between CCDs is 10':
    ccd_ras, ccd_decs = ccd_pointing(mount_ra, mount_dec,\
                                     width=xsi, height=ysi,\
                                     w_olap=10./60., h_olap=10./60.)
        
        #Make a source list for each CCD:
    for j in np.arange(ccd_ras.size):
        ccd = "{0:02}".format(j+1)
        lname="GOTO_"+ccd+"_"+date+"_"+visit
        makelist(ccd_ras[j], ccd_decs[j], ccd, date+"/reg/"+lname, variability=True)
            
        #Three exposures per pointing:
    for h in np.arange(3):
        expo = "{0:02}".format(h+1)
        fwhm_p = 10.**np.random.normal(loc=np.log10(fwhm_n), scale=0.04)

            #Loop over the CCDs:
        for j in np.arange(ccd_ras.size):
            ccd = "{0:02}".format(j+1)
            fname = "GOTO_"+ccd+"_"+date+"_"+visit+"_"+expo

                #Write the necessary data into the .conf file: 
            with open('goto.conf') as infile, \
                    open(date+"/conf/"+fname+'.conf', 'w') as outfile:
                for line in infile:
                        #Filename of output:
                    line = line.replace('goto.fits', \
                                                date+"/fits/"+fname+'.fits')
                        #Seeing FWHM:    
                    line = line.replace('2.96065834', str(fwhm_p))
                    outfile.write(line)

                #Run SkyMaker with generated list and .conf:
            os.system('sky templist'+ccd+'.list -c '+ date+"/conf/"+fname+'.conf')   
            os.rename(date+'/fits/'+fname+'.list',date+'/list/'+fname+'.list')
                
                #Edit header of output .fits file:
                #Here, the "1." is just a placeholder for the FWHM.
            edit_header(date,ccd,visit,expo,\
                        ccd_ras[j],ccd_decs[j],\
                        mount_ra, mount_dec,\
                        fwhm_p)

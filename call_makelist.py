import numpy as np
from astropy import wcs
from astropy.io import fits
import os
import re
from datetime import datetime
from input_skymaker import makelist
from edit_header import edit_header
import glob
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting

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
date="20170803"
if not os.path.exists(date):
    os.makedirs(date)
    os.makedirs(date+'/fits')
    os.makedirs(date+'/conf')
    os.makedirs(date+'/list')
    os.makedirs(date+'/reg')
    os.makedirs(date+'/var')
    os.makedirs(date+'/radec')
    
#Define the "mean" FWHM for that night:
fwhm_n = 10.**np.random.normal(loc=np.log10(1.55), scale=0.12, size=1)
fwhm_n = np.clip(fwhm_n,0.8,10.)
#plt.hist(fwhm_n, 100)
#plt.show()
#quit()

#Width and height of each chip in degs:
xpix=8176
ypix=6132
xsi=1.24*xpix/3600.
ysi=1.24*ypix/3600.

#Width of each pointing (10' overlap between CCDs):
xpsi = 2.*xsi-(10./60.)
ypsi = 2.*ysi-(10./60.)


#Read in the positions:
ra = np.array([])
dec = np.array([])
files = glob.glob('pointings/15420.txt')
for filen in files:
    with open(filen,'r') as fin:
        for line in fin:
            if 'OBS' in line:
                m = re.search('ra: (.+?);', line)
                ra = np.append(ra, float(m.group(1)))
                m = re.search('dec: (.+?);', line)
                dec = np.append(dec, float(m.group(1)))

ra = np.array(ra[0:1])
dec = np.array(dec[0:1])

#Read in the calibration files:
flat = np.array([fits.getdata('calibs/flat/UT1_FLAT_2017-07-26.fits', 1),\
                 fits.getdata('calibs/flat/UT2_FLAT_2017-07-26.fits', 1),\
                 fits.getdata('calibs/flat/UT2_FLAT_2017-07-26.fits', 1),\
                 fits.getdata('calibs/flat/UT4_FLAT_2017-07-26.fits', 1)])
bias = np.array([fits.getdata('calibs/bias/UT1_BIAS_2017-07-26.fits', 1),\
                 fits.getdata('calibs/bias/UT2_BIAS_2017-07-26.fits', 1),\
                 fits.getdata('calibs/bias/UT2_BIAS_2017-07-26.fits', 1),\
                 fits.getdata('calibs/bias/UT4_BIAS_2017-07-26.fits', 1)])
dark = np.array([fits.getdata('calibs/dark/UT1_DARK_2017-07-26.fits', 1),\
                 fits.getdata('calibs/dark/UT2_DARK_2017-07-26.fits', 1),\
                 fits.getdata('calibs/dark/UT2_DARK_2017-07-26.fits', 1),\
                 fits.getdata('calibs/dark/UT4_DARK_2017-07-26.fits', 1)])
scan = np.array([fits.getdata('calibs/scan/BIAS-UT1-r0002212.fts', 0),\
                 fits.getdata('calibs/scan/BIAS-UT2-r0002212.fts', 0),\
                 fits.getdata('calibs/scan/BIAS-UT2-r0002212.fts', 0),\
                 fits.getdata('calibs/scan/BIAS-UT4-r0002212.fts', 0)])

y = np.mean(scan[1,1:,0:14],axis=1)
x = np.arange(y.size)+1
poly = models.Polynomial1D(degree=4)
fitfun = fitting.LevMarLSQFitter()
fit = fitfun(poly, x, y)
#Fit with a 4th order polynomial:

plt.plot(x, y)
plt.plot(fit(x))
plt.show()
quit()

i = 0
for m in np.arange(ra.size):
    i = i+1
    visit = "{0:04}".format(i)
        
    #Get (RA,Dec) of mount pointing:
    mount_ra, mount_dec = ra[m], dec[m]
                                             
    #Get (RA,Dec) of CCDs:
    #Overlap between CCDs is 10':
    ccd_ras, ccd_decs = ccd_pointing(mount_ra, mount_dec,\
                                     width=xsi, height=ysi,\
                                     w_olap=10./60., h_olap=10./60.)
       
    #Make a source list for each CCD:
    for j in np.arange(ccd_ras.size):
        print "CCD: ", j
        ccd = "{0:02}".format(j+1)
        lname="GOTO_"+ccd+"_"+date+"_"+visit
        makelist(ccd_ras[j], ccd_decs[j], ccd, date, lname, variable=0.03, transients=50)
            
    #Three exposures per pointing:
    for h in np.arange(3):
        expo = "{0:02}".format(h+1)
        fwhm_p = 10.**np.random.normal(loc=np.log10(fwhm_n), scale=0.04)
        fwhm_p = np.clip(fwhm_p,0.8,10.)
        fwhm_p = np.array(1.3)
        
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

            #Multiply through by appropriate flat:
            sci, head = fits.getdata(date+'/fits/'+fname+'.fits', header=True)
            sci = sci * flat[j]

            #Add the appropriate dark:
            sci = sci + dark[j]

            #Add the appropriate bias:
            sci = sci + bias[j]

            #Add the overscan:


            #Clip to saturation point
            sci = np.clip(sci,0,65000)

            #Write to file:
            fits.writeto(date+'/fits/'+fname+'.fits', sci, head, overwrite=True)

            #Edit header of output .fits file:
            edit_header(date,ccd,visit,expo,\
                        ccd_ras[j],ccd_decs[j],\
                        mount_ra, mount_dec,\
                        fwhm_p)

            

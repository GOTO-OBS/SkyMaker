from numpy import *
import numpy as np
import os
import csv

from astropy import wcs
import astropy.units as u
from astropy import coordinates as coord

from astroquery.sdss import SDSS
from astroquery.vizier import Vizier

#import matplotlib.pyplot as plt

def makelist(ra, dec, ccd, date, lname, variability=False):

    vname = date+'/var/'+lname+'.var'
    rname = date+'/reg/'+lname+'.reg'
    
    #Define size of camera in pixels:
    xsi = 8176.
    ysi = 6132.

    #Corresponding size in degs:
    wid = xsi*1.24/3600.
    hei = ysi*1.24/3600.

    #Define the columns to get from UCAC, magnitude filter and set unlimited rows:
    v = Vizier(columns=['_RAJ2000', '_DEJ2000', 'Vmag'], column_filters={"Vmag":"<17."})
    v.ROW_LIMIT=-1

    #Query UCAC4
    ucac = v.query_region(\
                            coord.SkyCoord(ra=ra, dec=dec,\
                            unit=(u.deg, u.deg),\
                            frame='icrs'),\
                            height=str(hei+1)+"d", width=str(wid+1)+"d",\
                            catalog=["UCAC4"])[0]

    #For SDSS, need to calculate bounding box:
    declim = dec+np.array([-hei,hei])/2.
    declim_r = np.pi*declim/180.
    wid_r = np.pi*wid/180.
    d_lim_r = 2.*np.arcsin(np.sin(wid_r/4.)/np.cos(declim_r))
    ralim = ra+np.array([-1,1])*np.max(180.*d_lim_r/np.pi)

    #Query SDSS:
    query = "SELECT p.objid, p.ra, p.dec, p.g, "+\
    "p.deVRad_g, p.deVPhi_g, p.deVAB_g, "+\
    "p.type, p.flags_g, (flags & dbo.fPhotoFlags('SATURATED')) as SAT "+\
    "FROM PhotoPrimary AS p "+\
    "WHERE  p.ra BETWEEN "+str(ralim[0])+" AND "+str(ralim[1])+" AND "+\
    "p.dec BETWEEN "+str(declim[0])+" AND "+str(declim[1])+" AND "+\
    "p.g BETWEEN 16 AND 22"# AND "+\
    #"p.htmid*37 & 0x000000000000FFFF < (650 * 10)"
    sdss = SDSS.query_sql(query)

    #Remove saturated sources and separate gals from stars:
    o = np.where(sdss['SAT'] == 0)
    sdss = sdss[o]

    #Find and remove matching sources:
    c = coord.SkyCoord(ra=sdss['ra']*u.degree, dec=sdss['dec']*u.degree)
    cs = coord.SkyCoord(ra=np.array(ucac['_RAJ2000'])*u.degree, \
                        dec=np.array(ucac['_DEJ2000'])*u.degree)
    idx, d2d, d3d = cs.match_to_catalog_sky(c)
    match = np.where(d2d.value*3600. > 1.0)
    ucac = ucac[match]

    #Group stars from SDSS and UCAC together: 
    o = np.where(sdss['type'] == 6)
    sdss_stars = sdss[o]
    stars = np.append(np.array([ucac['_RAJ2000'], ucac['_DEJ2000'], ucac['Vmag']]),\
                      np.array([sdss_stars['ra'], sdss_stars['dec'], sdss_stars['g']]),\
                      axis=1)

    #Extract gals from SDSS:
    o = np.where(sdss['type'] != 6)
    sdss_gals = sdss[o]
    gals = np.array([sdss_gals['ra'], sdss_gals['dec'],\
                     sdss_gals['g'],\
                     sdss_gals['deVRad_g'], sdss_gals['deVAB_g'], sdss_gals['deVPhi_g']])

    #Generate WCS to convert (RA,Dec) to (x,y)
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = np.array([4088, 3066])+np.random.normal(scale=3.,size=2)
    w.wcs.cdelt = np.array([3.444e-4, 3.444e-4])
    w.wcs.crval = [ra, dec] #Pointing position of telescope.
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    pixcrds = np.column_stack((stars[0,:],stars[1,:]))
    pixcrdg = np.column_stack((gals[0,:], gals[1,:]))
    
    worldg =  w.wcs_world2pix(pixcrdg, 1)
    worlds =  w.wcs_world2pix(pixcrds, 1)
    
    gals[0,:] =  worldg[:,0]
    gals[1,:] =  worldg[:,1]
    
    stars[0,:] =  worlds[:,0]
    stars[1,:] =  worlds[:,1]

    ind = [-1]
    if variability:
        #Select 0.5% of m<19 stars and add random variation:
        varyfile = open(vname,'w')
        bright = np.squeeze(np.where((stars[2,:]<19) &\
                                      (stars[0,:]>0) & (stars[0,:]<xsi) &\
                                      (stars[1,:]>0) & (stars[1,:]<ysi)))
        n_vary = int(np.round(0.005*(np.size(bright))))
        if n_vary > 0:
            ind = np.random.choice(bright, n_vary, replace=False)
            orig = stars[2,ind]
            stars[2,ind] = stars[2,ind] + np.random.normal(0, 1, ind.size)
            for i in range(0,np.size(ind)):
                varyfile.write(str(stars[0,ind[i]])+', '+str(stars[1,ind[i]])+', '+\
                               str(orig[i])+', '+str(stars[2,ind[i]])+'\n')
            
    #Write stars to file:
    myfile = open('templist'+ccd+'.list','w')
    regfile = open(rname,'w')
    regfile.write('image\n')
    
    for i in range(0,(stars.shape)[1]):
        myfile.write((str(100)+' '+str(stars[0,i])+' '+str(stars[1,i])+' '+str(stars[2,i]))+'\n')
        if np.any(ind == i):
            regfile.write('circle('+str(stars[0,i])+','+str(stars[1,i])+',3) #color=blue\n')
        else:
            regfile.write('circle('+str(stars[0,i])+','+str(stars[1,i])+',3)\n')
            
    #Write galss to file:
    for i in range(0,(gals.shape)[1]):
        myfile.write((str(200)+' '+str(gals[0,i])+' '+str(gals[1,i])+' '+str(gals[2,i])+' ' +\
                      str(0)+' ' +str(0)+' ' +str(0)+' ' +str(0)+' ' + str(gals[3,i])+' '+\
                      str(gals[4,i])+' '+str(gals[5,i])+'\n'))
        regfile.write('circle('+str(gals[0,i])+','+str(gals[1,i])+',3) #color=red\n')
    
    myfile.close()
    regfile.close()
    
#What needs to be done now is:
# - Generate a list of pointings for the mount, covering ~100sq deg.
#   Note that SDSS only covers certain parts of the sky, so you need
#   to point only in these regions.
# - For each pointing, calculate the pointing for each telescope.
#   There will need to be some overlap between pointings.
# - Generate a filename `fname` for each telescope; see below for example.
#   The source list will have name `fname`.list.
# - Automatically edit the .conf file to include the correct output
#   filename `fname`.fits and its own PSF. Each simulation should
#   have its own .conf file called `fname`.conf
# - We also need to include telescope jitter.
# - Edit the header information to reflect the settings for each pointing.
# - Loop over the mount pointings, generating twelve 2-minute sim images
#   (three for each telescope) per mount pointing. We can
#   use the same PSF for each set of twelve. In fact, we can probably assume
#   the PSf remains stable for about 30 minutes (so 5 mount pointings).  

#fname="GOTO_01_20170323_00001"
#makelist(180.,30.,fname)

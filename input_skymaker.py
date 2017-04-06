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

def makelist(ra, dec, ccd, fname):

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
    result = v.query_region(\
                            coord.SkyCoord(ra=ra, dec=dec,\
                            unit=(u.deg, u.deg),\
                            frame='icrs'),\
                            height=str(hei+1)+"d", width=str(wid+1)+"d",\
                            catalog=["UCAC4"])[0]

    #Extract values from UCAC4
    xs = result['_RAJ2000']
    ys = result['_DEJ2000']
    ms = result['Vmag']

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
    "p.g BETWEEN 13 AND 21"# AND "+\
    #"p.htmid*37 & 0x000000000000FFFF < (650 * 10)"
    res = SDSS.query_sql(query)

    #Extract the values from the SDSS table:
    x = res['ra']
    y = res['dec']
    m = res['g']
    t = res['type']
    a = res['deVRad_g']
    ba = res['deVAB_g']
    pa = res['deVPhi_g']
    
    #Remove the saturated sources:
    ma = np.ma.masked_where(res['SAT']!=0, res['SAT'])
    x = x[~ma.mask].copy()
    y = y[~ma.mask].copy()
    m = m[~ma.mask].copy()
    t = t[~ma.mask].copy()
    a = a[~ma.mask].copy()
    ba = ba[~ma.mask].copy()
    pa = pa[~ma.mask].copy()

    #Find and remove matching sources:
    c = coord.SkyCoord(ra=x*u.degree, dec=y*u.degree)
    cs = coord.SkyCoord(ra=xs, dec=ys)
    idx, d2d, d3d = cs.match_to_catalog_sky(c)
    match = np.where(d2d.value*3600. > 1.0)
    xs = xs[match]
    ys = ys[match]
    ms = ms[match]
    
    #Generate WCS to convert (RA,Dec) to (x,y)
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = np.array([4088, 3066])+np.random.normal(scale=3.,size=2)
    w.wcs.cdelt = np.array([3.444e-4, 3.444e-4])
    w.wcs.crval = [ra, dec] #Pointing position of telescope.
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    
    pixcrd = np.column_stack((x,y))
    pixcrds = np.column_stack((xs,ys))
    
    world =  w.wcs_world2pix(pixcrd, 1)
    worlds =  w.wcs_world2pix(pixcrds, 1)
    
    xp =  world[:,0]
    yp =  world[:,1]
    
    xps =  worlds[:,0]
    yps =  worlds[:,1]
    
    #Write SDSS to file, if type==6, then it's a star, otherwise galaxy:
    myfile = open('templist'+ccd+'.list','w')
    regfile = open(fname+'.reg','w')
    regfile.write('fk5\n')
    for i in range(0,x.size):
        if t[i] == 6:
            myfile.write((str(100)+' '+str(xp[i])+' '+str(yp[i])+' '+str(m[i])+'\n'))
        else:
            myfile.write((str(200)+' '+str(xp[i])+' '+str(yp[i])+' '+str(m[i])+' ' +\
                        str(0)+' ' +str(0)+' ' +str(0)+' ' +str(0)+' ' + str(a[i])+' '+\
                        str(ba[i])+' '+str(pa[i])+'\n'))
        regfile.write('circle('+str(x[i])+','+str(y[i])+',5")\n')

    #Write UCAC4 to file (all stars):
    for i in range(0,xs.size):
            myfile.write((str(100)+' '+str(xps[i])+' '+str(yps[i])+' '+str(ms[i])+'\n'))
            regfile.write('circle('+str(xs[i])+','+str(ys[i])+',5") # color=red\n')              
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

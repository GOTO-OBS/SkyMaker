from numpy import *
import numpy as np
import os
import csv
from astropy import wcs
from astroquery.sdss import SDSS
from astropy import coordinates as coord
from astroquery.vizier import Vizier
import astropy.units as u
arr = np.array(np.loadtxt('result2.csv', delimiter=",", skiprows=1))
Vizier.ROW_LIMIT=-1
v = Vizier(columns=['_RAJ2000', '_DEJ2000', 'Vmag'],column_filters={"Vmag":"<16"})
v.ROW_LIMIT=-1
result = v.query_region(coord.SkyCoord(ra=0.0,dec=0.0,unit=(u.deg,u.deg),frame='icrs'),width='3d',catalog=['I/322A/out'])[0]

query = "SELECT p.objid,p.ra,p.dec,p.g, p.deVRad_g, p.deVPhi_g, p.deVAB_g, p.type, p.flags_g, (flags & dbo.fPhotoFlags('SATURATED')) as SAT  FROM PhotoPrimary AS p  WHERE  p.ra BETWEEN -1.4 AND 1.4  AND p.dec BETWEEN -1.1 AND 1.1 AND p.g BETWEEN 5 AND 20  AND p.htmid*37 & 0x000000000000FFFF < (650 * 10)"
res = SDSS.query_sql(query)
x = res['ra']

y = res['dec']
m = res['g']
t = res['type']
a = res['deVRad_g']
ba = res['deVAB_g']
pa = res['deVPhi_g']

    #mask  = (arr[i,9]==0)
ma = np.ma.masked_where(res['SAT']!=0, res['SAT'])

x = x[~ma.mask].copy()
y = y[~ma.mask].copy()
m = m[~ma.mask].copy()
t = t[~ma.mask].copy()
a = a[~ma.mask].copy()
ba = ba[~ma.mask].copy()
pa = pa[~ma.mask].copy()
xs = result['_RAJ2000']
ys = result['_DEJ2000']
ms = result['Vmag']
print y, ys
w = wcs.WCS(naxis=2)
w.wcs.crpix = [4088, 3066]
w.wcs.cdelt = np.array([3.444e-4, 3.444e-4])
w.wcs.crval = [0., 0.] #Pointing position of telescope.
w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

pixcrd = np.column_stack((x,y))
pixcrds = np.column_stack((xs,ys))

print shape(pixcrd)

world =  w.wcs_world2pix(pixcrd, 1)
worlds =  w.wcs_world2pix(pixcrds, 1)

x =  world[:,0]
y =  world[:,1]

xs =  worlds[:,0]
ys =  worlds[:,1]

myfile= open('config/list.list','w')
for i in range(0,x.size):
    if t[i] ==6:
            myfile.write((str(100)+' '+str(x[i])+' '+str(y[i])+' '+str(m[i])+'\n'))

    else:
        myfile.write((str(200)+' '+str(x[i])+' '+str(y[i])+' '+str(m[i])+' ' +str(0)+' ' +str(0)+' ' +str(0)+' ' +str(0)+' ' + str(a[i])+' ' +str(a[i]*ba[i])+' '+str(pa[i])+'\n'))
for i in range(0,xs.size):
            myfile.write((str(100)+' '+str(xs[i])+' '+str(ys[i])+' '+str(ms[i])+'\n'))
myfile.close()

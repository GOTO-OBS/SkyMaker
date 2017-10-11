from astropy.io import fits
import numpy as np

ccds = np.arange(4, dtype=int)+1
points = np.arange(2)+1
exps = np.arange(3)+1
date = '20170803'

for ccd in ccds:
    for point in points:
        stack = []
        print 'Combining: '+'GOTO_'+\
            "{:02d}".format(ccd)+'_'+\
            date+'_'+\
            "{:04d}".format(point)
        for exp in exps:
            fname = 'GOTO_'+\
                    "{:02d}".format(ccd)+'_'+\
                    date+'_'+\
                    "{:04d}".format(point)+'_'+\
                    "{:02d}".format(exp)+'.fits'
            
            #Read in the file and append:
            if exp == 1:
                head = fits.getheader(date+'/fits/'+fname)
            #head.set('EXPTIME', 3*head['EXPTIME'])
            stack.append(fits.getdata(date+'/fits/'+fname))

        #Median combine the stack:
        medcomb = np.median(stack, axis=0)
        oname = 'GOTO_'+\
                "{:02d}".format(ccd)+'_'+\
                date+'_'+\
                "{:04d}".format(point)+'.fits'
        out = fits.PrimaryHDU(medcomb, header=head)
        out.writeto(date+'/fits/'+oname, overwrite=True)
        

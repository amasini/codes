import numpy as np
from astropy.io import fits
import sys
import os

wd='/Users/alberto/Desktop/XBOOTES/'
#use binned expomap, not vignetting-corrected
binsize=1
if binsize==8:
  di='out'
elif binsize==1:
  di='out2'
#take bkg in 9-12 keV count rate (cts/pixel/s)
bkg912=np.genfromtxt('avg_bkg.dat',skip_header=1,unpack=True, usecols=2)

#from Figure 2 of Hickox&Markevitch 2006, compute count rate in 0.5-7 keV
bkg057=(0.12+0.27+1.1)*bkg912
k=0
for dirs in os.listdir(wd):
    if os.path.isdir(dirs) == True:
        print(dirs)
        for obsid in os.listdir(wd+dirs+'/'):
            if os.path.isdir(wd+dirs+'/'+obsid+'/') == True:
                if obsid != 'AAA':
                    k=k+1
                    print(k,obsid)
                    if len(obsid) == 4:
                        stem='0'+obsid
                    elif len(obsid) == 3:
                        stem='00'+obsid
                    elif len(obsid) == 5:
                        stem=obsid
                    else:
                        print('Something\'s wrong with '+dirs+'/'+obsid+'/')
               	        sys.exit()
                    #open expomap
                    dat=fits.open(dirs+'/'+obsid+'/repro/'+di+'/broad_thresh.expmap')
                    expo=dat[0].data
                    header=dat[0].header
                    #print(header['DATAMODE'])
                    dat.close()
                    #rescale expomap
                    new=expo*bkg057[k-1]*binsize**2
                    #write bkgmap
                    hdu = fits.PrimaryHDU(new,header=header)
                    hdu.writeto(dirs+'/'+obsid+'/repro/'+di+'/acisf'+stem+'_bkgmap.fits',overwrite=True)



import numpy as np
import sys
import os
from astropy.io import fits

wd='/Users/alberto/Desktop/XBOOTES/'
k=0
for dirs in os.listdir(wd):
    if os.path.isdir(dirs) == True:
        print(dirs)
        for obsid in os.listdir(wd+dirs+'/'):
            if os.path.isdir(wd+dirs+'/'+obsid+'/') == True:
                if obsid != 'A':
                    k=k+1
                    print(k,dirs,obsid)
                    if len(obsid) == 4:
                        stem='0'+obsid
                    elif len(obsid) == 3:
                        stem='00'+obsid
                    elif len(obsid) == 5:
                        stem=obsid
                    else:
                        print('Something\'s wrong with '+dirs+'/'+obsid+'/')
               	        sys.exit()
                    filename='acisf'+stem+'_repro_fov1.fits'

                    fovfits=fits.open(wd+dirs+'/'+obsid+'/repro/'+filename)

                    fov=fovfits[1].data
                    fovfits[1].data=fov[fov["CCD_ID"]<4]
                    fovfits.writeto(wd+dirs+'/'+obsid+'/repro/fov_acisI.fits',overwrite=True)
                    fovfits.close()


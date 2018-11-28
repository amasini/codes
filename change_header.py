# from 05/10/18, this script is deprecated. It has been incorporated into the "repro_new_obsid.py" script.

import sys
import numpy as np
from astropy.io import fits
import os

wd='/Users/alberto/Desktop/XBOOTES/'

for dirs in os.listdir(wd):
    if os.path.isdir(dirs) == True:
        for obsid in os.listdir(wd+dirs+'/'):
            if os.path.isdir(wd+dirs+'/'+obsid+'/') == True:
                if obsid=='6999':
                    if len(obsid) == 4:
                        stem='0'+obsid
                    elif len(obsid) == 3:
                        stem='00'+obsid
                    elif len(obsid) == 5:
                        stem=obsid
                    else:
                        print('Something\'s wrong with '+dirs+'/'+obsid+'/')
                        sys.exit()

                    data_head=fits.getheader(wd+dirs+'/'+obsid+'/repro/acisf'+stem+'_repro_evt2.fits', 1)
                    #data_expo=data_head['EXPOSURE']
                    #take evt file header and wcs
                    #head=fits.getheader(wd+dirs+'/'+obsid+'/repro/out2/acisf'+stem+'_expomap.fits', 0)
                    #head['EXPOSURE']=data_expo

                    #replace old bkgmap header with this one
                    data, oldheader = fits.getdata(wd+'sim_full/acisf'+stem+'_sim_poiss.fits', header=True)
                    hdu = fits.PrimaryHDU(data,header=data_head)
                    hdu.writeto(wd+'/sim_full/acisf'+stem+'_sim_poiss_newheader.fits',overwrite=True)              

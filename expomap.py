import numpy as np
import sys
import subprocess as s
import os
from astropy.io import fits

wd='/Users/alberto/Desktop/XBOOTES/'

for dirs in os.listdir(wd):
    if os.path.isdir(dirs) == True:
        print(dirs)
        for obsid in os.listdir(wd+dirs+'/'):
            if os.path.isdir(wd+dirs+'/'+obsid+'/') == True:
                if len(obsid) == 4:
                    stem='0'+obsid
                elif len(obsid) == 3:
                    stem='00'+obsid
                elif len(obsid) == 5:
                    stem=obsid
                else:
                    print('Something\'s wrong with '+dirs+'/'+obsid+'/')
               	    sys.exit()
                #produce bkg maps
                #filename='acisf'+stem+'_repro_9to12keV.fits'
                
                #s.call('fluximage \"'+wd+dirs+'/'+obsid+'/repro/'+filename+'\" '+wd+dirs+'/'+obsid+'/repro/out2/ unit=time binsize=1 clobber=yes',shell=True)

                #produce vignetting-corrected expomaps with units=s
                '''
                filename='acisf'+stem+'_repro_evt2.fits'
                
                s.call('fluximage \"'+wd+dirs+'/'+obsid+'/repro/'+filename+'[events][ccd_id<4]\" '+wd+dirs+'/'+obsid+'/repro/eff_area unit=area binsize=1 clobber=yes',shell=True)
                s.call('fluximage \"'+wd+dirs+'/'+obsid+'/repro/'+filename+'[events][ccd_id<4]\" '+wd+dirs+'/'+obsid+'/repro/expo binsize=1 clobber=yes',shell=True)

                dat=fits.open(dirs+'/'+obsid+'/repro/eff_area_broad_thresh.expmap')
                max_effarea=np.max(dat[0].data)
                dat.close()
                dat2=fits.open(dirs+'/'+obsid+'/repro/expo_broad_thresh.expmap')
                oldexpo=dat2[0].data
                oldhead=dat2[0].header
                dat2.close()
                expo=oldexpo/max_effarea
                hdu=fits.PrimaryHDU(expo,header=oldhead)
                hdu.writeto(dirs+'/'+obsid+'/repro/out2/acisf'+stem+'_expomap.fits',overwrite=True)
                '''

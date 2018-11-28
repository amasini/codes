import numpy as np
import sys
import subprocess as s
import os

wd='/Users/alberto/Desktop/XBOOTES/'
for dirs in os.listdir(wd):
    if os.path.isdir(dirs) == True:
        #if file != 'bootes2':
            #print('Reprocessing in '+wd+'/'+file+'/.')
            #s.call('chandra_repro "'+wd+'/'+file+'/*" outdir=" " verbose=0', shell=True)
        for obsid in os.listdir(wd+dirs+'/'):
            if os.path.isdir(wd+dirs+'/'+obsid+'/') == True:
                print(dirs,obsid)
                if len(obsid) == 4:
                    stem='0'+obsid
                elif len(obsid) == 3:
                    stem='00'+obsid
                elif len(obsid) == 5:
                    stem=obsid
                else:
                    print('Something\'s wrong with '+dirs+'/'+obsid+'/')
               	    sys.exit()
                filename='acisf'+stem+'_repro_'
                #remove='rm -f '+wd+dirs+'/'+obsid+'/repro/'+filename+'9to13keV.fits'
                #s.call(remove,shell=True)
                #command='dmcopy \"'+wd+dirs+'/'+obsid+'/repro/'+filename+'evt2.fits[events][ccd_id=0:3][energy=9000:12000]\" \"'+wd+dirs+'/'+obsid+'/repro/'+filename+'9to12keV.fits\" clobber=yes'
                #s.call(command,shell=True)
                s.call('dmcopy \"'+wd+dirs+'/'+obsid+'/repro/'+filename+'evt2.fits[events][ccd_id=0:3][energy=500:7000]\" \"'+wd+dirs+'/'+obsid+'/repro/'+filename+'05to7keV.fits\" clobber=yes',shell=True)

#[bin x=::1,y=::1] to bin 1:1 in dmcopy

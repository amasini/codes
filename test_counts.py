# THIS SCRIPT CHECKS THE SIMULATIONS' AND DATA COUNTS IN ORDER TO SEE IF THE BACKGROUND 
# MODELING IS REASONABLE
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
'''
wd='/Users/alberto/Desktop/XBOOTES/bootes2/6998/repro/out2/'
bkg=fits.open(wd+'acisf06998_bkgmap.fits')
back=bkg[0].data
sum0=np.sum(back)
bkg.close()

dat=fits.open('/Users/alberto/Desktop/XBOOTES/bootes2/6998/repro/acisf06998_repro_05to7keV.fits')
dat00=dat[1].data
sum00=len(dat00)
dat.close()

a=[]
for i in range(0,1000):
    poiss=np.random.poisson(back)
    a.append(np.sum(poiss))

plt.figure()
plt.hist(a)
plt.axvline(x=sum0,linestyle='dashed',color='red')
plt.axvline(x=sum00,linestyle='dashed',color='green')
plt.show()
'''
sum0,sum00,ratio=[],[],[]
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

                bkg=fits.open(wd+'/sim_full/acisf'+stem+'_sim_poiss.fits')
                sim=bkg[0].data
                sum0.append(np.sum(sim))
                bkg.close()

                dat=fits.open(wd+dirs+'/'+obsid+'/repro/acisf'+stem+'_repro_05to7keV.fits')
                dat00=dat[1].data
                
                sum00.append(len(dat00))
                dat.close()

                rat=(float(len(dat00))-float(np.sum(sim)))/float(np.sum(sim))
                if rat > 1.:
                    print('Check background and data.')
                else:
                    if (rat > 0.5) or (rat < -0.3):
                        print('Check background and data also here.')
                        ratio.append(rat)
                    else:
                        ratio.append(rat)

print(np.median(ratio))
plt.figure()
plt.hist(ratio,bins=20,histtype='step',color='b',linewidth=2)
plt.axvline(x=np.median(ratio),color='r',linestyle='dashed')
plt.xlabel('Data-Bgd/Bgd')
plt.ylabel('N')
plt.tight_layout()
plt.show()
#plt.savefig(wd+'data_back_no-outliers.png',format='png')


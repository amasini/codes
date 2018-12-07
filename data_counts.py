import matplotlib.pyplot as plt
import sys
import subprocess as s
import os
import time
import numpy as np
from astropy.io import fits

wd='/Users/alberto/Desktop/XBOOTES/'
w=open(wd+'data_counts.dat','w')
for dirs in os.listdir(wd):
	if os.path.isdir(wd+dirs) == True:
		if dirs=='data':
			for obsid in os.listdir(wd+dirs+'/'):
				if (os.path.isdir(wd+dirs+'/'+obsid+'/') == True) and (obsid!='junk_output'):
					
					if len(obsid) == 4:
						stem='0'+obsid
					elif len(obsid) == 3:
						stem='00'+obsid
					elif len(obsid) == 5:
						stem=obsid
					else:
						print('Something\'s wrong with '+dirs+'/'+obsid+'/')
						sys.exit()

					dat=fits.open(wd+dirs+'/'+obsid+'/repro_new_asol/acisf'+stem+'_repro_05to7keV.fits')
					datafull=dat[1].data
					dat.close()
					
					dat1=fits.open(wd+dirs+'/'+obsid+'/repro_new_asol/acisf'+stem+'_repro_05to2keV.img')
					datasoft=dat1[0].data
					dat1.close()
					
					dat2=fits.open(wd+dirs+'/'+obsid+'/repro_new_asol/acisf'+stem+'_repro_2to7keV.img')
					datahard=dat2[0].data
					dat2.close()
					
					print(obsid,len(datafull),np.sum(datasoft),np.sum(datahard))
					w.write(dirs+' \t '+obsid+' \t '+str(len(datafull))+' \t '+str(np.sum(datasoft))+' \t '+str(np.sum(datahard))+'\n')
w.close()

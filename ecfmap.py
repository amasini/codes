# Script to compute ECF for CDWFS start for single obsids computing
# for each pixel ECF*expo.
import numpy as np
import subprocess as s
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import time

wd='/Users/alberto/Desktop/XBOOTES/'

obsid,cy_obs=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=[1,6],dtype='str')
cy,ecff,ecfs,ecfh=np.genfromtxt(wd+'ecf_cycles_1.8.dat',unpack=True,skip_header=2)
tstart=time.time()
'''
#####################
###### PART 1 #######
#####################
band='hard'

# Define the ecf for the chosen band for all cycles, Gamma=1.4
if band == 'broad':
	ecflist=ecff
elif band == 'soft':
	ecflist=ecfs
else:
	ecflist=ecfh

# Create 3 psfmaps in the full band: r90 (binned 4x4), r90*expo, r90^2*expo (both binned native Chandra pixel)
for i in range(len(obsid)):
	if len(obsid[i]) == 4:
		stem='0'+obsid[i]
	elif len(obsid[i]) == 3:
		stem='00'+obsid[i]
	elif len(obsid[i]) == 5:
		stem=obsid[i]

	print(obsid[i])

	expomap=wd+'data/'+obsid[i]+'/repro_new_asol/out/acisf'+stem+'_'+band+'_expomap_4reb.fits'
	
	expofile=fits.open(expomap)
	exp=expofile[0].data
	expheader=expofile[0].header
	expofile.close()
	
	mask=np.zeros_like(exp,dtype=float)
	ecf=ecflist[cy==int(cy_obs[i])][0]
	
	mask[exp>0]=ecf
		
	output=mask*exp #use exp*ECF
	hdu1 = fits.PrimaryHDU(output,header=expheader)
	hdu1.writeto(wd+'ecfmaps/'+stem+'_'+band+'_ecf-x-expo_1.8_4reb.fits',overwrite=True)

	#s.call('dmcopy \"'+wd+'/ecfmaps/'+stem+'_broad_ecf-x-expo_1.8.fits[bin (x,y)=::4]\" outfile='+wd+'/ecfmaps/'+stem+'_broad_ecf-x-exp_1.8_4reb.fits clobber=yes',shell=True)
	
print((time.time()-tstart)/60.,'minutes for all obsids.')


sys.exit()
'''

'''
#####################
###### PART 2 #######
#####################
REBIN 4X4 WITH create_images_to_mosaic.py
MOSAIC WITH reproject_image "../ecfmaps/*_broad_ecf-x-exp_4reb.fits" matchfile=../new_mosaics_detection/cdwfs_broad_4reb.fits outfile=../new_mosaics_detection/cdwfs_broad_ecfmap-x-exp_4reb.fits
'''

#####################
###### PART 3 #######
#####################

#now take the r90**2 x expo map and divide it by the total exposure map (watch out for the 
# divisions by 0)
for band in ['broad','soft','hard']:
	mappa=fits.open(wd+'new_mosaics_detection/cdwfs_'+band+'_ecfmap-x-exp_1.8_4reb.fits')
	map=mappa[0].data
	map[np.isnan(map)]=0.0 #put nans to 0
	mappa.close()
	
	expmap=fits.open(wd+'new_mosaics_detection/cdwfs_'+band+'_expomap_4reb.fits')
	exp=expmap[0].data
	expheader=expmap[0].header
	exp[np.isnan(exp)]=0.0 #put nans to 0
	expmap.close()
	
	new=np.zeros_like(exp,dtype=float)
	for i in range(exp.shape[0]):
		for j in range(exp.shape[1]):
			if exp[i][j] != 0.0:
				new[i][j]=map[i][j]/(16.0*exp[i][j]) # need to multiply by 16 because expomaps have been divided by 16
			else:
				new[i][j]=0

	hdu1 = fits.PrimaryHDU(new,header=expheader)
	hdu1.writeto(wd+'new_mosaics_detection/cdwfs_'+band+'_ecfmap_1.8_4reb.fits',overwrite=True)

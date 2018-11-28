# Script to compute psfmap (r90) for XBOOTES start for single obsids computing
# for each pixel r90*expo. r90 can be computed from theta wich is a function of x and y 
# once I know the center of the image. Then I mosaic all the r90*expo and divide the full
# mosaic image by the full exposure image, in order to have the weigthed average of r90
import numpy as np
import subprocess as s
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import time

wd='/Users/alberto/Desktop/XBOOTES/'
obsid=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=1,dtype='str')
'''
###r90 maps 4x4 rebinned
obsid=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=1,dtype='str')
#obsid=['3596']
tstart=time.time()
for obs in obsid:
	if len(obs) == 4:
		stem='0'+obs
	elif len(obs) == 3:
		stem='00'+obs
	elif len(obs) == 5:
		stem=obs

	print(obs)
	#no need to go trough ra and dec if I use the correct transformation from x,y -> theta
	imamap=wd+'data/'+obs+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img'

	imafile=fits.open(imamap)
	ima=imafile[0].data
	imaheader=imafile[0].header
	crpix1=imaheader["CRPIX1"] #aimpoint pixel in x
	crpix2=imaheader["CRPIX2"] #aimpoint pixel in y

	r90=np.zeros_like(ima,dtype=float)
	for i in range(ima.shape[0]):
		for j in range(ima.shape[1]):
			#r90 = 1+theta^2/10. -> write file with r90*expo
			#r90[i][j]=1+(((i-crpix1)**2+(j-crpix2)**2)*6.724e-5)/10. #where 6.724e-5 is (0.492/60)^2
			r90[i][j]=1.+((((i-crpix1)**2+(j-crpix2)**2)*1.07584e-3)/10.) #where 1.07584e-3 is (1.968/60)^2

	output=r90
	hdu1 = fits.PrimaryHDU(output,header=imaheader)
	hdu1.writeto(wd+'psfmaps/'+stem+'_r90_4rebinned.fits',overwrite=True)

print((time.time()-tstart)/60.,'minutes for all obsids.')

sys.exit()


tstart=time.time()
for obs in obsid:
	if len(obs) == 4:
		stem='0'+obs
	elif len(obs) == 3:
		stem='00'+obs
	elif len(obs) == 5:
		stem=obs

	print(obs)
	#no need to go trough ra and dec if I use the correct transformation from x,y -> theta
	expomap=wd+'data/'+obs+'/repro_new_asol/out/acisf'+stem+'_broad_expomap.fits'

	expofile=fits.open(expomap)
	exp=expofile[0].data
	expheader=expofile[0].header
	crpix1=expheader["CRPIX1"] #aimpoint pixel in x
	crpix2=expheader["CRPIX2"] #aimpoint pixel in y

	r90=np.zeros_like(exp,dtype=float)
	for i in range(exp.shape[0]):
		for j in range(exp.shape[1]):
			#r90 = 1+theta^2/10. -> write file with r90*expo
			r90[i][j]=1+(((i-crpix1)**2+(j-crpix2)**2)*6.724e-5)/10.

	#output=r90*r90*exp #use exp*r90^2
	output=r90*exp #use exp*r90
	hdu1 = fits.PrimaryHDU(output,header=expheader)
	hdu1.writeto(wd+'psfmaps/'+stem+'_r90-x-expo.fits',overwrite=True)

print((time.time()-tstart)/60.,'minutes for all obsids.')

sys.exit()
'''
#mosaic of the single obsids! see mosaic.py code

#now take the r90**2 x expo map and divide it by the total exposure map (watch out for the 
# divisions by 0)
'''
for chunk_number in range(0,6):
	
	mappa=fits.open(wd+'chunks_of_mosaics_fullres/mosaic_broad_chunk_'+str(chunk_number)+'_r90.fits')
	map=mappa[0].data
	map[np.isnan(map)]=0.0 #put nans to 0

	expmap=fits.open(wd+'chunks_of_mosaics_fullres/mosaic_broad_chunk_'+str(chunk_number)+'_exp.fits')
	exp=expmap[0].data
	expheader=expmap[0].header
	exp[np.isnan(exp)]=0.0 #put nans to 0

	new=np.zeros_like(exp)
	for i in range(exp.shape[0]):
		for j in range(exp.shape[1]):
			if exp[i][j] != 0.0:
				new[i][j]=map[i][j]/exp[i][j]
			else:
				new[i][j]=0

	hdu1 = fits.PrimaryHDU(new,header=expheader)
	hdu1.writeto(wd+'chunks_of_mosaics_fullres/mosaic_broad_chunk_'+str(chunk_number)+'_averager90.fits',overwrite=True)

'''
#mappa=fits.open(wd+'murray_sens/mosaic_broad_r90squared-x-expo_4rebinned_murray.fits')
mappa=fits.open(wd+'mosaic_broad_r90-x-expo_4rebinned.fits')
map=mappa[0].data
map[np.isnan(map)]=0.0 #put nans to 0

expmap=fits.open(wd+'mosaic_broad_expomap_4rebinned.fits')
exp=expmap[0].data
expheader=expmap[0].header
exp[np.isnan(exp)]=0.0 #put nans to 0

new=np.zeros_like(exp,dtype=float)
for i in range(exp.shape[0]):
	for j in range(exp.shape[1]):
		if exp[i][j] != 0.0:
			new[i][j]=map[i][j]/exp[i][j]
		else:
			new[i][j]=0

hdu1 = fits.PrimaryHDU(new,header=expheader)
hdu1.writeto(wd+'average_r90map.fits',overwrite=True)
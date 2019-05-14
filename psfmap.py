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
#murray=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=5)
#obsid=obsid[murray==1]
#obsid=['3596']
tstart=time.time()

#####################
###### PART 1 #######
#####################
band='hard'
# Create 3 psfmaps in the full band: r90 (binned 4x4), r90*expo, r90^2*expo (both binned native Chandra pixel)
for obs in obsid:
	if len(obs) == 4:
		stem='0'+obs
	elif len(obs) == 3:
		stem='00'+obs
	elif len(obs) == 5:
		stem=obs

	print(obs)

	expomap=wd+'data/'+obs+'/repro_new_asol/out/acisf'+stem+'_'+band+'_expomap.fits'

	expofile=fits.open(expomap)
	exp=expofile[0].data
	expheader=expofile[0].header
	crpix1=expheader["CRPIX1"] #aimpoint pixel in x
	crpix2=expheader["CRPIX2"] #aimpoint pixel in y
	expofile.close()
	
	mask=np.zeros_like(exp,dtype=float)
	mask[exp>0]=1.

	r90=np.zeros_like(exp,dtype=float)
	for i in range(exp.shape[0]):
		for j in range(exp.shape[1]):
			#r90 = 1+theta^2/10. -> write file with r90*expo
			#r90[i][j]=1+(((i-crpix1)**2+(j-crpix2)**2)*6.724e-5)/10.
			r90[i][j]=1.8+(((i-crpix1)**2+(j-crpix2)**2)*6.724e-5)/10.

	output=r90*mask
	hdu1 = fits.PrimaryHDU(output,header=expheader)
	hdu1.writeto(wd+'psfmaps/'+stem+'_'+band+'_r90.fits',overwrite=True)
	
	output=r90*exp #use exp*r90
	hdu1 = fits.PrimaryHDU(output,header=expheader)
	hdu1.writeto(wd+'psfmaps/'+stem+'_'+band+'_r90-x-expo.fits',overwrite=True)
	
	output=r90*r90*exp #use exp*r90^2
	hdu1 = fits.PrimaryHDU(output,header=expheader)
	hdu1.writeto(wd+'psfmaps/'+stem+'_'+band+'_r90squared-x-expo.fits',overwrite=True)
	
print((time.time()-tstart)/60.,'minutes for all obsids.')
sys.exit()
'''
#####################
###### PART 2 #######
#####################
#mosaic the single obsids (see mosaic.py code)
#start from 4221 in a 7200x7200 square
w1=open(wd+'commands_r90.xco','w')
w2=open(wd+'commands_r90-x-expo.xco','w')
w3=open(wd+'commands_r90squared-x-expo.xco','w')

w1.write('read/fits/size=7200/rebin=4 '+wd+'psfmaps/04221_r90.fits\n')
w2.write('read/fits/size=7200/rebin=4 '+wd+'psfmaps/04221_r90-x-expo.fits\n')
w3.write('read/fits/size=7200/rebin=4 '+wd+'psfmaps/04221_r90squared-x-expo.fits\n')

w1.write('save_image\n')
w2.write('save_image\n')
w3.write('save_image\n')

for obs in obsid:
	if len(obs) == 4:
		stem='0'+obs
	elif len(obs) == 3:
		stem='00'+obs
	elif len(obs) == 5:
		stem=obs
	if obs!='4221':	
		w1.write('read/fits/size=7200/rebin=4 '+wd+'psfmaps/'+stem+'_r90.fits\n')
		w2.write('read/fits/size=7200/rebin=4 '+wd+'psfmaps/'+stem+'_r90-x-expo.fits\n')
		w3.write('read/fits/size=7200/rebin=4 '+wd+'psfmaps/'+stem+'_r90squared-x-expo.fits\n')
		
		w1.write('sum_image\n')
		w2.write('sum_image\n')
		w3.write('sum_image\n')
		
		w1.write('save_image\n')
		w2.write('save_image\n')
		w3.write('save_image\n')
	
w1.write('write/fits '+wd+'/r90_4rebinned.fits\n')
w2.write('write/fits '+wd+'/r90-x-expo_4rebinned.fits\n')
w3.write('write/fits '+wd+'/xbootes_r90squared-x-expo_4rebinned.fits\n')

w1.write('exit\n')
w2.write('exit\n')
w3.write('exit\n')

w1.close()
w2.close()
w3.close()

#s.call('ximage @'+wd+'commands_r90.xco',shell=True)
#.call('ximage @'+wd+'commands_r90-x-expo.xco',shell=True)
s.call('ximage @'+wd+'commands_r90squared-x-expo.xco',shell=True)

print((time.time()-tstart)/60.,'minutes for all mosaics.')
sys.exit()
'''
#####################
###### PART 3 #######
#####################

#now take the r90**2 x expo map and divide it by the total exposure map (watch out for the 
# divisions by 0)
mappa=fits.open(wd+'xbootes_r90squared-x-expo_4rebinned.fits')
map=mappa[0].data
map[np.isnan(map)]=0.0 #put nans to 0

expmap=fits.open(wd+'xbootes_broad_expomap_4rebinned.fits')
exp=expmap[0].data
expheader=expmap[0].header
exp[np.isnan(exp)]=0.0 #put nans to 0

new=np.zeros_like(exp,dtype=float)
for i in range(exp.shape[0]):
	for j in range(exp.shape[1]):
		if exp[i][j] != 0.0:
			new[i][j]=map[i][j]/(exp[i][j])
		else:
			new[i][j]=0

hdu1 = fits.PrimaryHDU(new,header=expheader)
hdu1.writeto(wd+'xbootes_r90squared_4rebinned.fits',overwrite=True)

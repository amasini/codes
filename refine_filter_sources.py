# Code to refine list of filtered sources (sources that are actually in the expomap mosaic FOV)
# MADE FOR sim_indep_new/ SIMULATIONS!!
import numpy as np
import sys
import subprocess as s
from astropy.io import fits
import time

wd = '/Users/alberto/Desktop/XBOOTES/'
band = 'hard'

# Take input sources in band (filtered with the stacked fov file)
(f_s,ra_s,dec_s)=np.genfromtxt(wd+'poiss_rand_'+band+'_filtered.dat',skip_header=1,unpack=True,usecols=[0,1,2])

# Define the corners of a inner area to avoid looping on all the sources
x0,x1 = 216.35,219.5
y0,y1 = 32.475,35.725

# This function defines the oblique part of the mosaic to be cut (line passing through two given points)
def eqy(x):
	x00,y00=218.06,32.30
	x11,y11=219.31,33.48
	return (y11-y00)/(x11-x00)*(x-x00)+y00

# Define the mask to speed up the loop
mask = (ra_s >= x0) & (ra_s <= x1) & (dec_s >= y0) & (dec_s <= y1) & (dec_s >= eqy(ra_s))

# These are the points outside the mask (e.g., on the edges), on which to loop
filt_x = ra_s[~mask]
filt_y = dec_s[~mask]
filt_f = f_s[~mask]

# Retain the sources inside the mask and write them in the output file
cpx = ra_s[mask]
cpy = dec_s[mask]
cpf = f_s[mask]
w=open(wd+'poiss_rand_'+band+'_filtered_new.dat','w')
w.write(band+' flux \t RA \t DEC \n')
for j in range(len(cpx)):
	w.write(str(cpf[j])+' \t '+str(cpx[j])+' \t '+str(cpy[j])+'\n')

# Open the flux to countrate map (obtained inverting the countrate to flux map)
ftcmap=wd+'new_mosaics_detection/cdwfs_'+band+'_ftcmap_4reb.fits'
ima=fits.open(ftcmap)
im2=ima[0].data
ima.close()

count0,count1=0,0
tin=time.time()

# Loop on the sources on the edges to find out how many are outside the map
for i in range(len(filt_x)):
	print(i+1,len(filt_x))
	
	# Convert the coordinates (ra,dec) -> logicalx, logicaly
	s.call('dmcoords '+ftcmap+' asol=none opt=cel celfmt=deg ra='+str(filt_x[i])+' dec='+str(filt_y[i])+'',shell=True)
	res=s.check_output('pget dmcoords logicalx logicaly',shell=True)
	logicalx,logicaly=res.splitlines()
	logicalx,logicaly=float(logicalx),float(logicaly)

	if int(round(logicaly)-1) < im2.shape[0]:
		if int(round(logicalx)-1) < im2.shape[1]:
			cf=im2[int(round(logicaly)-1),int(round(logicalx)-1)]
			if cf == 0:
				count0 += 1
			else:
				w.write(str(filt_f[i])+' \t '+str(filt_x[i])+' \t '+str(filt_y[i])+'\n')
		else:
			count1 += 1
	else:
		count1 += 1

w.close()

print(float((time.time()-tin)/3600.),'hours')	
print(len(ra_s),'total sources')
print(count0,'with CF = 0')
print(count1,'which fall out of map edges')
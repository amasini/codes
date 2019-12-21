# Code to refine list of filtered sources (sources that are actually in the expomap mosaic FOV)
# MADE FOR sim_indep_new/ SIMULATIONS!!
import numpy as np
import sys
import subprocess as s
from astropy.io import fits
import time
import matplotlib.pyplot as plt
import shapely
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
from shapely.geometry import Point

wd = '/Users/alberto/Desktop/XBOOTES/'
band = 'hard'

# Take input sources in band (filtered with the stacked fov file)
#(f_s,ra_s,dec_s)=np.genfromtxt(wd+'poiss_rand_'+band+'_filtered.dat',skip_header=1,unpack=True,usecols=[0,1,2])
(f_s,ra_s,dec_s,index,gamma)=np.genfromtxt(wd+'poiss_rand_'+band+'_sim_indep_06-Dec-19.dat',skip_header=1,unpack=True,usecols=[0,1,2,3,4])

vertici=[220.0267948,35.7377765,219.7303732,35.9314269,219.5519014,35.7636720,219.2934582,35.9283171,219.2480182,35.8839873,219.1305235,35.7738744,219.1265918,35.7692047,218.8681778,35.9321660,218.6972598,35.7732633,218.4379386,35.9353177,218.2620080,35.7734558,218.0124660,35.9385023,217.8257475,35.7744323,217.5753390,35.9354395,217.4136532,35.7747358,217.1498440,35.9355801,216.9726266,35.7686856,216.7137137,35.9356757,216.5450738,35.7721707,216.2834802,35.9342766,216.0094940,35.6718321,216.1050535,35.5993514,216.1582446,35.5429560,216.1406036,35.4836556,216.2139380,35.4706769,216.0278235,35.2815291,216.2841327,35.1127412,216.0394191,34.8720152,216.2616422,34.7199489,216.0542906,34.4580421,216.2579338,34.2809554,216.0516934,34.0597323,216.2785841,33.8974689,216.0555553,33.6435332,216.2831962,33.4747427,216.0496764,33.2349607,216.2888133,33.0801220,216.0440277,32.8356026,216.2778931,32.6636406,216.0557721,32.4172337,216.3413017,32.2306153,216.5089491,32.4037988,216.7680693,32.2223728,216.9392303,32.3991121,217.1984293,32.2286001,217.3688933,32.4047599,217.6289103,32.2274704,217.7981150,32.4003258,218.0599154,32.2239739,219.3012120,33.4021530,219.3117554,33.4130557,219.5643741,33.2478125,219.7989135,33.4917700,219.5572610,33.6551546,219.7987687,33.8940856,219.5556077,34.0559728,219.8020684,34.3054690,219.5590559,34.4651185,219.8019599,34.7066805,219.5620980,34.8753852,219.8048289,35.1172671,219.5664249,35.2783838]
vert_x = vertici[0::2]
vert_y = vertici[1::2]

vertici2=list((vert_x[i],vert_y[i]) for i in range(len(vert_x)))

# the 'bounding' polygon of the Bootes field
poly1 = shapely.geometry.Polygon(vertici2)

mask = [poly1.intersects(shapely.geometry.Point(x,y)) for x,y in zip(ra_s,dec_s)]

cpx = ra_s[mask]
cpy = dec_s[mask]
cpf = f_s[mask]
cpi = index[mask]
cpg = gamma[mask]

# Define the corners of a inner area to avoid looping on all the sources
x0,x1 = 216.35,219.5
y0,y1 = 32.475,35.725

# This function defines the oblique part of the mosaic to be cut (line passing through two given points)
def eqy(x):
	x00,y00=218.06,32.30
	x11,y11=219.31,33.48
	return (y11-y00)/(x11-x00)*(x-x00)+y00


# Define the mask to speed up the loop
mask2 = (cpx >= x0) & (cpx <= x1) & (cpy >= y0) & (cpy <= y1) & (cpy >= eqy(cpx))

# These are the points outside the mask (e.g., on the edges), on which to loop
filt_x = cpx[~mask2]
filt_y = cpy[~mask2]
filt_f = cpf[~mask2]
filt_i = cpi[~mask2]
filt_g = cpg[~mask2]

# Retain the sources inside the mask and write them in the output file
keep_x = cpx[mask2]
keep_y = cpy[mask2]
keep_f = cpf[mask2]
keep_i = cpi[mask2]
keep_g = cpg[mask2]

w=open(wd+'poiss_rand_'+band+'_sim_indep_06-Dec-19_filtered_new.dat','w')
w2=open(wd+'poiss_rand_'+band+'_sim_indep_06-Dec-19_filtered_new.reg','w')
w.write(band+' flux \t RA \t DEC \t Index \t Gamma\n')
for j in range(len(keep_x)):
	w.write(str(keep_f[j])+' \t '+str(keep_x[j])+' \t '+str(keep_y[j])+' \t '+str(keep_i[j])+' \t '+str(keep_g[j])+'\n')
	w2.write('circle('+str(keep_x[j])+'d,'+str(keep_y[j])+'d,2\") #color=yellow \n')

# Open the flux to countrate map (obtained inverting the countrate to flux map)
ftcmap=wd+'new_mosaics_detection/cdwfs_'+band+'_expomap_4reb.fits'
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
	
	if (int(round(logicaly)-1) < im2.shape[0]) and (int(round(logicaly)-1) >= 0):
		if (int(round(logicalx)-1) < im2.shape[1]) and (int(round(logicalx)-1) >= 0):
			cf=im2[int(round(logicaly)-1),int(round(logicalx)-1)]
			print(cf)
			if cf == 0:
				count0 += 1
			else:
				w.write(str(filt_f[i])+' \t '+str(filt_x[i])+' \t '+str(filt_y[i])+' \t '+str(filt_i[i])+' \t '+str(filt_g[i])+'\n')
				w2.write('circle('+str(filt_x[i])+'d,'+str(filt_y[i])+'d,2\") #color=yellow \n')
		else:
			count1 += 1
	else:
		count1 += 1
	
w.close()
w2.close()

print(float((time.time()-tin)/3600.),'hours')	
#print(len(ra_s),'total sources')
#print(count0,'with EXP = 0')
#print(count1,'which fall out of map edges')
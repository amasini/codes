#Script to run a new simulation (ie Poissonian realization) of XBOOTES, and to mosaic them in predefined chunks.
import numpy as np
import sys
from astropy.io import fits
import subprocess as s

wd='/Users/alberto/Desktop/XBOOTES/'

obsid=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=1,dtype=str)

for i in range(len(obsid)):
	if len(obsid[i]) == 4:
		stem='0'+obsid[i]
	elif len(obsid[i]) == 3:
		stem='00'+obsid[i]
	elif len(obsid[i]) == 5:
		stem=obsid[i]
	print(i+1)
    #takes simulation
	sim=fits.open(wd+'/sim_soft/acisf'+stem+'_sim.fits')
	simulation=sim[0].data
	backheader=sim[0].header
    
    #creates a new Poissonian realization of the image and saves it
	poiss_sim=np.random.poisson(simulation)
	hdu0 = fits.PrimaryHDU(poiss_sim,header=backheader)
	hdu0.writeto(wd+'/sim_soft/acisf'+stem+'_sim_poiss.fits',overwrite=True)
	poiss_sim2=poiss_sim.astype(float) #this prevents the simulation to have bitpix=64 causing issues with dmextract
	hdu1 = fits.PrimaryHDU(poiss_sim2,header=backheader)
	hdu1.writeto(wd+'/sim_soft/acisf'+stem+'_sim_poiss_bitpix-64.fits',overwrite=True)

sys.exit() #Stop here

#take list of obsids to mosaic in chunks (already ordered from closer to the center)
nchunks=6
for i in range(nchunks):
	obsid=np.genfromtxt(wd+'mosaic_chunk_'+str(i)+'_obsids.dat',unpack=True,dtype='str')
	w=open(wd+'commands.xco','w')
	for j in range(len(obsid)):
		if len(obsid[j]) == 4:
			stem='0'+obsid[j]
		elif len(obsid[j]) == 3:
			stem='00'+obsid[j]
		elif len(obsid[j]) == 5:
			stem=obsid[j]
		if j==0:
			w.write('read/fits/size=16000 '+wd+'sim_full_new/acisf'+stem+'_sim_poiss.fits\n')
			w.write('save_image\n')
		else:
			w.write('read/fits/size=3000 '+wd+'sim_full_new/acisf'+stem+'_sim_poiss.fits\n')
			w.write('sum_image\n')
			w.write('save_image\n')
	w.write('write/fits '+wd+'mosaic_broad_chunk_'+str(i)+'_new.fits\n')
	w.write('exit\n')
	w.close()

	s.call('ximage @'+wd+'commands.xco',shell=True)
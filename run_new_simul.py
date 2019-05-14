#Script to run a new simulation (ie Poissonian realization) of CDWFS
import numpy as np
import sys
from astropy.io import fits
import subprocess as s

wd='/Users/alberto/Desktop/XBOOTES/'

bands=['soft','hard']

for band in bands:
    #takes simulation
	sim=fits.open(wd+'/sim_all/cdwfs_'+band+'_sim_4reb.fits')
	simulation=sim[0].data
	backheader=sim[0].header
	sim.close()
	
	#creates a new Poissonian realization of the image and saves it
	poiss_sim=np.random.poisson(simulation)
	#hdu0 = fits.PrimaryHDU(poiss_sim,header=backheader)
	#hdu0.writeto(wd+'/sim/cdwfs_'+band+'_sim_poiss_4reb.fits',overwrite=True)
	poiss_sim2=poiss_sim.astype(float) #this prevents the simulation to have bitpix=64 causing issues with dmextract
	hdu1 = fits.PrimaryHDU(poiss_sim2,header=backheader)
	hdu1.writeto(wd+'/sim/cdwfs_'+band+'_sim_poiss_4reb.fits',overwrite=True)

# sum soft and hard to get full band
softsim=fits.open(wd+'/sim/cdwfs_soft_sim_poiss_4reb.fits')
softsimulation=softsim[0].data
backheader=softsim[0].header

hardsim=fits.open(wd+'/sim/cdwfs_hard_sim_poiss_4reb.fits')
hardsimulation=hardsim[0].data

fullsimulation=softsimulation+hardsimulation
hdu0 = fits.PrimaryHDU(fullsimulation,header=backheader)
hdu0.writeto(wd+'/sim/cdwfs_broad_sim_poiss_4reb.fits',overwrite=True)
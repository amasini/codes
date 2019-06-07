#Script to run a new simulation (ie Poissonian realization) of CDWFS
import numpy as np
import sys
from astropy.io import fits
import subprocess as s
import time

wd='/Users/alberto/Desktop/XBOOTES/'

nsim=10
tin=time.time()
for i in range(1,nsim):
	
	for band in ['soft','hard']:
    	#takes simulation
		sim=fits.open(wd+'/sim_all_new/cdwfs_'+band+'_sim_4reb.fits')
		simulation=sim[0].data
		backheader=sim[0].header
		sim.close()
	
		#creates a new Poissonian realization of the image and saves it
		poiss_sim=np.random.poisson(simulation)
		#hdu0 = fits.PrimaryHDU(poiss_sim,header=backheader)
		#hdu0.writeto(wd+'/sim/cdwfs_'+band+'_sim_poiss_4reb.fits',overwrite=True)
		poiss_sim2=poiss_sim.astype(float) #this prevents the simulation to have bitpix=64 causing issues with dmextract
		hdu1 = fits.PrimaryHDU(poiss_sim2,header=backheader)
		hdu1.writeto(wd+'/sim_all_new/'+str(i)+'cdwfs_'+band+'_sim_poiss_4reb.fits',overwrite=True)

	# sum soft and hard to get full band
	softsim=fits.open(wd+'/sim_all_new/'+str(i)+'cdwfs_soft_sim_poiss_4reb.fits')
	softsimulation=softsim[0].data
	backheader=softsim[0].header

	hardsim=fits.open(wd+'/sim_all_new/'+str(i)+'cdwfs_hard_sim_poiss_4reb.fits')
	hardsimulation=hardsim[0].data

	fullsimulation=softsimulation+hardsimulation
	hdu0 = fits.PrimaryHDU(fullsimulation,header=backheader)
	hdu0.writeto(wd+'/sim_all_new/'+str(i)+'cdwfs_broad_sim_poiss_4reb.fits',overwrite=True)
	
	for band in ['broad','soft','hard']:
		# Run wavedetect for each band
		s.call('wavdetect infile='+wd+'/sim_all_new/'+str(i)+'cdwfs_'+band+'_sim_poiss_4reb.fits bkginput='+wd+'/new_mosaics_detection/cdwfs_'+band+'_bkgmap_4reb.fits psffile= outfile='+wd+'/sim_all_new/'+str(i)+'cdwfs_'+band+'_sim_src.fits scellfile=out_scell.fits imagefile=out_ima.fits defnbkgfile=out_defnbkg.fits regfile='+wd+'/sim_all_new/'+str(i)+'cdwfs_'+band+'_sim_src.reg scales="1 2 4" sigthresh=5e-5 clobber=yes',shell=True)
		
		# call compute_fluxes_v2 script to do the first computation
		s.call('python '+wd+'/codes/compute_fluxes_data_v2.py '+band+' '+str(i), shell=True)
		
		# call clean_multiple_v2 script to do the cleaning from multiply-detected sources
		s.call('python '+wd+'/codes/clean_multiple_sources_v2.py '+band+' '+str(i), shell=True)
		
		# call reliability script to work out the reliability curves
		s.call('python '+wd+'/codes/reliability.py '+band+' '+str(i), shell=True)

print((time.time()-tin)/3600.,'hours for the whole run of ',nsim,'simulations.')
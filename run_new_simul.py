#Script to run a new simulation (ie Poissonian realization) of CDWFS
import numpy as np
import sys
from astropy.io import fits
import subprocess as s
import time
import os.path
import smtplib
from email.mime.text import MIMEText

wd='/Users/alberto/Desktop/XBOOTES/'

simfolder = 'sim_indep_12-Dec-19/'
simfolder2 = 'sim_indep_12-Dec-19/'

nsim=10
tin=time.time()
for i in range(0,nsim):
	
	bands=['broad','soft']
	for band in bands:
		
		if os.path.isfile(wd+simfolder+'cdwfs_'+band+'_sim_4reb.fits') == False:
			print('Mosaic of '+band+' simulation is missing. Doing it now.')
			s.call('reproject_image \"'+wd+simfolder2+'acisf*_'+band+'_sim_4reb.fits" matchfile='+wd+'/new_mosaics_detection/cdwfs_broad_4reb.fits outfile='+wd+simfolder+'cdwfs_'+band+'_sim_4reb.fits',shell=True)
			msg = MIMEText(str((time.time()-tin)/3600.)+' hours for the mosaic.')
			
			me = 'alberto311290@gmail.com'
			you = 'alberto.masini@dartmouth.edu'
			msg['Subject'] = 'Code has done the '+band+' mosaic.'
			msg['From'] = me
			msg['To'] = you

			# Send the message via our own SMTP server, but don't include the
			# envelope header.
			s0 = smtplib.SMTP('localhost')
			s0.sendmail(me, [you], msg.as_string())
			s0.quit()
		
			#sys.exit()
			
    	#takes simulation
		sim=fits.open(wd+simfolder+'cdwfs_'+band+'_sim_4reb.fits')
		simulation=sim[0].data
		backheader=sim[0].header
		sim.close()
	
		#creates a new Poissonian realization of the image and saves it
		poiss_sim=np.random.poisson(simulation)
		poiss_sim2=poiss_sim.astype(float) #this prevents the simulation to have bitpix=64 causing issues with dmextract
		hdu1 = fits.PrimaryHDU(poiss_sim2,header=backheader)
		hdu1.writeto(wd+simfolder+str(i)+'cdwfs_'+band+'_sim_poiss_4reb.fits',overwrite=True)
		
		# Run wavdetect to detect sources
		s.call('wavdetect infile='+wd+simfolder+str(i)+'cdwfs_'+band+'_sim_poiss_4reb.fits expfile='+wd+'/new_mosaics_detection/cdwfs_'+band+'_expomap_4reb.fits psffile='+wd+'/new_mosaics_detection/cdwfs_'+band+'_r90_4reb.fits outfile='+wd+simfolder+str(i)+'cdwfs_'+band+'_sim_src_exp-psf.fits scellfile=out_scell.fits imagefile=out_ima.fits defnbkgfile=out_defnbkg.fits regfile='+wd+simfolder+str(i)+'cdwfs_'+band+'_sim_src_exp-psf.reg expthresh=0.01 scales="1.4 2 4" sigthresh=5e-5 clobber=yes',shell=True)
		
		# call compute_fluxes_v2 script to do the aperture photometry
		s.call('python '+wd+'/codes/compute_fluxes_sim.py '+band+' '+str(i), shell=True)
		
		# call clean_multiple_v2 script to do the cleaning from multiply-detected sources (almost useless, but still)
		s.call('python '+wd+'/codes/clean_multiple_sources_sim.py '+band+' '+str(i), shell=True)
		
print((time.time()-tin)/3600.,'hours for the whole run of ',nsim,'simulations.')

# Send a mail when done
msg = MIMEText(str((time.time()-tin)/3600.)+' hours for the whole run of '+str(nsim)+' simulations.')

me = 'alberto311290@gmail.com'
you = 'alberto.masini@dartmouth.edu'
msg['Subject'] = 'Code is done!'
msg['From'] = me
msg['To'] = you

s0 = smtplib.SMTP('localhost')
s0.sendmail(me, [you], msg.as_string())
s0.quit()

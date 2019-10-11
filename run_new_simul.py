#Script to run a new simulation (ie Poissonian realization) of CDWFS
import numpy as np
import sys
from astropy.io import fits
import subprocess as s
import time
import os.path

# Import smtplib for the actual sending function
import smtplib

# Import the email modules we'll need
from email.mime.text import MIMEText

wd='/Users/alberto/Desktop/XBOOTES/'

simfolder = 'sim_indep/'

nsim=1
tin=time.time()
for i in range(0,nsim):
	'''
	bands=['soft']
	for band in bands:
		
		if os.path.isfile(wd+simfolder+'cdwfs_'+band+'_sim_4reb.fits') == False:
			print('Mosaic of '+band+' simulation is missing. Doing it now.')
			s.call('reproject_image \"'+wd+simfolder+'*_'+band+'_sim_4reb.fits" matchfile='+wd+'/new_mosaics_detection/cdwfs_broad_4reb.fits outfile='+wd+simfolder+'cdwfs_'+band+'_sim_4reb.fits',shell=True)
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
	'''
	'''
	if len(bands) == 2:
		# sum soft and hard to get full band
		softsim=fits.open(wd+simfolder+str(i)+'cdwfs_soft_sim_poiss_4reb.fits')
		softsimulation=softsim[0].data
		backheader=softsim[0].header

		hardsim=fits.open(wd+simfolder+str(i)+'cdwfs_hard_sim_poiss_4reb.fits')
		hardsimulation=hardsim[0].data

		fullsimulation=softsimulation+hardsimulation
		hdu0 = fits.PrimaryHDU(fullsimulation,header=backheader)
		hdu0.writeto(wd+simfolder+str(i)+'cdwfs_broad_sim_poiss_4reb.fits',overwrite=True)
	'''
	for band in ['hard']:
		# Run wavedetect for each band
		#s.call('wavdetect infile='+wd+simfolder+str(i)+'cdwfs_'+band+'_sim_poiss_4reb.fits bkginput='+wd+'/new_mosaics_detection/cdwfs_'+band+'_bkgmap_4reb.fits psffile=' '  outfile='+wd+simfolder+str(i)+'cdwfs_'+band+'_sim_src_exp-psf.fits scellfile=out_scell.fits imagefile=out_ima.fits defnbkgfile=out_defnbkg.fits regfile='+wd+simfolder+str(i)+'cdwfs_'+band+'_sim_src_exp-psf.reg scales="1 2 4" sigthresh=5e-5 clobber=yes',shell=True)
		s.call('wavdetect infile='+wd+simfolder+str(i)+'cdwfs_'+band+'_sim_poiss_4reb.fits  bkginput='+wd+'/new_mosaics_detection/cdwfs_'+band+'_bkgmap_4reb.fits expfile='+wd+'/new_mosaics_detection/cdwfs_'+band+'_expomap_4reb.fits psffile='+wd+'/new_mosaics_detection/cdwfs_'+band+'_r90_4reb.fits outfile='+wd+simfolder+str(i)+'cdwfs_'+band+'_sim_src_exp-psf.fits scellfile=out_scell.fits imagefile=out_ima.fits defnbkgfile=out_defnbkg.fits regfile='+wd+simfolder+str(i)+'cdwfs_'+band+'_sim_src_exp-psf.reg expthresh=0.01 scales="1.4 2 4" sigthresh=5e-5 clobber=yes',shell=True)
		
		# call compute_fluxes_v2 script to do the first computation
		s.call('python '+wd+'/codes/compute_fluxes_sim.py '+band+' '+str(i), shell=True)
		
		# call clean_multiple_v2 script to do the cleaning from multiply-detected sources
		s.call('python '+wd+'/codes/clean_multiple_sources_sim.py '+band+' '+str(i), shell=True)
		
		# call reliability script to work out the reliability curves
		#s.call('python '+wd+'/codes/reliability.py '+band+' '+str(i), shell=True)

print((time.time()-tin)/3600.,'hours for the whole run of ',nsim,'simulations.')

sys.exit()

# Open a plain text file for reading.  For this example, assume that
# the text file contains only ASCII characters.
#textfile='prova.txt'
#fp = open(textfile, 'rb')
# Create a text/plain message
msg = MIMEText(str((time.time()-tin)/3600.)+' hours for the whole run of '+str(nsim)+' simulations.')
#fp.close()

# me == the sender's email address
# you == the recipient's email address
me = 'alberto311290@gmail.com'
you = 'alberto.masini@dartmouth.edu'
msg['Subject'] = 'Code is done!'
msg['From'] = me
msg['To'] = you

# Send the message via our own SMTP server, but don't include the
# envelope header.
s0 = smtplib.SMTP('localhost')
s0.sendmail(me, [you], msg.as_string())
s0.quit()

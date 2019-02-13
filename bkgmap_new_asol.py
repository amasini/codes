# Create quiescent instrumental background from 9-12 keV (flare-filtered) evt files, following
# Hickox & Markevitch 06 conversions to F,S,H bands 
import numpy as np
import sys
import subprocess as s
from astropy.io import fits

wd='/Users/alberto/Desktop/XBOOTES/'

obsid=np.genfromtxt(wd+'data_counts.dat',unpack=True, usecols=1, dtype='str')#
band=['broad']

R=0.1375 #radius of FoV in degrees
dmax=0.0416667 #degrees; 2.5 arcmin
radius=dmax*3600 #150"

#w=open(wd+'avg_bkg_new.dat','w')
#w.write('obsid \t MJD \t sur_bri[cts/pixel2/sec] \t err_sur_bri\n')
bkg912file=np.genfromtxt(wd+'avg_bkg_new.dat',unpack=True,usecols=2,dtype='float',skip_header=1)
for i in range(len(obsid)):
	if len(obsid[i]) == 4:
		stem='0'+obsid[i]
	elif len(obsid[i]) == 3:
		stem='00'+obsid[i]
	elif len(obsid[i]) == 5:
		stem=obsid[i]

	print(obsid[i])
	
	filename='acisf'+stem+'_repro_'
	# Creates 9-12 keV evt file to get background
	#s.call('dmcopy \"'+wd+'data/'+obsid[i]+'/repro_new_asol/'+filename+'evt2.fits[events][ccd_id=0:3][energy=9000:12000]\" \"'+wd+'data/'+obsid[i]+'/repro_new_asol/'+filename+'9to12keV.fits\" clobber=yes',shell=True)

	filename='acisf'+stem+'_repro_9to12keV_cl.fits'
	# Extract useful info from file
	ra=s.check_output('dmkeypar '+wd+'data/'+obsid[i]+'/repro_new_asol/'+filename+' RA_PNT echo=yes',shell=True)
	dec=s.check_output('dmkeypar '+wd+'data/'+obsid[i]+'/repro_new_asol/'+filename+' DEC_PNT echo=yes',shell=True)
	exp=s.check_output('dmkeypar '+wd+'data/'+obsid[i]+'/repro_new_asol/'+filename+' LIVETIME echo=yes',shell=True)
	mjd=s.check_output('dmkeypar '+wd+'data/'+obsid[i]+'/repro_new_asol/'+filename+' MJD_OBS echo=yes',shell=True)
	bkg_sur_bri=[]
	# Extract surf brightness 100 times from random positions on the field
	w=open(wd+'instrback_region_'+obsid[i]+'.reg','w')
	for k in range(100):
		randra=np.random.uniform(float(ra)+(0.9*R-dmax),float(ra)-(0.9*R-dmax))
		randdec=np.random.uniform(float(dec)+(0.9*R-dmax),float(dec)-(0.9*R-dmax))
		w.write('circle('+str(randra)+'d, '+str(randdec)+'d,'+str(radius)+'\")\n')
		
		s.call('punlearn dmextract',shell=True)
		s.call('dmextract \''+wd+'data/'+obsid[i]+'/repro_new_asol/'+filename+'[bin sky=circle('+str(randra)+'d,'+str(randdec)+'d,'+str(radius)+'\")]\' outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
		out=s.check_output('dmlist \'counts.fits[col SUR_BRI]\' data,clean | grep -v SUR_',shell=True)
		print(float(out)/float(exp))
		bkg_sur_bri.append(float(out)/float(exp))
		#s.call('rm -f counts.fits',shell=True)
	w.close()
	sys.exit()
    # Use the median and the standard deviation as uncertainty           
	bkg912=np.median(bkg_sur_bri)
	e_bkg912=np.std(bkg_sur_bri)
	#print('='*15)
	#print('The median 9-12 keV bkg surface brightness for '+obsid[i]+' is')
	#print(bkg912,e_bkg912)
	#print('='*15)
	#w.write(obsid[i]+' \t '+str(float(mjd))+' \t '+str(bkg912)+' \t '+str(e_bkg912)+'\n')
	
	bkg912=bkg912file[i]
	# Creates the bkgmap from the 9-12 kev one
	for j in range(len(band)):
		if band[j] == 'broad':
			instrbkg=(0.13+0.27+1.12)*bkg912 #these are cts/s/pixel**2 (F,S,H)
			#coeff=0.19 #these are cts/s/chip; each acis-I chip has 1024x1024 pixels
		elif band[j] == '0.5-2':
			instrbkg=(0.13+0.27)*bkg912 #these are cts/s/pixel**2 (S)
			#coeff=0.07 #these are cts/s/chip; each acis-I chip has 1024x1024 pixels
		elif band[j] == 'hard':
			#instrbkg=1.15*bkg912 #these are cts/s/pixel**2 (H)
			instrbkg=1.12*bkg912 #these are cts/s/pixel**2 (H)
			#coeff=0.12 #these are cts/s/chip; each acis-I chip has 1024x1024 pixels
		# Take the NON-VIGNETTING CORRECTED expomap
		dat=fits.open(wd+'data/'+obsid[i]+'/repro_new_asol/out/'+band[j]+'_thresh.expmap')
		expo=dat[0].data
		header=dat[0].header
		dat.close()
		
		# Rescale expomap
		new=expo*instrbkg
		
		#totbkg=expo*(coeff/1024**2)
		#print(np.sum(totbkg),np.sum(new))
		#cxb=totbkg-new
		#if np.sum(cxb)>0.:
		#	totbkgmap=new+cxb
		#else:
		#	totbkgmap=new
		#	print('In obsid '+obsid[i]+' the cxb in '+band[j]+' band was not added due to high bkg')
		
		# Write bkgmap
		hdu0 = fits.PrimaryHDU(new,header=header)
		hdu0.writeto(wd+'data/'+obsid[i]+'/repro_new_asol/out/acisf'+stem+'_'+band[j]+'_bkgmap_2.fits',overwrite=True)
		#hdu = fits.PrimaryHDU(totbkgmap,header=header)
		#hdu.writeto(wd+'data/'+obsid[i]+'/repro_new_asol/out/acisf'+stem+'_'+band[j]+'_bkgmap_new.fits',overwrite=True)

#w.close()
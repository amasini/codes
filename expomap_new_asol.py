# THIS SCRIPT CREATES EXPOSURE MAPS USING FLUXIMAGE
############### 
# WARNING!!!! #
############### 
#THE SOFT BAND IS (0.5-1.2 keV) DIFFERENT FROM WHAT I WAS ASSUMING (0.5-2 keV) AND ALL MAPS 
# MUST BE DONE AGAIN FOR THIS BAND - FULL AND HARD ARE OK
import numpy as np
import sys
import subprocess as s
from astropy.io import fits

wd='/Users/alberto/Desktop/XBOOTES/'

obsid=np.genfromtxt(wd+'data_counts.dat',unpack=True, usecols=1, dtype='str')
#obsid=['19676','19777','19778','19779']

for i in range(len(obsid)):
	if len(obsid[i]) == 4:
		stem='0'+obsid[i]
	elif len(obsid[i]) == 3:
		stem='00'+obsid[i]
	elif len(obsid[i]) == 5:
		stem=obsid[i]
	
	#produce exposure map in seconds NOT vignetting-corrected
	s.call('fluximage \"'+wd+'data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_cl.fits[events][ccd_id<4]\" '+wd+'data/'+obsid[i]+'/repro_new_asol/out/ units=time bands="0.5:2:1.5" binsize=1 clobber=yes',shell=True)
	
	#produce eff area map, plus exposure map in default units vignetting-corrected
	s.call('fluximage \"'+wd+'data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_cl.fits[events][ccd_id<4]\" '+wd+'data/'+obsid[i]+'/repro_new_asol/eff_area units=area bands="0.5:2:1.5" binsize=1 clobber=yes',shell=True)
	s.call('fluximage \"'+wd+'data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_cl.fits[events][ccd_id<4]\" '+wd+'data/'+obsid[i]+'/repro_new_asol/expo bands="0.5:2:1.5" binsize=1 clobber=yes',shell=True)
	
	#divide the default expomap by the maximum eff area to get the vign corrected expomap in seconds
	for band in ['0.5-2']:
		dat=fits.open(wd+'data/'+obsid[i]+'/repro_new_asol/eff_area_'+band+'_thresh.expmap')
		max_effarea=np.max(dat[0].data)
		dat.close()
		dat2=fits.open(wd+'data/'+obsid[i]+'/repro_new_asol/expo_'+band+'_thresh.expmap')
		oldexpo=dat2[0].data
		oldhead=dat2[0].header
		dat2.close()
		expo=oldexpo/max_effarea
		hdu=fits.PrimaryHDU(expo,header=oldhead)
		hdu.writeto(wd+'data/'+obsid[i]+'/repro_new_asol/out/acisf'+stem+'_'+band+'_expomap.fits',overwrite=True)
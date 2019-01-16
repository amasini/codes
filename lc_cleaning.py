# Script to filter 9-12 keV evt files on a 3sigma clipping of the lightcurve 
import numpy as np
import sys
import subprocess as s
from lightcurves import *

wd='/Users/alberto/Desktop/XBOOTES/'

obs=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=1,dtype='str')

l0,l1=[],[]
for i in range(len(obs)):
	if len(obs[i]) == 5:
		stem=obs[i]
	elif len(obs[i]) == 4:
		stem='0'+obs[i]
	elif len(obs[i]) == 3:
		stem='00'+obs[i]
	
	print(obs[i])
	'''
	s.call('punlearn dmextract',shell=True)
	# Create the LC
	s.call('dmextract infile="'+wd+'/data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV.fits[bin time=::100]" opt="ltc1" outfile="'+wd+'/data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_lc_05to7keV.fits" clobber=yes',shell=True)
	
	# Filter the LC with a default 3sigma clipping and creates GTi file
	#s.call('deflare "'+wd+'/data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_lc_9to12keV.fits" outfile="'+wd+'/data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_lc_9to12keV.gti" method=sigma plot=false verbose=0',shell=True)
	if obs[i] != '19652':
		lc_sigma_clip(wd+'/data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_lc_05to7keV.fits', outfile=wd+'/data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_lc_05to7keV.gti', plot=False, verbose=0)
	else:
		lc_sigma_clip(wd+'/data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_lc_05to7keV.fits', outfile=wd+'/data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_lc_05to7keV.gti')
	# Filter the original file with the GTI file
	s.call('dmcopy "'+wd+'/data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV.fits[@'+wd+'/data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_lc_05to7keV.gti]" '+wd+'/data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_cl.fits clobber=yes',shell=True)
	'''
	# Compare before/after LIVETIMEs
	res=s.check_output('dmkeypar '+wd+'/data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV.fits LIVETIME echo=yes',shell=True)
	l0.append(float(res))
	
	res=s.check_output('dmkeypar '+wd+'/data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_cl.fits LIVETIME echo=yes',shell=True)
	l1.append(float(res))
	

sum0=np.sum(np.array(l0))
sum1=np.sum(np.array(l1))
print('Total exposure time (in Ms) of CDWFS before, after, difference (in ks)')
print(sum0/1e6,sum1/1e6,(sum0-sum1)/1e3)

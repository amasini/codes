import numpy as np
import os
import sys
import subprocess as s

wd='/Users/alberto/Desktop/XBOOTES/'

obsid=np.genfromtxt(wd+'data_counts.dat',usecols=1,unpack=True,dtype='str')

for i in range(len(obsid)):
	if len(obsid[i]) == 4:
		stem='0'+obsid[i]
	elif len(obsid[i]) == 3:
		stem='00'+obsid[i]
	elif len(obsid[i]) == 5:
		stem=obsid[i]
			
	#check to have all the files for the mosaics
	#if os.path.isfile(wd+'sim_full/acisf'+stem+'_sim_poiss.fits'):
	#	pass
	#else:
	#	print('Please create simulation for '+obsid[i]+'')
			
	if os.path.isfile(wd+'data/'+obsid[i]+'/repro_new_asol/out/acisf'+stem+'_0.5-2_bkgmap_total.fits'):
		pass
	else:
		print('Please create bkgmap for '+obsid[i]+'')
	'''	
	if os.path.isfile(wd+'data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV.img'):
		pass
	else:
		print('Creating 05-7 keV img for '+obsid[i]+'')
		s.call('dmcopy \"'+wd+'data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_evt2.fits[ccd_id=0:3][energy=500:7000]\" '+wd+'data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV.img opt=image clobber=yes',shell=True)
		
	if os.path.isfile(wd+'data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_05to2keV.img'):
		pass
	else:
		print('Creating 05-2 keV img for '+obsid[i]+'')
		s.call('dmcopy \"'+wd+'data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_evt2.fits[ccd_id=0:3][energy=500:2000]\" '+wd+'data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_05to2keV.img opt=image clobber=yes',shell=True)
	
	if os.path.isfile(wd+'data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_2to7keV.img'):
		pass
	else:
		print('Creating 2-7 keV img for '+obsid[i]+'')
		s.call('dmcopy \"'+wd+'data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_evt2.fits[ccd_id=0:3][energy=2000:7000]\" '+wd+'data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_2to7keV.img opt=image clobber=yes',shell=True)
	'''
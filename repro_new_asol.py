import numpy as np
import sys
import os.path
import subprocess as s

wd='/Users/alberto/Desktop/XBOOTES/'

#(dirs,obsid)=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=[0,1],dtype='str')
dirs=['CDWFS']
obsid=['19676','19777','19778','19779']

for i in range(len(obsid)):
	if len(obsid[i]) == 4:
		stem='0'+obsid[i]
	elif len(obsid[i]) == 3:
		stem='00'+obsid[i]
	elif len(obsid[i]) == 5:
		stem=obsid[i]
		
	if os.path.isfile(wd+'data/'+obsid[i]+'/acisf'+stem+'_new_asol.fits'):
		pass
	else:
		print('Error, new aspect file not found for obsid '+dirs[i]+'/'+obsid[i]+'.')
	
	s.call('punlearn ardlib',shell=True)
	if i!='A':
		filename=s.check_output('ls '+wd+'data/'+obsid[i]+'/primary/pcadf*_asol1.fits',shell=True)
		s.call('mv '+filename[:-1]+' '+filename[:-12]+'_original.fits',shell=True)
		s.call('cp '+wd+'/data/'+obsid[i]+'/acisf'+stem+'_new_asol.fits '+wd+'/data/'+obsid[i]+'/primary/ ',shell=True)
		s.call('mv '+wd+'/data/'+obsid[i]+'/primary/acisf'+stem+'_new_asol.fits '+wd+'/data/'+obsid[i]+'/primary/acisf'+stem+'_new_asol1.fits',shell=True)
		
		s.call('chandra_repro indir='+wd+'/data/'+obsid[i]+'/ outdir='+wd+'/data/'+obsid[i]+'/repro_new_asol check_vf_pha=yes > '+wd+'/data/'+obsid[i]+'/repro_new_asol.txt',shell=True)
		s.call('dmkeypar '+wd+'/data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_evt2.fits ASOLFILE echo+',shell=True)

	#creates full band evt file binned to native pixel scale
	s.call('dmcopy \"'+wd+'data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_evt2.fits[events][ccd_id=0:3][energy=500:7000]\" \"'+wd+'data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV.fits\" clobber=yes',shell=True)
	#creates full band *image* binned 4x4 to native pixel scale
	s.call('dmcopy \"'+wd+'data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV.fits[events][bin x=::4,y=::4]\" '+wd+'data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img opt=image clobber=yes',shell=True)
	#creates 9-12 keV image to get background
	s.call('dmcopy \"'+wd+'data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_evt2.fits[events][ccd_id=0:3][energy=9000:12000]\" \"'+wd+'data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_9to12keV.fits\" clobber=yes',shell=True)
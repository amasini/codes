# script to reproject full resolution background images to 4x4 0.5-7 keV images using
# reproject_image infile matchfile=../acisf03597_repro_05to7keV_4rebinned.img outfile=
import numpy as np
import sys
import subprocess as s

wd='/Users/alberto/Desktop/XBOOTES/'

band='broad'
band2='05to7keV'
obsid=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=1,dtype='str')

for obs in obsid:
	if len(obs) == 4:
		stem='0'+obs
	elif len(obs) == 3:
		stem='00'+obs
	elif len(obs) == 5:
		stem=obs
	
	infile=wd+'data/'+obs+'/repro_new_asol/out/acisf'+stem+'_'+band+'_bkgmap.fits'
	matchfile=wd+'data/'+obs+'/repro_new_asol/acisf'+stem+'_repro_'+band2+'_4rebinned.img'
	outfile=wd+'data/'+obs+'/repro_new_asol/out/acisf'+stem+'_'+band+'_bkgmap_4rebinned.img'
	
	print(obs)
	s.call('reproject_image '+infile+' matchfile='+matchfile+' outfile='+outfile+'',shell=True)
import numpy as np
import subprocess as s
import sys


wd='/Users/alberto/Desktop/XBOOTES/new_mosaics_detection/'
date = '200113'

# RUN WAVDETECT
for band in ['broad','soft','hard']:

	s.call('punlearn wavdetect', shell=True)
	infile=wd+'cdwfs_'+band+'_4reb.fits'
	#bkginput=wd+'cdwfs_'+band+'_bkgmap_4reb.fits'
	expfile=wd+'cdwfs_'+band+'_expomap_4reb.fits'
	psffile=wd+'cdwfs_'+band+'_r90_4reb.fits'
	
	
	outfile=wd+'cdwfs_'+band+'_src_'+date+'.fits'
	regfile=wd+'cdwfs_'+band+'_src_'+date+'.reg'

	scellfile=wd+'out_scell.fits'
	imagefile=wd+'out_ima.fits'
	defnbkgfile=wd+'out_defnbkg.fits'
	
	scales='\"1.4 2 4\"'
	expthresh='0.01'
	sigthresh='5e-5'
	clobber='yes'
	
	print('Running wavdetect in the '+band+' band...')
	s.call('wavdetect infile='+infile+' expfile='+expfile+' psffile='+psffile+' outfile='+outfile+' scellfile='+scellfile+' imagefile='+imagefile+' defnbkgfile='+defnbkgfile+' regfile='+regfile+' expthresh='+expthresh+' scales='+scales+' sigthresh='+sigthresh+' clobber='+clobber+'',shell=True)

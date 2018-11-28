#create psfmaps and merge them together following Option #3 of the thread: Minimum PSF size
import numpy as np
import subprocess as s
import sys

wd='/Users/alberto/Desktop/XBOOTES/'

(dirs,obsid)=np.genfromtxt(wd+'mosaic_chunk_1_obsids.dat',unpack=True,dtype='str')

for i in range(len(obsid)):
	if len(obsid[i]) == 4:
		stem='0'+obsid[i]
	elif len(obsid[i]) == 3:
		stem='00'+obsid[i]
	elif len(obsid[i]) == 5:
		stem=obsid[i]
	s.call('mkpsfmap '+wd+dirs[i]+'/'+obsid[i]+'/repro/expo_broad_thresh.img outfile='+wd+dirs[i]+'/'+obsid[i]+'/repro/out2/acisf'+stem+'_psfmap.fits energy=1.4967 ecf=0.393',shell=True)
	s.call('dmimgthresh '+wd+dirs[i]+'/'+obsid[i]+'/repro/out2/acisf'+stem+'_psfmap.fits '+wd+'/psfmaps/'+str(i)+'_chunk_1_psfmap.fits expfile='+wd+dirs[i]+'/'+obsid[i]+'/repro/out2/broad_thresh.expmap cut=1 value=INDEF clob+',shell=True)
	
s.call('dmimgfilt \"'+wd+'/psfmaps/*_chunk_1_psfmap.fits\" '+wd+'min_chunk_1.psfmap min \"point(0,0)\"',shell=True)
s.call('dmhedit '+wd+'min_chunk_1.psfmap file= op=add key=BUNIT value=\"arcsec\"',shell=True)
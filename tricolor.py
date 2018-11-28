#create exposure corrected images for tricolor mosaic, in 0.5-2 keV, 2-4.5 keV, 4.5-7 keV
# the effective energies at which the expo map is computed are quite random...
import numpy as np
import subprocess as s

wd='/Users/alberto/Desktop/XBOOTES/'

obsid=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=1,dtype='str')

for i in range(len(obsid)):
	if len(obsid[i]) == 4:
		stem='0'+obsid[i]
	elif len(obsid[i]) == 3:
		stem='00'+obsid[i]
	elif len(obsid[i]) == 5:
		stem=obsid[i]
	s.call('fluximage \"'+wd+'data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_evt2.fits[events][ccd_id<4]\" '+wd+'data/'+obsid[i]+'/repro_new_asol/expo bands="0.5:2.0:1.49,2.0:4.5:4.0,4.5:7.0:6.0" binsize=4 clobber=yes',shell=True)
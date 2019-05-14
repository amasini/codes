# Code to extract background spectra from data and stowed background in CDWFS - part2
import numpy as np
import sys
import subprocess as s
from astropy.io import fits
import os
import matplotlib.pyplot as plt

def distance(pointa, pointb):
	xx = np.cos(pointa[1]/180*3.141592)
	return np.sqrt(((pointa[0]-pointb[0])*xx*3600)**2 +((pointa[1]-pointb[1])*3600)**2)
	
wd='/Users/alberto/Desktop/XBOOTES/'

obsid=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=1,dtype='str')
#obsid=['10450']

'''/// THIS PART CANNOT BE EXECUTED IN A TERMINAL WITH CIAO ///'''
for obs in obsid:
	if len(obs) == 4:
		stem='0'+obs
	elif len(obs) == 3:
		stem='00'+obs
	elif len(obs) == 5:
		stem=obs
	print(obs)
	
	### LOAD STOWED BACK AND DATA BACK IN XSPEC AND PRINT OUT THE COUNT RATES IN 9-12 KEV BAND
	# Change working directory
	os.chdir(wd+'data/'+obs+'/repro_new_asol/extract/')
	'''
	w=open(wd+'data/'+obs+'/repro_new_asol/extract/stow-to-data-ratio.xcm','w')
	w.write('data 1:1 spec_broad.pi \n')
	w.write('data 2:2 '+stem+'_stowed_cp.pha \n')
	w.write('resp 2:2 '+stem+'_stowed.rmf \n')

	w.write('setpl ene \n')

	w.write('ign 1-2:**-9.0 12.0-** \n')
	#w.write('cpd /xw \n')
	#w.write('pl lda \n')

	w.write('set id [open stow-to-data-ratio.dat w]\n')
	w.write('puts $id "[tcloutr rate 1]"\n') # Write the data count rate
	w.write('puts $id "[tcloutr rate 2]"\n') # Write the stowed count rate
	w.write('close $id\n')

	w.write('exit\n')
	w.close()

	s.call('xspec - stow-to-data-ratio.xcm',shell=True)
	
	# Take the ratio from the file
	cr=np.genfromtxt('stow-to-data-ratio.dat',unpack=True,usecols=0)
	dat=float(cr[0])
	stow=float(cr[1])
	ratio=dat/stow

	### COPY STOWED SPECTRUM AND RENORMALIZE ITS COUNTS BY THE RATIO DATA/STOWED
	s.call('cp '+wd+'data/'+obs+'/repro_new_asol/extract/'+stem+'_stowed.pha '+wd+'data/'+obs+'/repro_new_asol/extract/'+stem+'_stowed_cp.pha',shell=True)
	with fits.open(wd+'data/'+obs+'/repro_new_asol/extract/'+stem+'_stowed_cp.pha', mode='update') as hdul:
		# Change something in hdul.
		hdul[1].data['COUNTS']=ratio*hdul[1].data['COUNTS']
		hdul.flush()  # changes are written back to '+stem+'_stowed_cp.pha
	'''
	### FIT THE SPECTRA IN XSPEC
	w=open(wd+'data/'+obs+'/repro_new_asol/extract/fit.xcm','w')
	w.write('data 1:1 spec_broad.pi \n')
	w.write('back 1:1 '+stem+'_stowed_cp.pha \n')

	w.write('ign 1-2:**-0.3 5.0-** \n')
	w.write('setpl ene \n')
	w.write('statistic cstat \n')

	w.write('mo apec+pha*pow \n')
	w.write('0.18 -1\n') # kT 1
	w.write('1 -1 \n') # Ab 2
	w.write('0 -1 \n') # z 3
	w.write('1e-4 \n') # Norm 4

	w.write('0.013 -1 \n') # Nh 5
	w.write('1.4 -1 \n') # Gam 6
	w.write('1e-5 \n') # Norm 7

	w.write('query yes \n')
	w.write('renorm \n')
	w.write('fit \n')

	w.write('set id [open fitparams.dat w]\n')
	w.write('puts $id "[tcloutr param 1]"\n') # Write the kT value
	w.write('puts $id "[tcloutr param 4]"\n') # Write the Norm APEC
	w.write('puts $id "[tcloutr param 7]"\n') # Write the Norm CXB

	w.write('new 7 0 -1 \n') # Shut off the CXB
	w.write('puts $id "[tcloutr rate 1]"\n') # Write the APEC count rate

	w.write('new 7 1e-5 \n') # Put back the CXB
	w.write('thaw 7 \n') # Put back the CXB
	w.write('fit \n')
	w.write('new 4 0 -1 \n') # Shut off the APEC
	w.write('puts $id "[tcloutr rate 1]"\n') # Write the CXB count rate

	w.write('close $id\n')
	w.write('exit\n')
	w.close()

	s.call('xspec - fit.xcm',shell=True)
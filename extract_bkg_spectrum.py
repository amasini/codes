# Code to extract background spectra in the F,S,H bands in the CDWFS
import numpy as np
import sys
import subprocess as s
from astropy.io import fits
import os

wd='/Users/alberto/Desktop/XBOOTES/'

obs=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=1,dtype='str')

for i in range(len(obs)):
	if len(obs[i]) == 5:
		stem=obs[i]
	elif len(obs[i]) == 4:
		stem='0'+obs[i]
	elif len(obs[i]) == 3:
		stem='00'+obs[i]
	
	for band in ['broad']:
		if band == 'broad':
			band2 = '05to7'
		elif band == 'soft':
			band2 = '05to2'
		else:
			band2 == '2to7'
		#'''
		# Create appropriate regionfile
		w=open(wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_'+band+'_src-bkg.reg','w')
		file1=open(wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_'+band+'_src.reg','r')
		for line in file1:
			w.write(line)
		file1.close()
		file2=open(wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_'+band+'_bkg.reg','r')
		for line in file2:
			w.write('-'+line)
		file2.close()
		w.close()
		
		# Extract spectrum with specextract
		s.call('specextract "'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_'+band2+'keV_cl.fits[sky=region('+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_'+band+'_src-bkg.reg)]" '+wd+'data/'+obs[i]+'/repro_new_asol/extract/spec_'+band+' grouptype=NONE binspec=NONE clobber=yes',shell=True)
		'''
		
		# Change working directory
		os.chdir(wd+'data/'+obs[i]+'/repro_new_asol/extract/')
		
		# Create XSPEC command file (fit with apec+powerlaw)
		w=open(wd+'data/'+obs[i]+'/repro_new_asol/extract/fit_bkg_spec_'+band+'.xcm','w')
		w.write('data 1:1 spec_'+band+'.pi\n')
		w.write('ign bad\n')
		w.write('setpl ene\n')
		if band == 'broad':
			w.write('ign 1:**-0.5 7.0-**\n')
		elif band == 'soft':
			w.write('ign 1:**-0.5 2.0-**\n')
		else:
			w.write('ign 1:**-2.0 7.0-**\n')
		
		w.write('model apec+pha*pow\n')
		w.write('0.18\n') # 1
		w.write('\n') # 2
		w.write('0 -1\n') # 3
		w.write('1e-4\n') # 4
		w.write('0.013 -1\n') # 5
		w.write('0.5\n') # 6
		w.write('1e-4\n') # 7
		
		w.write('set statistic cstat\n')	
		w.write('query yes\n')
		w.write('renorm\n')
		w.write('fit\n')
		w.write('set id [open out_'+band+'.dat w]\n')
		w.write('puts $id "[tcloutr rate 1]"\n') # Print the count rates of data, unc, model
		w.write('new 7 0\n')
		w.write('puts $id "[tcloutr rate 1]"\n') # Print the count rates of data, unc, model
		w.write('puts $id "[tcloutr expos 1]"\n') # Print the exposure of spectrum
		w.write('close $id\n') # Closes the output file
		#w.write('cpd /xw\n')
		#w.write('setpl area\n')
		#w.write('pl ldata del\n')
		
		w.write('exit\n')
		w.close()
		
		s.call('xspec - fit_bkg_spec_'+band+'.xcm',shell=True)
		
		# Change working directory
		os.chdir(wd+'codes/')
		
		# Get the results from the out.dat file
		(data_cr,model_cr)=np.genfromtxt(wd+'data/'+obs[i]+'/repro_new_asol/extract/out_'+band+'.dat',unpack=True,skip_footer=2,usecols=[0,2])
		apec_cr=np.genfromtxt(wd+'data/'+obs[i]+'/repro_new_asol/extract/out_'+band+'.dat',unpack=True,skip_footer=1,skip_header=1,usecols=2)
		exp=np.genfromtxt(wd+'data/'+obs[i]+'/repro_new_asol/extract/out_'+band+'.dat',unpack=True,skip_header=2,usecols=0)
		new_apec_cr=(data_cr/model_cr)*apec_cr
		print(data_cr,model_cr,apec_cr,exp,new_apec_cr,new_apec_cr*exp)
		#'''
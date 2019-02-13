# Code to extract background spectra in the F,S,H bands in the CDWFS
import numpy as np
import sys
import subprocess as s
from astropy.io import fits
import os
import matplotlib.pyplot as plt

wd='/Users/alberto/Desktop/XBOOTES/'

obs=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=1,dtype='str')
#obs=['10450']
apec=[]

w2=open(wd+'apec_cr_soft_with-gauss-line.dat','w')
w2.write('obsid \t apec_cr \t pl_cr \t exp \t area_pix \t apec_surf_bri-cts/s/pix \t apec_cr*exp \t pl_cr*exp\n')

#w2=open(wd+'IB_soft.dat','w')
#w2.write('obsid \t area_pix \t IB_cts\n')
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
		'''
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
		# Extract PIXEL Area
		s.call('dmextract infile="'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to2keV_4rebinned.img[bin sky=region('+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_broad_src-bkg.reg)]" outfile=prova_totbkg.fits opt=generic mode=h clobber=yes',shell=True)
		src_area=s.check_output('dmlist "prova_totbkg.fits[cols AREA]" data,clean | grep -v AREA',shell=True)
		
		# Extract INSTR BACK
		#s.call('dmextract infile="'+wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_0.5-2_bkgmap_instr.fits[bin sky=region('+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_broad_src-bkg.reg)]" outfile=prova_totbkg.fits opt=generic mode=h clobber=yes',shell=True)
		#src_area=s.check_output('dmlist "prova_totbkg.fits[cols AREA]" data,clean | grep -v AREA',shell=True)
		#src_cts=s.check_output('dmlist "prova_totbkg.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
		
		# Extract MJD; warning: XSPEC doesn't work if called with a CIAO tool
		#mjd=s.check_output('dmkeypar '+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_'+band2+'keV_cl.fits MJD_OBS echo=yes',shell=True)
		'''
		# Change working directory
		os.chdir(wd+'data/'+obs[i]+'/repro_new_asol/extract/')
		
		# Create XSPEC command file (fit with apec+powerlaw)
		w=open(wd+'data/'+obs[i]+'/repro_new_asol/extract/fit_bkg_spec_'+band+'.xcm','w')
		w.write('data 1:1 spec_'+band+'.pi\n')
		w.write('ign bad\n')
		w.write('setpl ene\n')

		w.write('ign 1:**-0.5 7.0-**\n')
		w.write('model apec+pha*pow+gauss\n')
		w.write('0.18\n') # 1
		w.write('\n') # 2
		w.write('0 -1\n') # 3
		w.write('1e-4\n') # 4
		w.write('0.013 -1\n') # 5
		w.write('0.5\n') # 6
		w.write('1e-4\n') # 7
		w.write('2.18\n') # 8
		w.write('0.1\n') # 9
		w.write('1e-5\n') # 10
		
		w.write('statistic cstat\n')	
		w.write('query yes\n')
		w.write('renorm\n')
		w.write('fit\n')
		
		w.write('fre 8 9 10\n')
		w.write('ign 1:**-0.5 2.0-**\n')
		w.write('fit\n')
		
		w.write('set id [open out_'+band+'.dat w]\n')
		w.write('puts $id "[tcloutr rate 1]"\n') # Print the count rates of data, unc, model
		
		# SHUT OFF PL - TAKE APEC
		w.write('new 7 0\n')
		w.write('puts $id "[tcloutr rate 1]"\n') # Print the count rates of data, unc, model
		
		# SHUT OFF APEC - TAKE PL
		w.write('new 7 1e-4\n')
		w.write('renorm\n')
		w.write('fit\n')
		w.write('new 4 0\n')
		w.write('puts $id "[tcloutr rate 1]"\n') # Print the count rates of data, unc, model
		
		w.write('puts $id "[tcloutr expos 1]"\n') # Print the exposure of spectrum
		w.write('close $id\n') # Closes the output file
		
		#w.write('cpd /xw\n')
		#w.write('setpl area\n')
		#w.write('pl ldata del\n')
		
		w.write('exit\n')
		w.close()
		
		s.call('xspec - fit_bkg_spec_'+band+'.xcm',shell=True)

		#
		# Change working directory
		os.chdir(wd+'codes/')
		'''
		# Get the results from the out.dat file
		(data_cr,model_cr)=np.genfromtxt(wd+'data/'+obs[i]+'/repro_new_asol/extract/out_'+band+'.dat',unpack=True,skip_footer=3,usecols=[0,2])
		apec_cr=np.genfromtxt(wd+'data/'+obs[i]+'/repro_new_asol/extract/out_'+band+'.dat',unpack=True,skip_footer=2,skip_header=1,usecols=2)
		pl_cr=np.genfromtxt(wd+'data/'+obs[i]+'/repro_new_asol/extract/out_'+band+'.dat',unpack=True,skip_footer=1,skip_header=2,usecols=2)
		exp=np.genfromtxt(wd+'data/'+obs[i]+'/repro_new_asol/extract/out_'+band+'.dat',unpack=True,skip_header=3,usecols=0)
		new_apec_cr=(data_cr/model_cr)*apec_cr
		new_pl_cr=(data_cr/model_cr)*pl_cr
		apec.append(new_apec_cr)
		print(data_cr,model_cr,apec_cr,exp,new_apec_cr,new_apec_cr*exp)
		w2.write(obs[i]+' \t '+str(new_apec_cr)+' \t '+str(new_pl_cr)+' \t '+str(exp)+' \t '+str(float(src_area))+' \t '+str(new_apec_cr/float(src_area))+' \t '+str(new_apec_cr*exp)+' \t '+str(new_pl_cr*exp)+'\n')
		
		#print(obs[i],src_cts)
		#w2.write(obs[i]+' \t '+str(float(src_area))+' \t '+str(float(src_cts))+'\n')

w2.close()
sys.exit()

plt.figure()
plt.hist(apec,bins=15)
plt.show()

time=np.genfromtxt(wd+'avg_bkg_new.dat',unpack=True,usecols=1,skip_header=1)

plt.figure()
plt.plot(time,apec,'go')
plt.yscale('log')
plt.show()

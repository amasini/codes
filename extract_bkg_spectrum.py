# Code to extract background spectra in the F,S,H bands in the CDWFS
import numpy as np
import sys
import subprocess as s
from astropy.io import fits
import os
import matplotlib.pyplot as plt

wd='/Users/alberto/Desktop/XBOOTES/'

####
### REPROJECTING THE STOWED BACKGROUND
#reproject_events infile=../../../stowed_2005.fits outfile=03130_stowed.fits aspect=acisf03130_new_asol1.fits match=acisf03130_repro_evt2.fits random=0 clobber=yes

### CHECK CALIBRATION BETWEEN DATA AND STOWED BACK FILES
#dmkeypar acisf03130_repro_evt2.fits GAINFILE echo+
#dmkeypar 03130_stowed.fits GAINFILE echo+

### FILTERING APPROPRIATE FOR VFAINT MODE
#dmcopy "03130_stowed.fits[status=0]" 03130_stowed_cleaned.fits

### EXTRACTING SPECTRUM
#dmextract "03130_stowed.fits[sky=region(acisf03130_broad_src-bkg.reg)][bin pi=1:1024:1]" 03130_stowed.pha wmap="det=8" clobber=yes

### MAKE THE RMF
#mkacisrmf infile=CALDB outfile=03130_stowed.rmf energy=0.3:12.0:0.01 channel=1:1024:1 chantype=PI wmap=03130_stowed.pha gain=CALDB asolfile=acisf03130_new_asol1.fits

#with fits.open(wd+'data/3130/repro_new_asol/03130_stowed_cp.pha', mode='update') as hdul:
#	# Change something in hdul.
#	hdul[1].data['COUNTS']=0.7087*hdul[1].data['COUNTS']
#	hdul.flush()  # changes are written back to original.fits
####


#obs=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=1,dtype='str')
obs=['3130']
apec=[]

#w2=open(wd+'apec_cr_soft_with-gauss-line_NEW.dat','w')
#w2.write('obsid \t apec_cr \t pl_cr \t exp \t area_pix \t apec_surf_bri-cts/s/pix \t apec_cr*exp \t pl_surf_bri-cts/s/pix \t pl_cr*exp\n')

w2=open(wd+'BACK_soft.dat','w')
w2.write('obsid \t CXB PL counrate \t APEC countrate\n')
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
		#s.call('dmextract infile="'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to2keV_4rebinned.img[bin sky=region('+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_broad_src-bkg.reg)]" outfile=prova_totbkg.fits opt=generic mode=h clobber=yes',shell=True)
		#src_area=s.check_output('dmlist "prova_totbkg.fits[cols AREA]" data,clean | grep -v AREA',shell=True)
		
		# Extract INSTR BACK
		#s.call('dmextract infile="'+wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_0.5-2_bkgmap_instr.fits[bin sky=region('+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_broad_src-bkg.reg)]" outfile=prova_totbkg.fits opt=generic mode=h clobber=yes',shell=True)
		#src_area=s.check_output('dmlist "prova_totbkg.fits[cols AREA]" data,clean | grep -v AREA',shell=True)
		#src_cts=s.check_output('dmlist "prova_totbkg.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
		
		# Extract MJD; warning: XSPEC doesn't work if called with a CIAO tool
		#mjd=s.check_output('dmkeypar '+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_'+band2+'keV_cl.fits MJD_OBS echo=yes',shell=True)
		
		# Change working directory
		os.chdir(wd+'data/'+obs[i]+'/repro_new_asol/extract/')
		
		# Create XSPEC command file (fit with apec+powerlaw)
		w=open(wd+'data/'+obs[i]+'/repro_new_asol/extract/fit_bkg_spec_'+band+'.xcm','w')
		w.write('data 1:1 spec_'+band+'.pi\n')
		w.write('ign bad\n')
		w.write('setpl ene\n')

		w.write('ign 1:**-0.5 6.9-**\n')
		#w.write('model apec+pha*pow\n')
		w.write('model apec+pha*pow+pow+gauss+gauss+gauss+gauss+gauss\n')
		w.write('0.18 -1\n') # 1  apec
		w.write('\n') # 2
		w.write('0 -1\n') # 3
		w.write('1e-4\n') # 4
		
		w.write('0.013 -1\n') # 5 pha
		w.write('1.4 -1\n') # 6 pow CXB
		w.write('1e-4\n') # 7
		
		w.write('-2.5\n') # 8  pow INSTR
		w.write('1e-4\n') # 9
		
		w.write('1.124 -1\n') # 10  lines INSTR
		w.write('0.136 -1\n') # 11
		w.write('1e-5\n') # 12
		
		w.write('1.486 -1\n') # 13
		w.write('0.0492 -1\n') # 14
		w.write('1e-5\n') # 15
		
		w.write('1.814 -1\n') # 16
		w.write('0.132 -1\n') # 17
		w.write('1e-5\n') # 18
		
		w.write('5.9 -1\n') # 19
		w.write('0.0331 -1\n') # 20
		w.write('1e-5\n') # 21
		
		w.write('2.19 -1\n') # 22
		w.write('0.1 -1\n') # 23
		w.write('1e-5\n') # 24
		
		w.write('statistic cstat\n')	
		w.write('query yes\n')
		w.write('renorm\n')
		w.write('fit\n')
		
		#############
		w.write('cpd /xw\n')
		w.write('setpl area\n')
		w.write('pl ldata del\n')
		#############

		#########################################################		
		# SHUT OFF PL - TAKE APEC
		'''
		w.write('set id [open out_'+band+'_apec.dat w]\n')
		w.write('puts $id "[tcloutr param 7]"\n') # Write the CXB PL norm (shouldn't vary)
		
		w.write('fre 7 8 9 12 15 18 21 24\n') # freeze powerlaws and lines
		w.write('ign 1:**-0.5 1.0-**\n')
		w.write('fit\n')
		
		w.write('puts $id "[tcloutr param 4]"\n') # Write the APEC norm after refitting in the 0.5-1 keV range
				
		w.write('delcomp 2-8\n')
		
		w.write('puts $id "[tcloutr rate 1]"\n') # Write the APEC count rate
		w.write('close $id\n') # Closes the output file
		
		w.write('exit\n')
		'''
		#########################################################
		# SHUT OFF APEC - TAKE PL
		'''
		w.write('set id [open out_'+band+'_pl.dat w]\n')
		w.write('puts $id "[tcloutr param 4]"\n') # Write the APEC norm
		
		w.write('fre 4 8 9 12 15 18 21 24\n') # freeze powerlaws and lines
		w.write('ign 1:**-0.5 1.0-**\n')
		w.write('fit\n')
		
		w.write('puts $id "[tcloutr param 7]"\n') # Write the CXB PL norm after refitting in the 0.5-1 keV range
				
		w.write('delcomp 1\n')
		w.write('delcomp 3-8\n')
		
		w.write('puts $id "[tcloutr rate 1]"\n') # Write the CXB count rate
		w.write('close $id\n') # Closes the output file
		'''
		#w.write('exit\n')

		w.close()
		
		s.call('xspec - fit_bkg_spec_'+band+'.xcm',shell=True)
		
		# Change working directory
		os.chdir(wd+'codes/')
		
		# Get the results from the out.dat file
		'''
		(data_cr,model_cr)=np.genfromtxt(wd+'data/'+obs[i]+'/repro_new_asol/extract/out_'+band+'.dat',unpack=True,skip_footer=2,usecols=[0,2])
		apec_cr=np.genfromtxt(wd+'data/'+obs[i]+'/repro_new_asol/extract/out_'+band+'.dat',unpack=True,skip_footer=1,skip_header=1,usecols=2)
		exp=np.genfromtxt(wd+'data/'+obs[i]+'/repro_new_asol/extract/out_'+band+'.dat',unpack=True,skip_header=2,usecols=0)
		pl_cr=np.genfromtxt(wd+'data/'+obs[i]+'/repro_new_asol/extract/out_'+band+'_pl.dat',unpack=True,skip_footer=1,skip_header=1,usecols=2)
		new_apec_cr=(data_cr/model_cr)*apec_cr
		new_pl_cr=(data_cr/model_cr)*pl_cr
		apec.append(new_apec_cr)
		#assume src_area=1649136 pixels
		src_area=1649136
		print(data_cr,model_cr,apec_cr,exp,new_apec_cr,new_apec_cr*exp)
		w2.write(obs[i]+' \t '+str(new_apec_cr)+' \t '+str(new_pl_cr)+' \t '+str(exp)+' \t '+str(float(src_area))+' \t '+str(new_apec_cr/float(src_area))+' \t '+str(new_apec_cr*exp)+' \t '+str(new_pl_cr/float(src_area))+' \t '+str(new_pl_cr*exp)+'\n')
		
		#print(obs[i],src_cts)
		#w2.write(obs[i]+' \t '+str(float(src_area))+' \t '+str(float(src_cts))+'\n')
		'''
		apeccr=np.genfromtxt(wd+'data/'+obs[i]+'/repro_new_asol/extract/out_'+band+'_apec.dat',unpack=True,skip_header=2,usecols=2)
		cxbcr=np.genfromtxt(wd+'data/'+obs[i]+'/repro_new_asol/extract/out_'+band+'_pl.dat',unpack=True,skip_header=2,usecols=2)
		w2.write(obs[i]+' \t '+str(cxbcr)+' \t '+str(apeccr)+'\n')
w2.close()
sys.exit()

time=np.genfromtxt(wd+'avg_bkg_new.dat',unpack=True,usecols=1,skip_header=1)

plt.figure()
plt.plot(time,apec,'go')
plt.yscale('log')
plt.show()

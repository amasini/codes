# Code to compute the CDWFS positional errors - THIS SCRIPT IS OBSOLETE
import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from ciao_contrib.region.check_fov import FOVFiles
import subprocess as s 
import os
from scipy.optimize import curve_fit

def func(x, a, b):
	return a * x**(b) + 0.4

wd = '/Users/alberto/Desktop/XBOOTES/'

# interpolate r50 from chandra proposer guide at 1.49 keV (fig 4.12)
(oa,r50)=np.genfromtxt(wd+'r50_1.49keV.txt',unpack=True) # oa is in arcminutes and r50 in arcsec
oa_grid=np.linspace(min(oa),max(oa),200)
r50_grid = np.interp(oa_grid, oa, r50)

poserr=[]

# take the cdwfs catalog
cat=fits.open(wd+'cdwfs_merged_cat1.fits')
ra=cat[1].data['RA']
dec=cat[1].data['DEC']
probf=cat[1].data['PROB_F']
probs=cat[1].data['PROB_S']
probh=cat[1].data['PROB_H']
cat.close()

# take the list of obsids
obs=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=1,dtype='str')
used_r90=0
# for each source, pick the obsid in which it is closest to aimpoint
for i in range(len(ra)):
	print(i+1)
	# Define the reference band (i.e., from which band the ra and dec come from?)
	if probf[i] != 9999:
		band = 'broad'
		band2 = '05to7keV'
	elif probs[i] != 9999:
		band = '0.5-2'
		band2 = '05to2keV'
	else:
		band = 'hard'
		band2 = '2to7keV'
	
	# Find the list of obsids which contain the source
	my_obs = FOVFiles(wd+'data/*/repro_new_asol/fov_acisI.fits')
	myobs = my_obs.inside(ra[i], dec[i])
	
	# if the source is contained is just one obsid, no need to choose: compute r50 and extract counts
	if len(myobs) == 1:
	
		obs=myobs[0][36:-30]
		if len(obs) == 4:
			stem='0'+obs
		elif len(obs) == 3:
			stem='00'+obs
		elif len(obs) == 5:
			stem=obs
		path=myobs[0][:-14]+'/out/acisf'+stem+'_broad_expomap.fits'
		s.call('punlearn dmcoords',shell=True)
		s.call('dmcoords '+path+' asol=none opt=cel celfmt=deg ra='+str(ra[i])+' dec='+str(dec[i])+'',shell=True)
		theta=s.check_output('pget dmcoords theta',shell=True)
		theta=float(theta)
		
		diff=abs(oa_grid-theta)
		r50_new=r50_grid[diff==min(diff)]
		r50_new=r50_new[0]
		
		# try to extract total and bkg counts from mosaics
		imagemap=wd+'mosaic_broad_4rebinned.fits'
		#imagemap=myobs[0][:-14]+'acisf'+stem+'_repro_05to7keV_4rebinned.img'
		backmap=wd+'mosaic_broad_bkgmap_4rebinned.fits'
		#backmap=myobs[0][:-14]+'out/acisf'+stem+'_broad_bkgmap_total.fits'

		s.call('dmextract infile="'+imagemap+'[bin pos=circle('+str(ra[i])+'d,'+str(dec[i])+'d,'+str(r50_new*0.000277778)+'d)]" mode=h bkg="'+backmap+'[bin pos=circle('+str(ra[i])+'d,'+str(dec[i])+'d,'+str(r50_new*0.000277778)+'d)]" outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
		cts=s.check_output('dmlist "counts.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
		bkg=s.check_output('dmlist "counts.fits[cols BG_COUNTS]" data,clean | grep -v BG_COUNTS',shell=True)
		cts,bkg=float(cts),float(bkg)
		net=cts-bkg
		if net > 0.:
			poserr.append(r50_new/np.sqrt(net))
		else:
			used_r90=used_r90+1
			print('#'*20)
			print('r50 not suitable here, try with r90')
			print('#'*20)
			
			if band != 'hard':
				r90=1+10*(theta/10.)**2
			else:
				r90=1.8+10*(theta/10.)**2
			
			s.call('dmextract infile="'+imagemap+'[bin pos=circle('+str(ra[i])+'d,'+str(dec[i])+'d,'+str(r90*0.000277778)+'d)]" mode=h bkg="'+backmap+'[bin pos=circle('+str(ra[i])+'d,'+str(dec[i])+'d,'+str(r90*0.000277778)+'d)]" outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
			cts=s.check_output('dmlist "counts.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
			bkg=s.check_output('dmlist "counts.fits[cols BG_COUNTS]" data,clean | grep -v BG_COUNTS',shell=True)
			cts,bkg=float(cts),float(bkg)
			net=cts-bkg
			if net > 0.:
				poserr.append(r90/np.sqrt(net))
			else:
				poserr.append(0.1)
				print('#'*20)
				print('Not even r90 is suitable here',ra[i],dec[i],obs,r90,cts,bkg,net)
				print('#'*20)
				#sys.exit()
		#if band != 'hard':
		#	r90=1+10*(theta/10.)**2
		#else:
		#	r90=1.8+10*(theta/10.)**2
		#r50=r90*5./9.
		
		'''
		imagemap=myobs[0][:-14]+'acisf'+stem+'_repro_'+band2+'_4rebinned.img'
		if os.path.isfile(myobs[0][:-14]+'/out/acisf'+stem+'_'+band+'_bkgmap_total.fits') == True:
			backmap=myobs[0][:-14]+'/out/acisf'+stem+'_'+band+'_bkgmap_total.fits'
		else:
			backmap=myobs[0][:-14]+'/out/acisf'+stem+'_'+band+'_bkgmap_instr.fits'
		
		s.call('dmextract infile="'+imagemap+'[bin pos=circle('+str(ra[i])+'d,'+str(dec[i])+'d,'+str(r50*0.000277778)+'d)]" mode=h bkg="'+backmap+'[bin pos=circle('+str(ra[i])+'d,'+str(dec[i])+'d,'+str(r50*0.000277778)+'d)]" outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
		cts=s.check_output('dmlist "counts.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
		bkg=s.check_output('dmlist "counts.fits[cols BG_COUNTS]" data,clean | grep -v BG_COUNTS',shell=True)
		cts,bkg=float(cts),float(bkg)
		net=cts-bkg
		
		if net <= 0:
			s.call('dmextract infile="'+imagemap+'[bin pos=circle('+str(ra[i])+'d,'+str(dec[i])+'d,'+str(r90*0.000277778)+'d)]" mode=h bkg="'+backmap+'[bin pos=circle('+str(ra[i])+'d,'+str(dec[i])+'d,'+str(r90*0.000277778)+'d)]" outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
			cts=s.check_output('dmlist "counts.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
			bkg=s.check_output('dmlist "counts.fits[cols BG_COUNTS]" data,clean | grep -v BG_COUNTS',shell=True)
			cts,bkg=float(cts),float(bkg)
			net=cts-bkg
			if net <= 0:
				poserr.append(0.1)
			else:
				poserr.append(r90/np.sqrt(net))
		
		else:
			poserr.append(r50/np.sqrt(net))
		'''
	else:
		sigmas=[]
		for j in range(len(myobs)):
			obs=myobs[j][36:-30]
			if len(obs) == 4:
				stem='0'+obs
			elif len(obs) == 3:
				stem='00'+obs
			elif len(obs) == 5:
				stem=obs
			path=myobs[j][:-14]+'/out/acisf'+stem+'_broad_expomap.fits'
	
			s.call('punlearn dmcoords',shell=True)
			s.call('dmcoords '+path+' asol=none opt=cel celfmt=deg ra='+str(ra[i])+' dec='+str(dec[i])+'',shell=True)
			theta=s.check_output('pget dmcoords theta',shell=True)
			#thetas.append([obs,float(theta)])
		
			diff=abs(oa_grid-float(theta))
			r50_new=r50_grid[diff==min(diff)]
			r50_new=r50_new[0]
			
			imagemap=wd+'data/'+obs+'/repro_new_asol/acisf'+stem+'_repro_'+band2+'_4rebinned.img'
			if os.path.isfile(wd+'data/'+obs+'/repro_new_asol/out/acisf'+stem+'_'+band+'_bkgmap_total.fits') == True:
				backmap=wd+'data/'+obs+'/repro_new_asol/out/acisf'+stem+'_'+band+'_bkgmap_total.fits'
			else:
				backmap=wd+'data/'+obs+'/repro_new_asol/out/acisf'+stem+'_'+band+'_bkgmap_instr.fits'
		
			s.call('dmextract infile="'+imagemap+'[bin pos=circle('+str(ra[i])+'d,'+str(dec[i])+'d,'+str(r50_new*0.000277778)+'d)]" mode=h bkg="'+backmap+'[bin pos=circle('+str(ra[i])+'d,'+str(dec[i])+'d,'+str(r50_new*0.000277778)+'d)]" outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
			cts=s.check_output('dmlist "counts.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
			bkg=s.check_output('dmlist "counts.fits[cols BG_COUNTS]" data,clean | grep -v BG_COUNTS',shell=True)
			cts,bkg=float(cts),float(bkg)
			net=cts-bkg
			if net <= 0:
				sigmas.append(10000.)
			else:
				sigmas.append(r50_new/np.sqrt(net))
		
		sigmas=np.array(sigmas)
		inv=1./sigmas**2
		poserr.append(1./np.sqrt(np.sum(inv)))
		
		'''
		thetas_ar=np.array(thetas)
		thetas_sorted = thetas_ar[thetas_ar[:,1].argsort()]
		
		obsid=thetas_sorted[0][0]
		if len(obsid) == 4:
			stem='0'+obsid
		elif len(obsid) == 3:
			stem='00'+obsid
		elif len(obsid) == 5:
			stem=obsid
			
		theta=float(thetas_sorted[0][1])
		
		#if band != 'hard':
		#	r90=1+10*(theta/10.)**2
		#else:
		#	r90=1.8+10*(theta/10.)**2
		#r50=r90*5./9.
		'''
		
		'''
		imagemap=wd+'data/'+obsid+'/repro_new_asol/acisf'+stem+'_repro_'+band2+'_4rebinned.img'
		if os.path.isfile(wd+'data/'+obsid+'/repro_new_asol/out/acisf'+stem+'_'+band+'_bkgmap_total.fits') == True:
			backmap=wd+'data/'+obsid+'/repro_new_asol/out/acisf'+stem+'_'+band+'_bkgmap_total.fits'
		else:
			backmap=wd+'data/'+obsid+'/repro_new_asol/out/acisf'+stem+'_'+band+'_bkgmap_instr.fits'
		
		s.call('dmextract infile="'+imagemap+'[bin pos=circle('+str(ra[i])+'d,'+str(dec[i])+'d,'+str(r50*0.000277778)+'d)]" mode=h bkg="'+backmap+'[bin pos=circle('+str(ra[i])+'d,'+str(dec[i])+'d,'+str(r50*0.000277778)+'d)]" outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
		cts=s.check_output('dmlist "counts.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
		bkg=s.check_output('dmlist "counts.fits[cols BG_COUNTS]" data,clean | grep -v BG_COUNTS',shell=True)
		cts,bkg=float(cts),float(bkg)
		net=cts-bkg
		
		if net <= 0:
			s.call('dmextract infile="'+imagemap+'[bin pos=circle('+str(ra[i])+'d,'+str(dec[i])+'d,'+str(r90*0.000277778)+'d)]" mode=h bkg="'+backmap+'[bin pos=circle('+str(ra[i])+'d,'+str(dec[i])+'d,'+str(r90*0.000277778)+'d)]" outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
			cts=s.check_output('dmlist "counts.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
			bkg=s.check_output('dmlist "counts.fits[cols BG_COUNTS]" data,clean | grep -v BG_COUNTS',shell=True)
			cts,bkg=float(cts),float(bkg)
			net=cts-bkg
			if net <= 0:
				poserr.append(0.1)
			else:
				poserr.append(r90/np.sqrt(net))
		else:
			poserr.append(r50/np.sqrt(net))
		'''
		

	'''
	if net <= 0:
		s.call('dmextract infile="'+imagemap+'[bin pos=circle('+str(ra[i])+'d,'+str(dec[i])+'d,'+str(r90*0.000277778)+'d)]" mode=h bkg="'+backmap+'[bin pos=circle('+str(ra[i])+'d,'+str(dec[i])+'d,'+str(r90*0.000277778)+'d)]" outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
		cts=s.check_output('dmlist "counts.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
		bkg=s.check_output('dmlist "counts.fits[cols BG_COUNTS]" data,clean | grep -v BG_COUNTS',shell=True)
		cts,bkg=float(cts),float(bkg)
		net=cts-bkg
		if net <= 0:
			poserr.append(0.1)
		else:
			poserr.append(r90/np.sqrt(net))
	
	else:
		poserr.append(r50/np.sqrt(net))
	'''
poserr=np.array(poserr)
poserr[poserr<0.1]=0.1
print(used_r90,'times used R90')
bins=np.logspace(np.log10(min(poserr)),np.log10(max(poserr)),100)

plt.figure()
plt.hist(poserr,bins=bins)
plt.xscale('log')
plt.show()

w=open(wd+'cdwfs_poserr_r50.dat','w')
for i in range(len(poserr)):
	w.write(str(poserr[i])+'\n')
w.close()
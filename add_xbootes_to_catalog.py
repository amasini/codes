import numpy as np
import sys
import subprocess as s
import scipy.stats.distributions
import time
import os
from astropy.io import fits
from astropy.table import Table

wd = '/Users/alberto/Desktop/XBOOTES/'

#### Part 0: WRITE OUT THE XBOOTES SOURCES WE MISS AND THE ONES WITH MULTIPLE CDWFS COUNTERPART
'''
kcat=fits.open(wd+'xbootes_kenter+05.fits')
ra_k=kcat[1].data['RAJ2000']
dec_k=kcat[1].data['DEJ2000']
name_k=kcat[1].data['CXOXB']
kcat.close()

cat=fits.open(wd+'new_mosaics_detection/cdwfs_merged_cat1_exp-psf.fits')
xb_count = cat[1].data['XB_ID']
cat.close()

xb_count0 = xb_count[xb_count != '0']
print(len(xb_count0),'total XB counterparts')
xb_unique = np.unique(xb_count0)
print(len(xb_unique),'unique XB counterparts')

w=open(wd+'xbootes_missing.dat','w')
w.write('CXOXB \t RA \t DEC\n')
for i in range(len(ra_k)):	
	if len(xb_count[xb_count == name_k[i]]) == 0:
		w.write(str(name_k[i])+' \t '+str(ra_k[i])+' \t '+str(dec_k[i])+'\n')
w.close()
sys.exit()
'''

#### Part 1: filter XBOOTES to create list of sources to be added to the catalog
'''
# This file contains the N XBOOTES sources missing from CDWFS
rak,deck=np.genfromtxt(wd+'xbootes_missing.dat',unpack=True, skip_header=1,usecols=[1,2])
namek=np.genfromtxt(wd+'xbootes_missing.dat',unpack=True, skip_header=1,usecols=0, dtype='str')

# Probability thresholds
cutf, cuts, cuth = 10**(-4.63),10**(-4.57),10**(-4.40)

w=open(wd+'xbootes_tbadded.dat','w')
w.write('CXOXB \t RA \t DEC \t PROB_F \t R90_F \t TOT_F \t BKG_F \t PROB_S \t R90_S \t TOT_S \t BKG_S \t PROB_H \t R90_H \t TOT_H \t BKG_H\n')
keep, used_r90=0,0
tin=time.time()

for i in range(len(rak)):
	print(i+1,len(rak))
	
	r90,prob,tot,back=[],[],[],[]
	for band in ['broad','soft','hard']:
			
		path=wd+'new_mosaics_detection/cdwfs_'+band+'_r90_4reb.fits'
		ima=fits.open(path)
		im=ima[0].data
		ima.close()
	
		imagemap=wd+'new_mosaics_detection/cdwfs_'+band+'_4reb.fits'
		backmap=wd+'new_mosaics_detection/cdwfs_'+band+'_bkgmap_4reb.fits'
		
		# Extract R90
		s.call('dmcoords '+path+' asol=none opt=cel celfmt=deg ra='+str(rak[i])+' dec='+str(deck[i])+'',shell=True)
		res=s.check_output('pget dmcoords logicalx logicaly',shell=True)
		logicalx,logicaly=res.splitlines()
		logicalx,logicaly=float(logicalx),float(logicaly)
		av_r90=im[int(round(logicaly)-1),int(round(logicalx)-1)]
		
		#extract counts and background using r90
		s.call('dmextract infile="'+imagemap+'[bin pos=circle('+str(rak[i])+'d,'+str(deck[i])+'d,'+str(av_r90*0.000277778)+'d)]" mode=h bkg="'+backmap+'[bin pos=circle('+str(rak[i])+'d,'+str(deck[i])+'d,'+str(av_r90*0.000277778)+'d)]" outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
		cts=s.check_output('dmlist "counts.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
		bkg=s.check_output('dmlist "counts.fits[cols BG_COUNTS]" data,clean | grep -v BG_COUNTS',shell=True)

		cts,bkg=float(cts),float(bkg)
		
		r90.append(av_r90)
		prob.append(scipy.stats.distributions.poisson.pmf(cts,bkg))
		tot.append(cts)
		back.append(bkg)
		
	# If the source is above threshold in at least one band,
	if (prob[0]<=cutf) | (prob[1]<=cuts) | (prob[2]<=cuth):
		keep=keep+1		
		
		w.write(str(namek[i])+' \t '+str(rak[i])+' \t '+str(deck[i])+' \t '+str(prob[0])+' \t '+str(r90[0])+' \t '+str(tot[0])+' \t '+str(back[0])+' \t '+str(prob[1])+' \t '+str(r90[1])+' \t '+str(tot[1])+' \t '+str(back[1])+' \t '+str(prob[2])+' \t '+str(r90[2])+' \t '+str(tot[2])+' \t '+str(back[2])+'\n')

w.close()
print(keep,'XBOOTES sources which satisfy reliability cut.')
print((time.time()-tin)/60.,'minutes for the match.')
sys.exit()
'''


#### Part 2: compute all the quantities (with errors and upperlimits)

namek,rak,deck,probf,r90f,totf,bkgf,probs,r90s,tots,bkgs,probh,r90h,toth,bkgh=np.genfromtxt(wd+'xbootes_tbadded.dat', skip_header=1, unpack=True)
namek=np.genfromtxt(wd+'xbootes_tbadded.dat',unpack=True, skip_header=1,usecols=0, dtype='str')

# Probability thresholds
cutf, cuts, cuth = 10**(-4.63),10**(-4.57),10**(-4.40)

cat0 = fits.open(wd+'new_mosaics_detection/cdwfs_merged_cat1_exp-psf.fits')
howmany = len(cat0[1].data['RA'])+1
cat0.close()

id,ra_u,dec_u,poserr,netf,e_netf_up,e_netf_lo,expf,crf,e_crf_up,e_crf_lo,fluxf,e_fluxf_up,e_fluxf_lo,nets,e_nets_up,e_nets_lo,exps,crs,e_crs_up,e_crs_lo,fluxs,e_fluxs_up,e_fluxs_lo,neth,e_neth_up,e_neth_lo,exph,crh,e_crh_up,e_crh_lo,fluxh,e_fluxh_up,e_fluxh_lo,hr,e_hr_up,e_hr_lo=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
used_r90=0
for i in range(len(rak)):
	print(i+1,len(rak))
	
	id.append(i+howmany)
	ra_u.append(rak[i])
	dec_u.append(deck[i])
	
	prob_array=[probf[i],probs[i],probh[i]]
	
	# Choose the POSERR based on most significant band
	if min(prob_array) == probf[i]: 
		band = 'broad'
		r50=r90f[i]*5./9.
		r90band=r90f[i]
		netband=totf[i]-bkgf[i]
	elif min(prob_array) == probs[i]:
		band = 'soft'
		r50=r90s[i]*5./9.
		r90band=r90s[i]
		netband=tots[i]-bkgs[i]
	else:
		band = 'hard'
		r50=r90h[i]*5./9.
		r90band=r90h[i]
		netband=toth[i]-bkgh[i]
	
	imagemap=wd+'new_mosaics_detection/cdwfs_'+band+'_4reb.fits'
	backmap=wd+'new_mosaics_detection/cdwfs_'+band+'_bkgmap_4reb.fits'
	
	s.call('dmextract infile="'+imagemap+'[bin pos=circle('+str(rak[i])+'d,'+str(deck[i])+'d,'+str(r50*0.000277778)+'d)]" mode=h bkg="'+backmap+'[bin pos=circle('+str(rak[i])+'d,'+str(deck[i])+'d,'+str(r50*0.000277778)+'d)]" outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
	cts=s.check_output('dmlist "counts.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
	bkg=s.check_output('dmlist "counts.fits[cols BG_COUNTS]" data,clean | grep -v BG_COUNTS',shell=True)
	cts,bkg=float(cts),float(bkg)
	net=cts-bkg
	if net > 0:
		poserr.append(r50/np.sqrt(net))
	else:
		used_r90 += 1
			
		if netband > 0:
			poserr.append(r90band/np.sqrt(netband))
		else:
			print('problem here:',rak[i],deck[i],r90f[i],band,totf[i],bkgf[i],netf[i])
			sys.exit()	
	
	
	#extract exposure,src and bkg counts BROAD BAND
	expomap=wd+'new_mosaics_detection/cdwfs_broad_expomap_4reb.fits'

	s.call('dmextract infile="'+expomap+'[bin pos=circle('+str(rak[i])+'d,'+str(deck[i])+'d,'+str(r90f[i]*0.000277778)+'d)]" mode=h outfile=expo2.fits opt=generic mode=h clobber=yes',shell=True)
	s.call('dmlist "expo2.fits[cols COUNTS, AREA]" data,clean | grep -v COUNTS > expo2.dat',shell=True)
	(totexpo,npix)=np.genfromtxt('expo2.dat',unpack=True)
	av_exp=16.0*totexpo/npix # npix is the number of NATIVE CHANDRA pixels, so need to divide it by 16!
	
	path=wd+'new_mosaics_detection/cdwfs_broad_ecfmap_4reb.fits'
	ima=fits.open(path)
	im2=ima[0].data
	ima.close()
	
	s.call('dmcoords '+path+' asol=none opt=cel celfmt=deg ra='+str(rak[i])+' dec='+str(deck[i])+'',shell=True)
	res=s.check_output('pget dmcoords logicalx logicaly',shell=True)
	logicalx,logicaly=res.splitlines()
	logicalx,logicaly=float(logicalx),float(logicaly)

	av_ecf=im2[int(round(logicaly)-1),int(round(logicalx)-1)]
	
	if probf[i] <= cutf: # Source is significant in F band
		
		e_cts_p=1+np.sqrt(totf[i]+0.75) # Gehrels+86 1sigma errors
		e_cts_n=np.sqrt(totf[i]-0.25)
		e_bkg_p=1+np.sqrt(bkgf[i]+0.75)
		if bkg >= 0.25:
			e_bkg_n=np.sqrt(bkgf[i]-0.25)
		else:
			e_bkg_n=0
		
		net = totf[i]-bkgf[i]
		e_net_p=np.sqrt(e_cts_p**2+e_bkg_p**2) # Propagate the errors
		e_net_n=np.sqrt(e_cts_n**2+e_bkg_n**2)

		cr=net*1.1/av_exp
		e_cr_p=e_net_p*1.1/av_exp # Propagate the errors
		e_cr_n=e_net_n*1.1/av_exp
		flux=cr*av_ecf
		e_flux_p=e_cr_p*av_ecf
		e_flux_n=e_cr_n*av_ecf
		
		netf.append(net)
		e_netf_up.append(e_net_p)
		e_netf_lo.append(e_net_n)
		expf.append(av_exp)
		crf.append(cr)
		e_crf_up.append(e_cr_p)
		e_crf_lo.append(e_cr_n)
		fluxf.append(flux)
		e_fluxf_up.append(e_flux_p)
		e_fluxf_lo.append(e_flux_n)
		
	else:
		
		netff = totf[i]-bkgf[i]
		netf2=netff+3*np.sqrt(netff+1)+(11./3.) #3sigma upper limit on net counts
		bkgf2=bkgf[i]+3*np.sqrt(bkgf[i]+1)+(11./3.) #3sigma upper limit on background counts
		
		net=np.sqrt(netf2**2+bkgf2**2) #propagate to get final 3sigma (hope it's correct) - check with Ryan
		
		cr=net*1.1/av_exp

		flux=cr*av_ecf
		
		netf.append(net)
		e_netf_up.append(0)
		e_netf_lo.append(0)
		expf.append(av_exp)
		crf.append(cr)
		e_crf_up.append(0)
		e_crf_lo.append(0)
		fluxf.append(flux)
		e_fluxf_up.append(0)
		e_fluxf_lo.append(0)
		
		
	#extract exposure,src and bkg counts SOFT BAND
	expomap=wd+'new_mosaics_detection/cdwfs_soft_expomap_4reb.fits'

	s.call('dmextract infile="'+expomap+'[bin pos=circle('+str(rak[i])+'d,'+str(deck[i])+'d,'+str(r90s[i]*0.000277778)+'d)]" mode=h outfile=expo2.fits opt=generic mode=h clobber=yes',shell=True)
	s.call('dmlist "expo2.fits[cols COUNTS, AREA]" data,clean | grep -v COUNTS > expo2.dat',shell=True)
	(totexpo,npix)=np.genfromtxt('expo2.dat',unpack=True)
	av_exp=16.0*totexpo/npix # npix is the number of NATIVE CHANDRA pixels, so need to divide it by 16!
	
	path=wd+'new_mosaics_detection/cdwfs_soft_ecfmap_4reb.fits'
	ima=fits.open(path)
	im2=ima[0].data
	ima.close()
	
	s.call('dmcoords '+path+' asol=none opt=cel celfmt=deg ra='+str(rak[i])+' dec='+str(deck[i])+'',shell=True)
	res=s.check_output('pget dmcoords logicalx logicaly',shell=True)
	logicalx,logicaly=res.splitlines()
	logicalx,logicaly=float(logicalx),float(logicaly)

	av_ecf=im2[int(round(logicaly)-1),int(round(logicalx)-1)]
	
	if probs[i] <= cuts: # Source is significant in S band
		
		e_cts_p=1+np.sqrt(tots[i]+0.75) # Gehrels+86 1sigma errors
		e_cts_n=np.sqrt(tots[i]-0.25)
		e_bkg_p=1+np.sqrt(bkgs[i]+0.75)
		if bkg >= 0.25:
			e_bkg_n=np.sqrt(bkgs[i]-0.25)
		else:
			e_bkg_n=0
		
		net = tots[i]-bkgs[i]
		e_net_p=np.sqrt(e_cts_p**2+e_bkg_p**2) # Propagate the errors
		e_net_n=np.sqrt(e_cts_n**2+e_bkg_n**2)

		cr=net*1.1/av_exp
		e_cr_p=e_net_p*1.1/av_exp # Propagate the errors
		e_cr_n=e_net_n*1.1/av_exp
		flux=cr*av_ecf
		e_flux_p=e_cr_p*av_ecf
		e_flux_n=e_cr_n*av_ecf
		
		nets.append(net)
		e_nets_up.append(e_net_p)
		e_nets_lo.append(e_net_n)
		exps.append(av_exp)
		crs.append(cr)
		e_crs_up.append(e_cr_p)
		e_crs_lo.append(e_cr_n)
		fluxs.append(flux)
		e_fluxs_up.append(e_flux_p)
		e_fluxs_lo.append(e_flux_n)
		
	else:
		
		netff = tots[i]-bkgs[i]
		netf2=netff+3*np.sqrt(netff+1)+(11./3.) #3sigma upper limit on net counts
		bkgf2=bkgs[i]+3*np.sqrt(bkgs[i]+1)+(11./3.) #3sigma upper limit on background counts
		
		net=np.sqrt(netf2**2+bkgf2**2) #propagate to get final 3sigma (hope it's correct) - check with Ryan
		
		cr=net*1.1/av_exp

		flux=cr*av_ecf
		
		nets.append(net)
		e_nets_up.append(0)
		e_nets_lo.append(0)
		exps.append(av_exp)
		crs.append(cr)
		e_crs_up.append(0)
		e_crs_lo.append(0)
		fluxs.append(flux)
		e_fluxs_up.append(0)
		e_fluxs_lo.append(0)
	
	
	#extract exposure,src and bkg counts HARD BAND
	expomap=wd+'new_mosaics_detection/cdwfs_hard_expomap_4reb.fits'

	s.call('dmextract infile="'+expomap+'[bin pos=circle('+str(rak[i])+'d,'+str(deck[i])+'d,'+str(r90h[i]*0.000277778)+'d)]" mode=h outfile=expo2.fits opt=generic mode=h clobber=yes',shell=True)
	s.call('dmlist "expo2.fits[cols COUNTS, AREA]" data,clean | grep -v COUNTS > expo2.dat',shell=True)
	(totexpo,npix)=np.genfromtxt('expo2.dat',unpack=True)
	av_exp=16.0*totexpo/npix # npix is the number of NATIVE CHANDRA pixels, so need to divide it by 16!
	
	path=wd+'new_mosaics_detection/cdwfs_hard_ecfmap_4reb.fits'
	ima=fits.open(path)
	im2=ima[0].data
	ima.close()
	
	s.call('dmcoords '+path+' asol=none opt=cel celfmt=deg ra='+str(rak[i])+' dec='+str(deck[i])+'',shell=True)
	res=s.check_output('pget dmcoords logicalx logicaly',shell=True)
	logicalx,logicaly=res.splitlines()
	logicalx,logicaly=float(logicalx),float(logicaly)

	av_ecf=im2[int(round(logicaly)-1),int(round(logicalx)-1)]
	
	if probh[i] <= cuth: # Source is significant in H band
		
		e_cts_p=1+np.sqrt(toth[i]+0.75) # Gehrels+86 1sigma errors
		e_cts_n=np.sqrt(toth[i]-0.25)
		e_bkg_p=1+np.sqrt(bkgh[i]+0.75)
		if bkg >= 0.25:
			e_bkg_n=np.sqrt(bkgh[i]-0.25)
		else:
			e_bkg_n=0
		
		net = toth[i]-bkgh[i]
		e_net_p=np.sqrt(e_cts_p**2+e_bkg_p**2) # Propagate the errors
		e_net_n=np.sqrt(e_cts_n**2+e_bkg_n**2)

		cr=net*1.1/av_exp
		e_cr_p=e_net_p*1.1/av_exp # Propagate the errors
		e_cr_n=e_net_n*1.1/av_exp
		flux=cr*av_ecf
		e_flux_p=e_cr_p*av_ecf
		e_flux_n=e_cr_n*av_ecf
		
		neth.append(net)
		e_neth_up.append(e_net_p)
		e_neth_lo.append(e_net_n)
		exph.append(av_exp)
		crh.append(cr)
		e_crh_up.append(e_cr_p)
		e_crh_lo.append(e_cr_n)
		fluxh.append(flux)
		e_fluxh_up.append(e_flux_p)
		e_fluxh_lo.append(e_flux_n)
		
	else:
		
		netff = toth[i]-bkgh[i]
		netf2=netff+3*np.sqrt(netff+1)+(11./3.) #3sigma upper limit on net counts
		bkgf2=bkgh[i]+3*np.sqrt(bkgh[i]+1)+(11./3.) #3sigma upper limit on background counts
		
		net=np.sqrt(netf2**2+bkgf2**2) #propagate to get final 3sigma (hope it's correct) - check with Ryan
		
		cr=net*1.1/av_exp

		flux=cr*av_ecf
		
		neth.append(net)
		e_neth_up.append(0)
		e_neth_lo.append(0)
		exph.append(av_exp)
		crh.append(cr)
		e_crh_up.append(0)
		e_crh_lo.append(0)
		fluxh.append(flux)
		e_fluxh_up.append(0)
		e_fluxh_lo.append(0)
		
	r90ratio=(r90h[i]/r90s[i])**2 # This takes into account that R90 is larger in hard band into the HR computation
	# Compute HR with BEHR
	os.chdir('/Users/alberto/BEHR/') 
	behr='./BEHR softsrc='+str(int(tots[i]))+' hardsrc='+str(int(toth[i]))+' softbkg='+str(int(round(bkgs[i])))+' hardbkg='+str(int(round(bkgh[i])))+' softarea=1 hardarea='+str(r90ratio)+' output='+str(i+howmany)+' outputHR=true'
	s.call(behr,shell=True)
	
	a,b,c=np.genfromtxt('/Users/alberto/BEHR/'+str(i+howmany)+'_HR.txt',unpack=True,usecols=[2,3,4]) # Use median
	hr.append(a)
	e_hr_lo.append(a-b)
	e_hr_up.append(c-a)
		
#write catalog
cat=Table([id,ra_u,dec_u,poserr,probf,r90f,totf,bkgf,netf,e_netf_up,e_netf_lo,expf,crf,e_crf_up,e_crf_lo,fluxf,e_fluxf_up,e_fluxf_lo,probs,r90s,tots,bkgs,nets,e_nets_up,e_nets_lo,exps,crs,e_crs_up,e_crs_lo,fluxs,e_fluxs_up,e_fluxs_lo,probh,r90h,toth,bkgh,neth,e_neth_up,e_neth_lo,exph,crh,e_crh_up,e_crh_lo,fluxh,e_fluxh_up,e_fluxh_lo,hr,e_hr_up,e_hr_lo,namek],names=('ID','RA','DEC','POS_ERR','PROB_F','R90_F','TOT_F','BKG_F','NET_F','E_NET_F_+','E_NET_F_-','EXP_F','CR_F','E_CR_F_+','E_CR_F_-','FLUX_F','E_FLUX_F_+','E_FLUX_F_-','PROB_S','R90_S','TOT_S','BKG_S','NET_S','E_NET_S_+','E_NET_S_-','EXP_S','CR_S','E_CR_S_+','E_CR_S_-','FLUX_S','E_FLUX_S_+','E_FLUX_S_-','PROB_H','R90_H','TOT_H','BKG_H','NET_H','E_NET_H_+','E_NET_H_-','EXP_H','CR_H','E_CR_H_+','E_CR_H_-','FLUX_H','E_FLUX_H_+','E_FLUX_H_-','HR','E_HR_+','E_HR_-','XB_ID'))
cat.write(wd+'xbootes_tbadded.fits',format='fits',overwrite=True)

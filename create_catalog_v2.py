# This script takes the merged F,S,H catalogs and applies probability cuts based on simulations, computes HR, and matches with Kenter+05 catalog.
import numpy as np
from astropy.io import fits
from astropy.table import Table
import sys
import matplotlib.pyplot as plt
import subprocess as s
import os
import scipy
import scipy.stats.distributions

def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*xx*3600)**2 +((pointa[1]-pointb[1])*3600)**2)

wd='/Users/alberto/Desktop/XBOOTES/'
date = '200113'

cat=fits.open(wd+'new_mosaics_detection/cdwfs_merged_cat0_'+date+'.fits')
raf=cat[1].data['RA_F']
decf=cat[1].data['DEC_F']
probf=cat[1].data['PROB_F']
r90f=cat[1].data['R90_F']
totf=cat[1].data['TOT_F']
bkgf=cat[1].data['BKG_F']
netf=cat[1].data['NET_F']
e_netf_up=cat[1].data['E_NET_F_+']
e_netf_lo=cat[1].data['E_NET_F_-']
expf=cat[1].data['EXP_F']
crf=cat[1].data['CR_F']
e_crf_up=cat[1].data['E_CR_F_+']
e_crf_lo=cat[1].data['E_CR_F_-']
fluxf=cat[1].data['FLUX_F']
e_fluxf_up=cat[1].data['E_FLUX_F_+']
e_fluxf_lo=cat[1].data['E_FLUX_F_-']

ras=cat[1].data['RA_S']
decs=cat[1].data['DEC_S']
probs=cat[1].data['PROB_S']
r90s=cat[1].data['R90_S']
tots=cat[1].data['TOT_S']
bkgs=cat[1].data['BKG_S']
nets=cat[1].data['NET_S']
e_nets_up=cat[1].data['E_NET_S_+']
e_nets_lo=cat[1].data['E_NET_S_-']
exps=cat[1].data['EXP_S']
crs=cat[1].data['CR_S']
e_crs_up=cat[1].data['E_CR_S_+']
e_crs_lo=cat[1].data['E_CR_S_-']
fluxs=cat[1].data['FLUX_S']
e_fluxs_up=cat[1].data['E_FLUX_S_+']
e_fluxs_lo=cat[1].data['E_FLUX_S_-']

rah=cat[1].data['RA_H']
dech=cat[1].data['DEC_H']
probh=cat[1].data['PROB_H']
r90h=cat[1].data['R90_H']
toth=cat[1].data['TOT_H']
bkgh=cat[1].data['BKG_H']
neth=cat[1].data['NET_H']
e_neth_up=cat[1].data['E_NET_H_+']
e_neth_lo=cat[1].data['E_NET_H_-']
exph=cat[1].data['EXP_H']
crh=cat[1].data['CR_H']
e_crh_up=cat[1].data['E_CR_H_+']
e_crh_lo=cat[1].data['E_CR_H_-']
fluxh=cat[1].data['FLUX_H']
e_fluxh_up=cat[1].data['E_FLUX_H_+']
e_fluxh_lo=cat[1].data['E_FLUX_H_-']

# Probability cuts in F,S,H bands
cutf,cuts,cuth=10**(-4.63),10**(-4.57),10**(-4.40)

# Undetected sources have probability of being spurious 9999
probf[raf==0.0]=9999
probs[ras==0.0]=9999
probh[rah==0.0]=9999

print(len(raf[probf<=cutf]))
print(len(ras[probs<=cuts]))
print(len(rah[probh<=cuth]))
print('we should have a total of',len(raf[(probf<=cutf) | (probs<=cuts) | (probh<=cuth)]),'sources.')

# Consider only those sources that are to be included in the final catalog
my_sample = cat[1].data[(probf<=cutf) | (probs<=cuts) | (probh<=cuth)]

# Choose best coordinates based on most significant detection
probs = np.array([my_sample['PROB_F'],my_sample['PROB_S'],my_sample['PROB_H']])
index = np.argmin(probs.T, axis=1)

my_raf = my_sample['RA_F']
my_ras = my_sample['RA_S']
my_rah = my_sample['RA_H']
my_decf = my_sample['DEC_F']
my_decs = my_sample['DEC_S']
my_dech = my_sample['DEC_H']

ra_u = my_raf.copy()
ra_u[index==1] = my_ras[index==1]
ra_u[index==2] = my_rah[index==2]

dec_u = my_decf.copy()
dec_u[index==1] = my_decs[index==1]
dec_u[index==2] = my_dech[index==2]

# Positional errors - use appropriate R90 (F for index=0, S for index=1 and H for index=2)
my_r90f = my_sample['R90_F']
my_r90s = my_sample['R90_S']
my_r90h = my_sample['R90_H']

# Compute R50 from R90
r90 = my_r90f.copy()
r90[index == 1] = my_r90s[index==1]
r90[index == 2] = my_r90h[index==2]

my_netf = my_sample['NET_F']
my_nets = my_sample['NET_S']
my_neth = my_sample['NET_H']
net_u = my_netf.copy()
net_u[index == 1] = my_nets[index==1]
net_u[index == 2] = my_neth[index==2]

poserr_r90 = r90/np.sqrt(net_u)

r50 = r90*5./9.

# dmextract has to be called on the right band (F for index=0, S for index=1 and H for index=2)
w0=open(wd+'new_mosaics_detection/cdwfs_broad_src_'+date+'_forR50.reg','w')
w1=open(wd+'new_mosaics_detection/cdwfs_soft_src_'+date+'_forR50.reg','w')
w2=open(wd+'new_mosaics_detection/cdwfs_hard_src_'+date+'_forR50.reg','w')
for i in range(len(r50)):
	if index[i] == 0:
		w0.write('circle('+str(ra_u[i])+'d, '+str(dec_u[i])+'d, '+str(r50[i])+'\")\n')
	elif index[i] == 1:
		w1.write('circle('+str(ra_u[i])+'d, '+str(dec_u[i])+'d, '+str(r50[i])+'\")\n')
	else:
		w2.write('circle('+str(ra_u[i])+'d, '+str(dec_u[i])+'d, '+str(r50[i])+'\")\n')
w0.close()
w1.close()
w2.close()

# Then, run dmextract to get the average R90 and average ECF
for band in ['broad','soft','hard']:
	s.call('dmextract \''+wd+'new_mosaics_detection/cdwfs_'+band+'_4reb.fits[bin sky=@'+wd+'new_mosaics_detection/cdwfs_'+band+'_src_'+date+'_forR50.reg]\' outfile='+wd+'new_mosaics_detection/cdwfs_'+band+'_src_'+date+'_R50AP.fits bkg=\''+wd+'new_mosaics_detection/cdwfs_'+band+'_bkgmap_4reb.fits[bin sky=@'+wd+'new_mosaics_detection/cdwfs_'+band+'_src_'+date+'_forR50.reg]\' clobber=yes', shell=True)

	cat = fits.open(wd+'new_mosaics_detection/cdwfs_'+band+'_src_'+date+'_R50AP.fits')
	tot = cat[1].data['COUNTS']
	bkg = cat[1].data['BG_COUNTS']
	cat.close()
	
	net = tot-bkg
	if band == 'broad':
		net0 = net
	elif band == 'soft':
		net1 = net
	else:
		net2 = net

poserr=np.zeros(len(r50))
poserr[index == 0] = r50[index == 0]/np.sqrt(net0)
poserr[index == 1] = r50[index == 1]/np.sqrt(net1)
poserr[index == 2] = r50[index == 2]/np.sqrt(net2)

print(len(poserr[np.isnan(poserr)==True]),'NaNs')
print(len(poserr[np.isinf(poserr)==True]),'Infs')
poserr[np.isnan(poserr)==True] = poserr_r90[np.isnan(poserr)==True]
poserr[np.isinf(poserr)==True] = poserr_r90[np.isinf(poserr)==True]

# Aperture photometry on undetected sources, with correct treatment of new detections and UL
my_probf = my_sample['PROB_F']
my_probs = my_sample['PROB_S']
my_probh = my_sample['PROB_H']
print(len(my_raf[my_raf==0.0]),'Non-detections in the F band')
print(len(my_ras[my_ras==0.0]),'Non-detections in the S band')
print(len(my_rah[my_rah==0.0]),'Non-detections in the H band')

sys.exit()



# Write out the region file with circles and 2" radius to get R90 and ECF
w=open(wd+'cdwfs_'+band+'_src_'+date+'_forAP.reg','w')
for i in range(len(ra)):
	w.write('circle('+str(ra[i])+'d, '+str(dec[i])+'d, '+str(r50[i])+'\")\n')
w.close()

count=0
ra_u,dec_u=[],[]
hr,e_hr_up,e_hr_lo=[],[],[]
out_name,unique,id,poserr=[],[],[],[]
os.chdir('/Users/alberto/BEHR/') 

kenter_match=0
kcat=fits.open(wd+'xbootes_kenter+05.fits')
ra_k=kcat[1].data['RAJ2000']
dec_k=kcat[1].data['DEJ2000']
name_k=kcat[1].data['CXOXB']
used_r90=0
for j in range(len(probf)):
	if (probf[j] <= cutf or probs[j] <= cuts or probh[j] <= cuth): # Source is above threshold in at least one band
		count=count+1
		id.append(count)
		unique.append(True)
		prob_array=[probf[j],probs[j],probh[j]]
		
		# Fix coordinates of the source, and compute POS_ERR
		if min(prob_array) == probf[j]: # Choose the ra,dec based on most significant band
			ra_u.append(raf[j])
			dec_u.append(decf[j])
			
			#compute positional error here as R50/sqrt(Net_R50)
			band='broad'
			imagemap=wd+'new_mosaics_detection/cdwfs_'+band+'_4reb.fits'
			backmap=wd+'new_mosaics_detection/cdwfs_'+band+'_bkgmap_4reb.fits'
			r50=r90f[j]*5./9.
			s.call('dmextract infile="'+imagemap+'[bin pos=circle('+str(raf[j])+'d,'+str(decf[j])+'d,'+str(r50*0.000277778)+'d)]" mode=h bkg="'+backmap+'[bin pos=circle('+str(raf[j])+'d,'+str(decf[j])+'d,'+str(r50*0.000277778)+'d)]" outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
			cts=s.check_output('dmlist "counts.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
			bkg=s.check_output('dmlist "counts.fits[cols BG_COUNTS]" data,clean | grep -v BG_COUNTS',shell=True)
			cts,bkg=float(cts),float(bkg)
			net=cts-bkg
			if net > 0:
				poserr.append(r50/np.sqrt(net))
			else:
				used_r90 += 1
				if netf[j] > 0:
					poserr.append(r90f[j]/np.sqrt(netf[j]))
				else:
					print('problem here:',raf[j],decf[j],r90,band,cts,bkg,net)
					sys.exit()
			
		elif min(prob_array) == probs[j]:
			ra_u.append(ras[j])
			dec_u.append(decs[j])
			
			#compute positional error here
			band='soft'
			imagemap=wd+'new_mosaics_detection/cdwfs_'+band+'_4reb.fits'
			backmap=wd+'new_mosaics_detection/cdwfs_'+band+'_bkgmap_4reb.fits'
			r50=r90s[j]*5./9.
			s.call('dmextract infile="'+imagemap+'[bin pos=circle('+str(ras[j])+'d,'+str(decs[j])+'d,'+str(r50*0.000277778)+'d)]" mode=h bkg="'+backmap+'[bin pos=circle('+str(ras[j])+'d,'+str(decs[j])+'d,'+str(r50*0.000277778)+'d)]" outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
			cts=s.check_output('dmlist "counts.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
			bkg=s.check_output('dmlist "counts.fits[cols BG_COUNTS]" data,clean | grep -v BG_COUNTS',shell=True)
			cts,bkg=float(cts),float(bkg)
			net=cts-bkg
			if net > 0:
				poserr.append(r50/np.sqrt(net))
			else:
				used_r90 += 1
				if nets[j] > 0:
					poserr.append(r90s[j]/np.sqrt(nets[j]))
				else:
					print('problem here:',ras[j],decs[j],r90,band,cts,bkg,net)
					sys.exit()
			
		else:
			ra_u.append(rah[j])
			dec_u.append(dech[j])
			
			#compute positional error here
			band='hard'
			imagemap=wd+'new_mosaics_detection/cdwfs_'+band+'_4reb.fits'
			backmap=wd+'new_mosaics_detection/cdwfs_'+band+'_bkgmap_4reb.fits'
			r50=r90h[j]*5./9.
			s.call('dmextract infile="'+imagemap+'[bin pos=circle('+str(rah[j])+'d,'+str(dech[j])+'d,'+str(r50*0.000277778)+'d)]" mode=h bkg="'+backmap+'[bin pos=circle('+str(rah[j])+'d,'+str(dech[j])+'d,'+str(r50*0.000277778)+'d)]" outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
			cts=s.check_output('dmlist "counts.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
			bkg=s.check_output('dmlist "counts.fits[cols BG_COUNTS]" data,clean | grep -v BG_COUNTS',shell=True)
			cts,bkg=float(cts),float(bkg)
			net=cts-bkg
			if net > 0:
				poserr.append(r50/np.sqrt(net))
			else:
				used_r90 += 1
				if neth[j] > 0:
					poserr.append(r90h[j]/np.sqrt(neth[j]))
				else:
					print('problem here:',rah[j],dech[j],r90,band,cts,bkg,net)
					sys.exit()

totf_u = totf[(probf[j] <= cutf) | (probs[j] <= cuts) | (probh[j] <= cuth)]	
ra_u = np.array(ra_u)
# For non-detections: compute aperture photometry and decide which are significant and which not
print(len(ra_u[totf_u == 0.0]))
sys.exit()
		
if totf[j] == 0.0: # if source is not detected in F band, extract aperture photometry 
	band='broad'
	
	path=wd+'new_mosaics_detection/cdwfs_'+band+'_r90_4reb.fits'
	s.call('dmcoords '+path+' asol=none opt=cel celfmt=deg ra='+str(ra_u[-1])+' dec='+str(dec_u[-1])+'',shell=True)
	res=s.check_output('pget dmcoords logicalx logicaly',shell=True)
	logicalx,logicaly=res.splitlines()
	logicalx,logicaly=float(logicalx),float(logicaly)
	
	ima=fits.open(path)
	im=ima[0].data
	ima.close()
	
	ima2=fits.open(wd+'new_mosaics_detection/cdwfs_'+band+'_ecfmap_4reb.fits')
	im2=ima2[0].data
	ima2.close()
	
	imagemap=wd+'new_mosaics_detection/cdwfs_'+band+'_4reb.fits'
	backmap=wd+'new_mosaics_detection/cdwfs_'+band+'_bkgmap_4reb.fits'
	expomap=wd+'new_mosaics_detection/cdwfs_'+band+'_expomap_4reb.fits'
	# if I have r90s, use that; otherwise, compute fresh new r90...
	if r90s[j] != 0.0:
		r90=r90s[j]
	else:
		#compute r90 from the psfmap
		r90=im[int(round(logicaly)-1),int(round(logicalx)-1)]
	
	r90f[j]=r90
	s.call('dmextract infile="'+imagemap+'[bin pos=circle('+str(ra_u[-1])+'d,'+str(dec_u[-1])+'d,'+str(r90*0.000277778)+'d)]" mode=h bkg="'+backmap+'[bin pos=circle('+str(ra_u[-1])+'d,'+str(dec_u[-1])+'d,'+str(r90*0.000277778)+'d)]" outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
	cts=s.check_output('dmlist "counts.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
	bkg=s.check_output('dmlist "counts.fits[cols BG_COUNTS]" data,clean | grep -v BG_COUNTS',shell=True)
	totf[j],bkgf[j]=float(cts),float(bkg)
	probf[j]=scipy.stats.distributions.poisson.pmf(totf[j],bkgf[j])
	
	#extract exposure from vignetting-corrected expomap
	s.call('dmextract infile="'+expomap+'[bin pos=circle('+str(ra_u[-1])+'d,'+str(dec_u[-1])+'d,'+str(r90*0.000277778)+'d)]" mode=h outfile=expo.fits opt=generic mode=h clobber=yes',shell=True)
	s.call('dmlist "expo.fits[cols COUNTS, AREA]" data,clean | grep -v COUNTS > expo.dat',shell=True)
	(totexpo,npix)=np.genfromtxt('expo.dat',unpack=True)
	av_exp=16.0*totexpo/npix # npix is the number of NATIVE CHANDRA pixels, so need divide it by 16!
	expf[j]=av_exp
	
	netff=totf[j]-bkgf[j]
	netf2=netff+3*np.sqrt(netff+1)+(11./3.) #3sigma upper limit on net counts
	bkgf2=bkgf[j]+3*np.sqrt(bkgf[j]+1)+(11./3.) #3sigma upper limit on background counts
	
	netf[j]=np.sqrt(netf2**2+bkgf2**2) #propagate to get final 3sigma (hope it's correct) - check with Ryan
	
	#3sigma upper limit on CR
	crf[j]=1.1*netf[j]/expf[j] 
	cf=im2[int(round(logicaly)-1),int(round(logicalx)-1)]
	
	#3sigma upper limit on flux
	fluxf[j]=cf*crf[j] 

if tots[j] == 0.0: # if source is not detected in S band, extract 3sigma upper limit and define counts for the HR
	band='soft'
	
	path=wd+'new_mosaics_detection/cdwfs_broad_r90_4reb.fits'
	s.call('dmcoords '+path+' asol=none opt=cel celfmt=deg ra='+str(ra_u[-1])+' dec='+str(dec_u[-1])+'',shell=True)
	res=s.check_output('pget dmcoords logicalx logicaly',shell=True)
	logicalx,logicaly=res.splitlines()
	logicalx,logicaly=float(logicalx),float(logicaly)
	
	ima=fits.open(path)
	im=ima[0].data
	ima.close()
	
	ima2=fits.open(wd+'new_mosaics_detection/cdwfs_'+band+'_ecfmap_4reb.fits')
	im2=ima2[0].data
	ima2.close()
	
	imagemap=wd+'new_mosaics_detection/cdwfs_'+band+'_4reb.fits'
	backmap=wd+'new_mosaics_detection/cdwfs_'+band+'_bkgmap_4reb.fits'
	expomap=wd+'new_mosaics_detection/cdwfs_'+band+'_expomap_4reb.fits'
	# if I have r90f, use that; otherwise, compute fresh new r90...
	if r90f[j] != 0.0:
		r90=r90f[j]
	else:
		#compute r90 from the psfmap
		r90=im[int(round(logicaly)-1),int(round(logicalx)-1)]
	
	r90s[j]=r90
	s.call('dmextract infile="'+imagemap+'[bin pos=circle('+str(ra_u[-1])+'d,'+str(dec_u[-1])+'d,'+str(r90*0.000277778)+'d)]" mode=h bkg="'+backmap+'[bin pos=circle('+str(ra_u[-1])+'d,'+str(dec_u[-1])+'d,'+str(r90*0.000277778)+'d)]" outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
	cts=s.check_output('dmlist "counts.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
	bkg=s.check_output('dmlist "counts.fits[cols BG_COUNTS]" data,clean | grep -v BG_COUNTS',shell=True)
	tots[j],bkgs[j]=float(cts),float(bkg)
	probs[j]=scipy.stats.distributions.poisson.pmf(tots[j],bkgs[j])
	
	#extract exposure from vignetting-corrected expomap
	s.call('dmextract infile="'+expomap+'[bin pos=circle('+str(ra_u[-1])+'d,'+str(dec_u[-1])+'d,'+str(r90*0.000277778)+'d)]" mode=h outfile=expo.fits opt=generic mode=h clobber=yes',shell=True)
	s.call('dmlist "expo.fits[cols COUNTS, AREA]" data,clean | grep -v COUNTS > expo.dat',shell=True)
	(totexpo,npix)=np.genfromtxt('expo.dat',unpack=True)
	av_exp=totexpo/npix
	exps[j]=av_exp
	
	netss=tots[j]-bkgs[j]
	nets2=netss+3*np.sqrt(netss+1)+(11./3.) #3sigma upper limit on net counts
	bkgs2=bkgs[j]+3*np.sqrt(bkgs[j]+1)+(11./3.) #3sigma upper limit on background counts
	
	nets[j]=np.sqrt(nets2**2+bkgs2**2) #propagate to get final 3sigma (hope it's correct) - check with Ryan
	
	#3sigma upper limit on CR
	crs[j]=1.1*nets[j]/exps[j] 
	cf=im2[int(round(logicaly)-1),int(round(logicalx)-1)]
	
	#3sigma upper limit on flux
	fluxs[j]=cf*crs[j]
	
if toth[j] == 0.0: # if source is not detected in H band, extract 3sigma upper limit and define counts for the HR
	band='hard'
	
	path=wd+'new_mosaics_detection/cdwfs_'+band+'_r90_4reb.fits'
	s.call('dmcoords '+path+' asol=none opt=cel celfmt=deg ra='+str(ra_u[-1])+' dec='+str(dec_u[-1])+'',shell=True)
	res=s.check_output('pget dmcoords logicalx logicaly',shell=True)
	logicalx,logicaly=res.splitlines()
	logicalx,logicaly=float(logicalx),float(logicaly)
	
	ima=fits.open(path)
	im=ima[0].data
	ima.close()
	
	ima2=fits.open(wd+'new_mosaics_detection/cdwfs_'+band+'_ecfmap_4reb.fits')
	im2=ima2[0].data
	ima2.close()
	
	imagemap=wd+'new_mosaics_detection/cdwfs_'+band+'_4reb.fits'
	backmap=wd+'new_mosaics_detection/cdwfs_'+band+'_bkgmap_4reb.fits'
	expomap=wd+'new_mosaics_detection/cdwfs_'+band+'_expomap_4reb.fits'
	# Here I have to compute fresh new r90... CREATE THE HARD BAND PSFMAP!
	#compute r90 from the psfmap
	r90=im[int(round(logicaly)-1),int(round(logicalx)-1)]
	
	r90h[j]=r90
	s.call('dmextract infile="'+imagemap+'[bin pos=circle('+str(ra_u[-1])+'d,'+str(dec_u[-1])+'d,'+str(r90*0.000277778)+'d)]" mode=h bkg="'+backmap+'[bin pos=circle('+str(ra_u[-1])+'d,'+str(dec_u[-1])+'d,'+str(r90*0.000277778)+'d)]" outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
	cts=s.check_output('dmlist "counts.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
	bkg=s.check_output('dmlist "counts.fits[cols BG_COUNTS]" data,clean | grep -v BG_COUNTS',shell=True)
	toth[j],bkgh[j]=float(cts),float(bkg)
	probh[j]=scipy.stats.distributions.poisson.pmf(toth[j],bkgh[j])
	
	#extract exposure from vignetting-corrected expomap
	s.call('dmextract infile="'+expomap+'[bin pos=circle('+str(ra_u[-1])+'d,'+str(dec_u[-1])+'d,'+str(r90*0.000277778)+'d)]" mode=h outfile=expo.fits opt=generic mode=h clobber=yes',shell=True)
	s.call('dmlist "expo.fits[cols COUNTS, AREA]" data,clean | grep -v COUNTS > expo.dat',shell=True)
	(totexpo,npix)=np.genfromtxt('expo.dat',unpack=True)
	av_exp=totexpo/npix
	exph[j]=av_exp
	
	nethh=toth[j]-bkgh[j]
	neth2=nethh+3*np.sqrt(nethh+1)+(11./3.) #3sigma upper limit on net counts
	bkgh2=bkgh[j]+3*np.sqrt(bkgh[j]+1)+(11./3.) #3sigma upper limit on background counts
	
	neth[j]=np.sqrt(neth2**2+bkgh2**2) #propagate to get final 3sigma (hope it's correct) - check with Ryan
	
	#3sigma upper limit on CR
	crh[j]=1.1*neth[j]/exph[j] 
	cf=im2[int(round(logicaly)-1),int(round(logicalx)-1)]
	
	#3sigma upper limit on flux
	fluxh[j]=cf*crh[j]

r90ratio=(r90h[j]/r90s[j])**2 # This takes into account that R90 is larger in hard band into the HR computation
# Compute HR with BEHR
behr='./BEHR softsrc='+str(int(tots[j]))+' hardsrc='+str(int(toth[j]))+' softbkg='+str(int(round(bkgs[j])))+' hardbkg='+str(int(round(bkgh[j])))+' softarea=1 hardarea='+str(r90ratio)+' output='+str(count)+' outputHR=true'
s.call(behr,shell=True)

a,b,c=np.genfromtxt('/Users/alberto/BEHR/'+str(count)+'_HR.txt',unpack=True,usecols=[2,3,4]) # Use median
hr.append(a)
e_hr_lo.append(a-b)
e_hr_up.append(c-a)

## Match with Kenter+05
cdwfs=[ra_u[-1],dec_u[-1]]

found=False
ra_cont,dec_cont,d,name=[],[],[],[]
for k in range(len(ra_k)): # Look up into Kenter catalog
	kenter=[ra_k[k],dec_k[k]]
	if distance(cdwfs,kenter) < 1.1*r90_u[-1]:
		if found==False: # This is the first match
			kenter_match=kenter_match+1
			found=True
			ra_cont.append(ra_k[k])
			dec_cont.append(dec_k[k])
			d.append(distance(cdwfs,kenter))
			name.append(name_k[k])
		elif found==True:
			ra_cont.append(ra_k[k])
			dec_cont.append(dec_k[k])
			d.append(distance(cdwfs,kenter))
			name.append(name_k[k])
			print('Finding a double countepart for '+str(ra_u[-1])+', '+str(dec_u[-1])+'')
for l in range(len(ra_cont)): # retain the closest counterpart only
	if min(d)==d[l]:
		out_name.append(name[l])
if found==False: # No kenter counterpart was found
	out_name.append(0)

else:
	unique.append(False)
kcat.close()

##### r90 used 122/7234 times = ~1.7% (2.4% as of 3-Oct-19)

# Write final, clean catalog with HR and Murray_ID
#define output stuff
unique=np.array(unique)
raf=raf[unique==True]
decf=decf[unique==True]
probf=probf[unique==True]
r90f=r90f[unique==True]
totf=totf[unique==True]
bkgf=bkgf[unique==True]
netf=netf[unique==True]
e_netf_up=e_netf_up[unique==True]
e_netf_lo=e_netf_lo[unique==True]
expf=expf[unique==True]
crf=crf[unique==True]
e_crf_up=e_crf_up[unique==True]
e_crf_lo=e_crf_lo[unique==True]
fluxf=fluxf[unique==True]
e_fluxf_up=e_fluxf_up[unique==True]
e_fluxf_lo=e_fluxf_lo[unique==True]

ras=ras[unique==True]
decs=decs[unique==True]
probs=probs[unique==True]
r90s=r90s[unique==True]
tots=tots[unique==True]
bkgs=bkgs[unique==True]
nets=nets[unique==True]
e_nets_up=e_nets_up[unique==True]
e_nets_lo=e_nets_lo[unique==True]
exps=exps[unique==True]
crs=crs[unique==True]
e_crs_up=e_crs_up[unique==True]
e_crs_lo=e_crs_lo[unique==True]
fluxs=fluxs[unique==True]
e_fluxs_up=e_fluxs_up[unique==True]
e_fluxs_lo=e_fluxs_lo[unique==True]

rah=rah[unique==True]
dech=dech[unique==True]
probh=probh[unique==True]
r90h=r90h[unique==True]
toth=toth[unique==True]
bkgh=bkgh[unique==True]
neth=neth[unique==True]
e_neth_up=e_neth_up[unique==True]
e_neth_lo=e_neth_lo[unique==True]
exph=exph[unique==True]
crh=crh[unique==True]
e_crh_up=e_crh_up[unique==True]
e_crh_lo=e_crh_lo[unique==True]
fluxh=fluxh[unique==True]
e_fluxh_up=e_fluxh_up[unique==True]
e_fluxh_lo=e_fluxh_lo[unique==True]

#write catalog
cat=Table([id,ra_u,dec_u,poserr,probf,r90f,totf,bkgf,netf,e_netf_up,e_netf_lo,expf,crf,e_crf_up,e_crf_lo,fluxf,e_fluxf_up,e_fluxf_lo,probs,r90s,tots,bkgs,nets,e_nets_up,e_nets_lo,exps,crs,e_crs_up,e_crs_lo,fluxs,e_fluxs_up,e_fluxs_lo,probh,r90h,toth,bkgh,neth,e_neth_up,e_neth_lo,exph,crh,e_crh_up,e_crh_lo,fluxh,e_fluxh_up,e_fluxh_lo,hr,e_hr_up,e_hr_lo,out_name],names=('ID','RA','DEC','POS_ERR','PROB_F','R90_F','TOT_F','BKG_F','NET_F','E_NET_F_+','E_NET_F_-','EXP_F','CR_F','E_CR_F_+','E_CR_F_-','FLUX_F','E_FLUX_F_+','E_FLUX_F_-','PROB_S','R90_S','TOT_S','BKG_S','NET_S','E_NET_S_+','E_NET_S_-','EXP_S','CR_S','E_CR_S_+','E_CR_S_-','FLUX_S','E_FLUX_S_+','E_FLUX_S_-','PROB_H','R90_H','TOT_H','BKG_H','NET_H','E_NET_H_+','E_NET_H_-','EXP_H','CR_H','E_CR_H_+','E_CR_H_-','FLUX_H','E_FLUX_H_+','E_FLUX_H_-','HR','E_HR_+','E_HR_-','XB_ID'))
cat.write(wd+'new_mosaics_detection/cdwfs_merged_cat1_200113.fits',format='fits',overwrite=True)

print('For the positional error, R90 was used ',used_r90,' times')

hr=np.array(hr)
hr_clean=hr[hr>-98]
print(count,len(raf))
print(kenter_match, 'Kenter')
plt.figure()
plt.hist(hr_clean,bins=40)
plt.xlabel('HR')
plt.show()

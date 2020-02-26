# This script takes the merged F,S,H catalogs and applies probability cuts based on 
# simulations, computes final aperture photometry and the HR with BEHR.
import numpy as np
from astropy.io import fits
from astropy.table import Table
import sys
import subprocess as s
import os
import scipy.stats.distributions
from scipy.stats import poisson
from scipy import optimize

# Function of which we need to find the zero - solving equation 1 of Gehrels86
def f(x):
	return poisson.cdf(k,x)-(1.-CL)
	
wd='/Users/alberto/Desktop/XBOOTES/'
date = '200113'

cat=fits.open(wd+'new_mosaics_detection/cdwfs_merged_cat0_'+date+'.fits')

raf=cat[1].data['RA_F']
decf=cat[1].data['DEC_F']
probf=cat[1].data['PROB_F']

ras=cat[1].data['RA_S']
decs=cat[1].data['DEC_S']
probs=cat[1].data['PROB_S']

rah=cat[1].data['RA_H']
dech=cat[1].data['DEC_H']
probh=cat[1].data['PROB_H']

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
sys.exit()
# Consider only those sources that are to be included in the final catalog
my_sample = cat[1].data[(probf<=cutf) | (probs<=cuts) | (probh<=cuth)]
cat.close()

my_raf = my_sample['RA_F']
my_decf = my_sample['DEC_F']
my_probf = my_sample['PROB_F']
my_totf = my_sample['TOT_F']
my_bkgf = my_sample['BKG_F']
my_enetfp = my_sample['E_NET_F_+']
my_enetfn = my_sample['E_NET_F_-']
my_expf = my_sample['EXP_F']
my_crf = my_sample['CR_F']
my_ecrfp = my_sample['E_CR_F_+']
my_ecrfn = my_sample['E_CR_F_-']
my_fluxf = my_sample['FLUX_F']
my_efluxfp = my_sample['E_FLUX_F_+']
my_efluxfn = my_sample['E_FLUX_F_-']

my_ras = my_sample['RA_S']
my_decs = my_sample['DEC_S']
my_probs = my_sample['PROB_S']
my_tots = my_sample['TOT_S']
my_bkgs = my_sample['BKG_S']
my_enetsp = my_sample['E_NET_S_+']
my_enetsn = my_sample['E_NET_S_-']
my_exps = my_sample['EXP_S']
my_crs = my_sample['CR_S']
my_ecrsp = my_sample['E_CR_S_+']
my_ecrsn = my_sample['E_CR_S_-']
my_fluxs = my_sample['FLUX_S']
my_efluxsp = my_sample['E_FLUX_S_+']
my_efluxsn = my_sample['E_FLUX_S_-']

my_rah = my_sample['RA_H']
my_dech = my_sample['DEC_H']
my_probh = my_sample['PROB_H']
my_toth = my_sample['TOT_H']
my_bkgh = my_sample['BKG_H']
my_enethp = my_sample['E_NET_H_+']
my_enethn = my_sample['E_NET_H_-']
my_exph = my_sample['EXP_H']
my_crh = my_sample['CR_H']
my_ecrhp = my_sample['E_CR_H_+']
my_ecrhn = my_sample['E_CR_H_-']
my_fluxh = my_sample['FLUX_H']
my_efluxhp = my_sample['E_FLUX_H_+']
my_efluxhn = my_sample['E_FLUX_H_-']

id = np.arange(1,len(my_sample)+1)

# Choose best coordinates based on most significant detection 
# (WATCH OUT! NEED TO KEEP INTO ACCOUNT THE THRESHOLDS, OTHERWISE RISKING OF GETTING 
# COORDINATES OF HIGHER PROBABILITY, BUT BELOW THRESHOLD - and using the new coordinates, 
# losing also the original significant band detection)
probs = np.array([my_sample['PROB_F'],my_sample['PROB_S'],my_sample['PROB_H']])
mask = [cutf, cuts, cuth]
index=[]
for i in range(len(probs.T)):
	if len(probs.T[i][probs.T[i] < mask]) > 1: # Source is above threshold in more than one band, pick the most significant one
		index.append(np.argmin(probs.T[i]))
	else: # source is above threshold in just one band, pick that one regardless of the probability value
		index.append(np.where(probs.T[i] < mask)[0][0])

index=np.array(index)
#index = np.argmin(probs.T, axis=1)
'''
final_prob = my_probf.copy()
final_prob[index==1] = my_probs[index==1]
final_prob[index==2] = my_probh[index==2]

cut = np.full_like(final_prob, cutf)
cut[index==1] = cuts
cut[index==2] = cuth

print(len(final_prob[final_prob<=cut]))

sys.exit()
'''

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

# Change the NaNs and Infinites using poserr_R90
poserr[np.isnan(poserr)==True] = poserr_r90[np.isnan(poserr)==True]
poserr[np.isinf(poserr)==True] = poserr_r90[np.isinf(poserr)==True]

# With the new unique coordinates, perform aperture photometry on all sources, 
# with correct treatment of new detections,  and upper limits
print(len(ra_u[my_raf==0.0]),'Non detections in F band')
print(len(ra_u[my_ras==0.0]),'Non detections in S band')
print(len(ra_u[my_rah==0.0]),'Non detections in H band')

for band in ['broad','soft','hard']:

	# Take the coordinates
	'''
	if band == 'broad':
		ra_nd = ra_u[my_raf==0.0]
		dec_nd = dec_u[my_raf==0.0]
	elif band == 'soft':
		ra_nd = ra_u[my_ras==0.0]
		dec_nd = dec_u[my_ras==0.0]
	else:
		ra_nd = ra_u[my_rah==0.0]
		dec_nd = dec_u[my_rah==0.0]
	'''
	ra_nd = ra_u
	dec_nd = dec_u
	
	# Write out the new region file with circles and 2" radius to get R90 and ECF
	w=open(wd+'new_mosaics_detection/cdwfs_'+band+'_src_'+date+'_forR90.reg','w')
	for i in range(len(ra_nd)):
		w.write('circle('+str(ra_nd[i])+'d, '+str(dec_nd[i])+'d, 2\")\n')
	w.close()

	# Then, run dmextract to get the average R90 and average ECF
	s.call('dmextract \''+wd+'new_mosaics_detection/cdwfs_'+band+'_r90_4reb.fits[bin sky=@'+wd+'new_mosaics_detection/cdwfs_'+band+'_src_'+date+'_forR90.reg]\' outfile='+wd+'new_mosaics_detection/cdwfs_'+band+'_src_'+date+'_R90.fits clobber=yes', shell=True)
	s.call('dmextract \''+wd+'new_mosaics_detection/cdwfs_'+band+'_ecfmap_4reb.fits[bin sky=@'+wd+'new_mosaics_detection/cdwfs_'+band+'_src_'+date+'_forR90.reg]\' outfile='+wd+'new_mosaics_detection/cdwfs_'+band+'_src_'+date+'_ECF.fits clobber=yes', shell=True)

	cat = fits.open(wd+'new_mosaics_detection/cdwfs_'+band+'_src_'+date+'_R90.fits')
	r90_nd = 16.*cat[1].data['COUNTS']/cat[1].data['AREA']
	cat.close()

	cat = fits.open(wd+'new_mosaics_detection/cdwfs_'+band+'_src_'+date+'_ECF.fits')
	ecf_nd = 16.*cat[1].data['COUNTS']/cat[1].data['AREA']
	cat.close()
	
	# Write out the new region file with circles and R90 radius
	w=open(wd+'new_mosaics_detection/cdwfs_'+band+'_src_'+date+'_forAP.reg','w')
	for i in range(len(ra_nd)):
		w.write('circle('+str(ra_nd[i])+'d, '+str(dec_nd[i])+'d, '+str(r90_nd[i])+'\")\n')
	w.close()

	# Then, run dmextract to get aperture photometry (TOT, BKG, average EXP)
	s.call('dmextract \''+wd+'new_mosaics_detection/cdwfs_'+band+'_4reb.fits[bin sky=@'+wd+'new_mosaics_detection/cdwfs_'+band+'_src_'+date+'_forAP.reg]\' outfile='+wd+'new_mosaics_detection/cdwfs_'+band+'_src_'+date+'_AP.fits bkg=\''+wd+'new_mosaics_detection/cdwfs_'+band+'_bkgmap_4reb.fits[bin sky=@'+wd+'new_mosaics_detection/cdwfs_'+band+'_src_'+date+'_forAP.reg]\' clobber=yes', shell=True)
	s.call('dmextract \''+wd+'new_mosaics_detection/cdwfs_'+band+'_expomap_4reb.fits[bin sky=@'+wd+'new_mosaics_detection/cdwfs_'+band+'_src_'+date+'_forAP.reg]\' outfile='+wd+'new_mosaics_detection/cdwfs_'+band+'_src_'+date+'_EXP.fits clobber=yes', shell=True)

	cat = fits.open(wd+'new_mosaics_detection/cdwfs_'+band+'_src_'+date+'_EXP.fits')
	exp_nd = 16.*cat[1].data['COUNTS']/cat[1].data['AREA']
	cat.close()

	cat = fits.open(wd+'new_mosaics_detection/cdwfs_'+band+'_src_'+date+'_AP.fits')
	tot_nd = cat[1].data['COUNTS']
	bkg_nd = cat[1].data['BG_COUNTS']
	cat.close()

	# This is the probability of having EXACTLY cts counts given bkg, not of having
	# AT LEAST cts counts given bkg.
	prob_nd=scipy.stats.distributions.poisson.pmf(tot_nd,bkg_nd)

	if band == 'broad':
		cut = cutf
		# This is a dummy mask, here in case I need it in the future;
		# I was doing aperture photometry on the non-detections only,
		# but this introduced some misplacements in the final catalog
		# for sources with more than one band detection. So I concluded
		# that doing aperture photometry again for ALL the sources is
		# more reasonable.
		mask = my_raf>=0.0 
	elif band == 'soft':
		cut = cuts
		mask = my_ras>=0.0
	else: 
		cut = cuth
		mask = my_rah>=0.0
	
	e_cts_p = np.zeros_like(tot_nd)
	e_cts_n = np.zeros_like(tot_nd)
	
	e_bkg_p = np.zeros_like(tot_nd)
	e_bkg_n = np.zeros_like(tot_nd)
	
	net_nd = np.zeros_like(tot_nd)
	e_net_p = np.zeros_like(tot_nd)
	e_net_n = np.zeros_like(tot_nd)
	
	cr_nd = np.zeros_like(tot_nd)
	e_cr_p = np.zeros_like(tot_nd)
	e_cr_n = np.zeros_like(tot_nd)
	
	flux_nd = np.zeros_like(tot_nd)
	e_flux_p = np.zeros_like(tot_nd)
	e_flux_n = np.zeros_like(tot_nd)
	
	# If the probability is below the threshold, compute normal aperture photometry
	e_cts_p[prob_nd <= cut] = 1+np.sqrt(tot_nd[prob_nd <= cut]+0.75)
	e_cts_n[prob_nd <= cut] = np.sqrt(tot_nd[prob_nd <= cut]-0.25)

	e_bkg_p[prob_nd <= cut]=1+np.sqrt(bkg_nd[prob_nd <= cut]+0.75)
	e_bkg_n[(prob_nd <= cut) & (bkg_nd >= 0.25)] = np.sqrt(bkg_nd[(prob_nd <= cut) & (bkg_nd >= 0.25)]-0.25)
	
	net_nd[prob_nd <= cut] = tot_nd[prob_nd <= cut]-bkg_nd[prob_nd <= cut]
	e_net_p[prob_nd <= cut] = np.sqrt(e_cts_p[prob_nd <= cut]**2+e_bkg_p[prob_nd <= cut]**2) # Propagate the errors
	e_net_n[prob_nd <= cut] = np.sqrt(e_cts_n[prob_nd <= cut]**2+e_bkg_n[prob_nd <= cut]**2)

	cr_nd[prob_nd <= cut] = net_nd[prob_nd <= cut]*1.1/exp_nd[prob_nd <= cut]
	e_cr_p[prob_nd <= cut] = e_net_p[prob_nd <= cut]*1.1/exp_nd[prob_nd <= cut] # Propagate the errors
	e_cr_n[prob_nd <= cut] = e_net_n[prob_nd <= cut]*1.1/exp_nd[prob_nd <= cut]
	
	flux_nd[prob_nd <= cut] = cr_nd[prob_nd <= cut]*ecf_nd[prob_nd <= cut]
	e_flux_p[prob_nd <= cut] = e_cr_p[prob_nd <= cut]*ecf_nd[prob_nd <= cut]
	e_flux_n[prob_nd <= cut] = e_cr_n[prob_nd <= cut]*ecf_nd[prob_nd <= cut]

	# Otherwise, compute 3sigma UL
	kk = tot_nd[prob_nd > cut] # Observed total counts
	bkg_k = bkg_nd[prob_nd > cut] # Expected background counts, accurately measured
	CL = 0.9987 # Three sigma is 0.9987
	s0 = []
	for l in range(len(kk)):
		k = kk[l]
		mu = optimize.newton(f, k, maxiter=500) # using k as a guess
		s0.append(mu-bkg_k[l]) # s0 is the 3sigma upperlimit on net cts
	s0 = np.array(s0)
	
	net_nd[prob_nd > cut] = s0
	e_net_p[prob_nd > cut] = -99
	e_net_n[prob_nd > cut] = -99
	
	cr_nd[prob_nd > cut] = net_nd[prob_nd > cut]*1.1/exp_nd[prob_nd > cut]
	e_cr_p[prob_nd > cut] = -99
	e_cr_n[prob_nd > cut] = -99
	
	flux_nd[prob_nd > cut]=cr_nd[prob_nd > cut]*ecf_nd[prob_nd > cut]
	e_flux_p[prob_nd > cut] = -99
	e_flux_n[prob_nd > cut] = -99
	
	# Output stuff
	if band == 'broad':
		my_probf[mask] = prob_nd
		my_r90f[mask] = r90_nd
		my_totf[mask] = tot_nd
		my_bkgf[mask] = bkg_nd
		my_netf[mask] = net_nd
		my_enetfp[mask] = e_net_p
		my_enetfn[mask] = e_net_n
		my_expf[mask] = exp_nd
		my_crf[mask] = cr_nd
		my_ecrfp[mask] = e_cr_p
		my_ecrfn[mask] = e_cr_n
		my_fluxf[mask] = flux_nd
		my_efluxfp[mask] = e_flux_p
		my_efluxfn[mask] = e_flux_n
	elif band == 'soft':
		my_probs[mask] = prob_nd
		my_r90s[mask] = r90_nd
		my_tots[mask] = tot_nd
		my_bkgs[mask] = bkg_nd
		my_nets[mask] = net_nd
		my_enetsp[mask] = e_net_p
		my_enetsn[mask] = e_net_n
		my_exps[mask] = exp_nd
		my_crs[mask] = cr_nd
		my_ecrsp[mask] = e_cr_p
		my_ecrsn[mask] = e_cr_n
		my_fluxs[mask] = flux_nd
		my_efluxsp[mask] = e_flux_p
		my_efluxsn[mask] = e_flux_n
	else:
		my_probh[mask] = prob_nd
		my_r90h[mask] = r90_nd
		my_toth[mask] = tot_nd
		my_bkgh[mask] = bkg_nd
		my_neth[mask] = net_nd
		my_enethp[mask] = e_net_p
		my_enethn[mask] = e_net_n
		my_exph[mask] = exp_nd
		my_crh[mask] = cr_nd
		my_ecrhp[mask] = e_cr_p
		my_ecrhn[mask] = e_cr_n
		my_fluxh[mask] = flux_nd
		my_efluxhp[mask] = e_flux_p
		my_efluxhn[mask] = e_flux_n


# Compute Hardness Ratios
hr,e_hr_lo,e_hr_up = [],[],[]
r90ratio=(my_r90h/my_r90s)**2 # This takes into account that R90 is larger in hard band into the HR computation
os.chdir('/Users/alberto/BEHR/') 
for j in range(len(r90ratio)):

	# Compute HR with BEHR
	behr='./BEHR softsrc='+str(int(my_tots[j]))+' hardsrc='+str(int(my_toth[j]))+' softbkg='+str(int(round(my_bkgs[j])))+' hardbkg='+str(int(round(my_bkgh[j])))+' softarea=1 hardarea='+str(r90ratio[j])+' output='+str(j)+' outputHR=true'
	s.call(behr,shell=True)

	a,b,c=np.genfromtxt('/Users/alberto/BEHR/'+str(j)+'_HR.txt',unpack=True,usecols=[2,3,4]) # Use median
	hr.append(a)
	e_hr_lo.append(a-b)
	e_hr_up.append(c-a)

#write catalog
cat=Table([id,ra_u,dec_u,poserr,my_probf,my_r90f,my_totf,my_bkgf,my_netf,my_enetfp,my_enetfn,my_expf,my_crf,my_ecrfp,my_ecrfn,my_fluxf,my_efluxfp,my_efluxfn,my_probs,my_r90s,my_tots,my_bkgs,my_nets,my_enetsp,my_enetsn,my_exps,my_crs,my_ecrsp,my_ecrsn,my_fluxs,my_efluxsp,my_efluxsn,my_probh,my_r90h,my_toth,my_bkgh,my_neth,my_enethp,my_enethn,my_exph,my_crh,my_ecrhp,my_ecrhn,my_fluxh,my_efluxhp,my_efluxhn,hr,e_hr_up,e_hr_lo],names=('ID','RA','DEC','POS_ERR','PROB_F','R90_F','TOT_F','BKG_F','NET_F','E_NET_F_+','E_NET_F_-','EXP_F','CR_F','E_CR_F_+','E_CR_F_-','FLUX_F','E_FLUX_F_+','E_FLUX_F_-','PROB_S','R90_S','TOT_S','BKG_S','NET_S','E_NET_S_+','E_NET_S_-','EXP_S','CR_S','E_CR_S_+','E_CR_S_-','FLUX_S','E_FLUX_S_+','E_FLUX_S_-','PROB_H','R90_H','TOT_H','BKG_H','NET_H','E_NET_H_+','E_NET_H_-','EXP_H','CR_H','E_CR_H_+','E_CR_H_-','FLUX_H','E_FLUX_H_+','E_FLUX_H_-','HR','E_HR_+','E_HR_-'))
cat.write(wd+'new_mosaics_detection/cdwfs_merged_cat1_'+date+'_provvisorio2.fits',format='fits',overwrite=True)
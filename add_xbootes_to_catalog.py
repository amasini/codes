import numpy as np
import sys
import subprocess as s
import scipy.stats.distributions
import time
import os
from astropy.io import fits
from astropy.table import Table
from astropy.table import Column
from scipy.stats import poisson
from scipy import optimize

# Function of which we need to find the zero - solving equation 1 of Gehrels86
def f(x):
	return poisson.cdf(k,x)-(1.-CL)
	
def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*xx*3600)**2 +((pointa[1]-pointb[1])*3600)**2)
    
wd = '/Users/alberto/Desktop/XBOOTES/'
date = '200113'

#### Part 1: MATCH XBOOTES WITH CDWFS
'''
tin = time.time()

kcat=fits.open(wd+'xbootes_kenter+05.fits')
ra_k=kcat[1].data['RAJ2000']
dec_k=kcat[1].data['DEJ2000']
name_k=kcat[1].data['CXOXB']
kcat.close()

cdwfs_cat = fits.open(wd+'new_mosaics_detection/cdwfs_merged_cat1_'+date+'_provvisorio2.fits')
id_c = cdwfs_cat[1].data['ID']
ra_c = cdwfs_cat[1].data['RA']
dec_c = cdwfs_cat[1].data['DEC']
r90_c = cdwfs_cat[1].data['R90_F']
cdwfs_cat.close()

w=open(wd+'new_mosaics_detection/cdwfs_merged_cat1_'+date+'_provvisorio2.reg','w')
for i in range(len(ra_c)):
	w.write('circle('+str(ra_c[i])+'d,'+str(dec_c[i])+'d,'+str(r90_c[i])+'\") \n')
w.close()

# Look for CDWFS sources into XB
out_name, matched, c_x_d = [],[],[]
for j in range(len(ra_c)): # Look up into CDWFS catalog
	
	cdwfs = [ra_c[j],dec_c[j]]
	found=False
	
	ra_cont,dec_cont,d,name=[],[],[],[]
	for k in range(len(ra_k)):
	
		kenter=[ra_k[k],dec_k[k]]
		
		if distance(cdwfs,kenter) < 1.1*r90_c[j]:
		
			if found == False: # This is the first match
			
				found = True
				ra_cont.append(ra_k[k])
				dec_cont.append(dec_k[k])
				d.append(distance(cdwfs,kenter))
				name.append(name_k[k])
				
			elif found == True:
			
				ra_cont.append(ra_k[k])
				dec_cont.append(dec_k[k])
				d.append(distance(cdwfs,kenter))
				name.append(name_k[k])
				print('Finding a double counterpart for '+str(id_c[j])+' ('+str(round(ra_c[j],5))+', '+str(round(dec_c[j],5))+'), IDs',name)
	
	if found == True:
	
		matched.append(1)
		if len(name) > 1: # If more than one CDWFS counterpart, keep the closest
		
			for l in range(len(name)): 
				if min(d) == d[l]:
					out_name.append(name[l])
					c_x_d.append(d[l])
		else:
			out_name.append(name[0])
			c_x_d.append(d[0])
	
	else:
		matched.append(0)
		out_name.append('-99')
		c_x_d.append(-99)
		
matched = np.array(matched)
out_name = np.array(out_name)
c_x_d = np.array(c_x_d)
c_unique, unq_idx, unq_cnt = np.unique(out_name[out_name != '-99'], return_inverse=True, return_counts=True)

dup_ids = c_unique[unq_cnt > 1]
unique_ids = c_unique[unq_cnt == 1]

# Return a boolean array of which elements are in unique array (True) and which are not (False)
# We want to loop on the Falses in order to see if there's something to be added to the 
# catalog
c = np.in1d(name_k,unique_ids)

ra_tbc = ra_k[c==False]
dec_tbc = dec_k[c==False]
name_tbc = name_k[c==False]
print(len(ra_tbc),'XBOOTES sources are missed, and are to be checked.')
w=open(wd+'xbootes_missing.dat','w')
w1=open(wd+'xbootes_missing.reg','w')
w.write('CXOXB \t RA \t DEC\n')
for i in range(len(ra_tbc)):
	w.write(str(name_tbc[i])+' \t '+str(ra_tbc[i])+' \t '+str(dec_tbc[i])+'\n')
	w1.write('circle('+str(ra_tbc[i])+'d, '+str(dec_tbc[i])+'d, 2\")\n')
w.close()
w1.close()

print(len(ra_c[matched == 1]),'CDWFS Matched')
print(len(ra_c[matched == 0]),'CDWFS Unmatched')
print(len(out_name[out_name != '-99']),'XB counterparts')
print(len(unique_ids),'unique XB counterparts')
print((time.time()-tin)/60.,'minutes for the match.')

# Add XB column to catalog and write out FITS file
table = Table.read(wd+'new_mosaics_detection/cdwfs_merged_cat1_'+date+'_provvisorio2.fits', format='fits')
#table.remove_column('XB_ID')
#table.remove_column('XB_CDWFS_D')

t2 = Column(name='XB_ID', data=out_name)
t3 = Column(name='XB_CDWFS_D', data=c_x_d)
table.add_column(t2)
table.add_column(t3)
table.write(wd+'new_mosaics_detection/cdwfs_merged_cat1_'+date+'_provvisorio2.fits', format='fits', overwrite='True')

sys.exit()
'''

#### Part 2: filter XBOOTES to create list of sources to be added to the catalog
'''
tin=time.time()

# This file contains the N XBOOTES sources missing from CDWFS
rak,deck=np.genfromtxt(wd+'xbootes_missing.dat',unpack=True, skip_header=1,usecols=[1,2])
namek=np.genfromtxt(wd+'xbootes_missing.dat',unpack=True, skip_header=1,usecols=0, dtype='str')

# Probability thresholds
cutf, cuts, cuth = 10**(-4.63),10**(-4.57),10**(-4.40)

r90,prob,tot,bkg=[],[],[],[]
for band in ['broad','soft','hard']:
	
	# Take the coordinates
	ra_nd = rak
	dec_nd = deck

	# Then, run dmextract to get the average R90
	s.call('dmextract \''+wd+'new_mosaics_detection/cdwfs_'+band+'_r90_4reb.fits[bin sky=@'+wd+'xbootes_missing.reg]\' outfile='+wd+'new_mosaics_detection/cdwfs_'+band+'_src_'+date+'_R90.fits clobber=yes', shell=True)

	cat = fits.open(wd+'new_mosaics_detection/cdwfs_'+band+'_src_'+date+'_R90.fits')
	r90_nd = 16.*cat[1].data['COUNTS']/cat[1].data['AREA']
	cat.close()

	# Write out the new region file with circles and R90 radius
	w=open(wd+'xbootes_missing_AP.reg','w')
	for i in range(len(ra_nd)):
		w.write('circle('+str(ra_nd[i])+'d, '+str(dec_nd[i])+'d, '+str(r90_nd[i])+'\")\n')
	w.close()

	# Then, run dmextract to get aperture photometry (TOT, BKG)
	s.call('dmextract \''+wd+'new_mosaics_detection/cdwfs_'+band+'_4reb.fits[bin sky=@'+wd+'xbootes_missing_AP.reg]\' outfile='+wd+'new_mosaics_detection/cdwfs_'+band+'_src_'+date+'_AP.fits bkg=\''+wd+'new_mosaics_detection/cdwfs_'+band+'_bkgmap_4reb.fits[bin sky=@'+wd+'xbootes_missing_AP.reg]\' clobber=yes', shell=True)
	
	cat = fits.open(wd+'new_mosaics_detection/cdwfs_'+band+'_src_'+date+'_AP.fits')
	tot_nd = cat[1].data['COUNTS']
	bkg_nd = cat[1].data['BG_COUNTS']
	cat.close()

	# This is the probability of having EXACTLY cts counts given bkg, not of having
	# AT LEAST cts counts given bkg.
	prob_nd=scipy.stats.distributions.poisson.pmf(tot_nd,bkg_nd)

	r90.append(r90_nd)
	prob.append(prob_nd)
	tot.append(tot_nd)
	bkg.append(bkg_nd)
	
w=open(wd+'xbootes_tbadded.dat','w')
w.write('CXOXB \t RA \t DEC \t PROB_F \t R90_F \t TOT_F \t BKG_F \t PROB_S \t R90_S \t TOT_S \t BKG_S \t PROB_H \t R90_H \t TOT_H \t BKG_H \n')
keep = 0
for i in range(len(namek)):

	# If the source is above threshold in at least one band,
	if (prob[0][i]<=cutf) | (prob[1][i]<=cuts) | (prob[2][i]<=cuth):
		keep=keep+1		
	
		w.write(str(namek[i])+' \t '+str(rak[i])+' \t '+str(deck[i])+' \t '+str(prob[0][i])+' \t '+str(r90[0][i])+' \t '+str(tot[0][i])+' \t '+str(bkg[0][i])+' \t '+str(prob[1][i])+' \t '+str(r90[1][i])+' \t '+str(tot[1][i])+' \t '+str(bkg[1][i])+' \t '+str(prob[2][i])+' \t '+str(r90[2][i])+' \t '+str(tot[2][i])+' \t '+str(bkg[2][i])+'\n')

w.close()
print(keep,'XBOOTES sources which satisfy reliability cut.')
print((time.time()-tin),'seconds for this part.')
sys.exit()
'''

#### Part 3: compute all the quantities (with errors and upperlimits)
'''
tin = time.time()

namek,rak,deck,probf,r90f,totf,bkgf,probs,r90s,tots,bkgs,probh,r90h,toth,bkgh=np.genfromtxt(wd+'xbootes_tbadded.dat', skip_header=1, unpack=True)
namek=np.genfromtxt(wd+'xbootes_tbadded.dat',unpack=True, skip_header=1,usecols=0, dtype='str')
dd = np.zeros_like(rak) # Array of distance between XB and CDWFS source - it's zero by definitiion here

# Probability thresholds
cutf, cuts, cuth = 10**(-4.63),10**(-4.57),10**(-4.40)

cat0 = fits.open(wd+'new_mosaics_detection/cdwfs_merged_cat1_'+date+'_provvisorio2-cp.fits')
howmany = len(cat0[1].data['RA'])+1
cat0.close()

# Consider only those sources that are to be included in the final catalog

id = np.arange(howmany,howmany+len(namek))

# Choose best coordinates based on most significant detection - NOT NEEDED!
probb = np.array([probf,probs,probh])
index = np.argmin(probb.T, axis=1)

ra_u = rak
dec_u = deck

# Positional errors - use appropriate R90 (F for index=0, S for index=1 and H for index=2)

# Compute R50 from R90
r90 = r90f.copy()
r90[index == 1] = r90s[index==1]
r90[index == 2] = r90h[index==2]

netf = totf - bkgf
nets = tots - bkgs
neth = toth - bkgh

net_u = netf.copy()
net_u[index == 1] = nets[index==1]
net_u[index == 2] = neth[index==2]

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

# Then, run dmextract to get the NET counts within R50
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
for band in ['broad','soft','hard']:

	# Take the coordinates

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
	elif band == 'soft':
		cut = cuts
	else: 
		cut = cuth
	
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
		my_probf = prob_nd
		my_r90f = r90_nd
		my_totf = tot_nd
		my_bkgf = bkg_nd
		my_netf = net_nd
		my_enetfp = e_net_p
		my_enetfn = e_net_n
		my_expf = exp_nd
		my_crf = cr_nd
		my_ecrfp = e_cr_p
		my_ecrfn = e_cr_n
		my_fluxf = flux_nd
		my_efluxfp = e_flux_p
		my_efluxfn = e_flux_n
	elif band == 'soft':
		my_probs = prob_nd
		my_r90s = r90_nd
		my_tots = tot_nd
		my_bkgs = bkg_nd
		my_nets = net_nd
		my_enetsp = e_net_p
		my_enetsn = e_net_n
		my_exps = exp_nd
		my_crs = cr_nd
		my_ecrsp = e_cr_p
		my_ecrsn = e_cr_n
		my_fluxs = flux_nd
		my_efluxsp = e_flux_p
		my_efluxsn = e_flux_n
	else:
		my_probh = prob_nd
		my_r90h = r90_nd
		my_toth = tot_nd
		my_bkgh = bkg_nd
		my_neth = net_nd
		my_enethp = e_net_p
		my_enethn = e_net_n
		my_exph = exp_nd
		my_crh = cr_nd
		my_ecrhp = e_cr_p
		my_ecrhn = e_cr_n
		my_fluxh = flux_nd
		my_efluxhp = e_flux_p
		my_efluxhn = e_flux_n

# Compute Hardness Ratios
hr,e_hr_lo,e_hr_up = [],[],[]
r90ratio=(r90h/r90s)**2 # This takes into account that R90 is larger in hard band into the HR computation
os.chdir('/Users/alberto/BEHR/') 
for j in range(len(r90ratio)):

	# Compute HR with BEHR
	behr='./BEHR softsrc='+str(int(my_tots[j]))+' hardsrc='+str(int(my_toth[j]))+' softbkg='+str(int(round(my_bkgs[j])))+' hardbkg='+str(int(round(my_bkgh[j])))+' softarea=1 hardarea='+str(r90ratio[j])+' output='+str(j+howmany)+' outputHR=true'
	s.call(behr,shell=True)

	a,b,c=np.genfromtxt('/Users/alberto/BEHR/'+str(j+howmany)+'_HR.txt',unpack=True,usecols=[2,3,4]) # Use median
	hr.append(a)
	e_hr_lo.append(a-b)
	e_hr_up.append(c-a)

#write catalog
cat=Table([id,ra_u,dec_u,poserr,my_probf,my_r90f,my_totf,my_bkgf,my_netf,my_enetfp,my_enetfn,my_expf,my_crf,my_ecrfp,my_ecrfn,my_fluxf,my_efluxfp,my_efluxfn,my_probs,my_r90s,my_tots,my_bkgs,my_nets,my_enetsp,my_enetsn,my_exps,my_crs,my_ecrsp,my_ecrsn,my_fluxs,my_efluxsp,my_efluxsn,my_probh,my_r90h,my_toth,my_bkgh,my_neth,my_enethp,my_enethn,my_exph,my_crh,my_ecrhp,my_ecrhn,my_fluxh,my_efluxhp,my_efluxhn,hr,e_hr_up,e_hr_lo,namek,dd],names=('ID','RA','DEC','POS_ERR','PROB_F','R90_F','TOT_F','BKG_F','NET_F','E_NET_F_+','E_NET_F_-','EXP_F','CR_F','E_CR_F_+','E_CR_F_-','FLUX_F','E_FLUX_F_+','E_FLUX_F_-','PROB_S','R90_S','TOT_S','BKG_S','NET_S','E_NET_S_+','E_NET_S_-','EXP_S','CR_S','E_CR_S_+','E_CR_S_-','FLUX_S','E_FLUX_S_+','E_FLUX_S_-','PROB_H','R90_H','TOT_H','BKG_H','NET_H','E_NET_H_+','E_NET_H_-','EXP_H','CR_H','E_CR_H_+','E_CR_H_-','FLUX_H','E_FLUX_H_+','E_FLUX_H_-','HR','E_HR_+','E_HR_-','XB_ID','XB_CDWFS_D'))
cat.write(wd+'xbootes_tbadded.fits',format='fits',overwrite=True)

print((time.time()-tin),'seconds for this part.')
'''
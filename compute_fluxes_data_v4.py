#  APERTURE PHOTOMETRY INDEPENDENTLY FOR THE THREE BANDS

import numpy as np
from astropy.io import fits
from astropy.table import Table
import scipy
import scipy.stats.distributions
import subprocess as s
import sys

wd="/Users/alberto/Desktop/XBOOTES/"

band='hard'
date='200113'

wd = '/Users/alberto/Desktop/XBOOTES/new_mosaics_detection/'

# Open the output from wavdetect
file = fits.open(wd+'cdwfs_'+band+'_src_'+date+'.fits')
data = file[1].data
file.close()

# Take the coordinates
ra = data['RA']
dec = data['DEC']

# Write out the new region file with circles and 2" radius to get R90 and ECF
w=open(wd+'cdwfs_'+band+'_src_'+date+'_forR90.reg','w')
for i in range(len(ra)):
	w.write('circle('+str(ra[i])+'d, '+str(dec[i])+'d, 2\")\n')
w.close()

# Then, run dmextract to get the average R90 and average ECF
s.call('dmextract \''+wd+'cdwfs_'+band+'_r90_4reb.fits[bin sky=@'+wd+'cdwfs_'+band+'_src_'+date+'_forR90.reg]\' outfile='+wd+'cdwfs_'+band+'_src_'+date+'_R90.fits clobber=yes', shell=True)
s.call('dmextract \''+wd+'cdwfs_'+band+'_ecfmap_4reb.fits[bin sky=@'+wd+'cdwfs_'+band+'_src_'+date+'_forR90.reg]\' outfile='+wd+'cdwfs_'+band+'_src_'+date+'_ECF.fits clobber=yes', shell=True)

cat = fits.open(wd+'cdwfs_'+band+'_src_'+date+'_R90.fits')
r90 = 16.*cat[1].data['COUNTS']/cat[1].data['AREA']
cat.close()

cat = fits.open(wd+'cdwfs_'+band+'_src_'+date+'_ECF.fits')
ecf = 16.*cat[1].data['COUNTS']/cat[1].data['AREA']
cat.close()

# Write out the new region file with circles and R90 radius
w=open(wd+'cdwfs_'+band+'_src_'+date+'_forAP.reg','w')
for i in range(len(ra)):
	w.write('circle('+str(ra[i])+'d, '+str(dec[i])+'d, '+str(r90[i])+'\")\n')
w.close()

# Then, run dmextract to get aperture photometry (TOT, BKG, average EXP)
s.call('dmextract \''+wd+'cdwfs_'+band+'_4reb.fits[bin sky=@'+wd+'cdwfs_'+band+'_src_'+date+'_forAP.reg]\' outfile='+wd+'cdwfs_'+band+'_src_'+date+'_AP.fits bkg=\''+wd+'cdwfs_'+band+'_bkgmap_4reb.fits[bin sky=@'+wd+'cdwfs_'+band+'_src_'+date+'_forAP.reg]\' clobber=yes', shell=True)
s.call('dmextract \''+wd+'cdwfs_'+band+'_expomap_4reb.fits[bin sky=@'+wd+'cdwfs_'+band+'_src_'+date+'_forAP.reg]\' outfile='+wd+'cdwfs_'+band+'_src_'+date+'_EXP.fits clobber=yes', shell=True)

cat = fits.open(wd+'cdwfs_'+band+'_src_'+date+'_EXP.fits')
exp = 16.*cat[1].data['COUNTS']/cat[1].data['AREA']
cat.close()

cat = fits.open(wd+'cdwfs_'+band+'_src_'+date+'_AP.fits')
tot = cat[1].data['COUNTS']
bkg = cat[1].data['BG_COUNTS']
cat.close()

e_cts_p=1+np.sqrt(tot+0.75) # Gehrels+86 1sigma errors
e_cts_n=np.sqrt(tot-0.25)
e_bkg_p=1+np.sqrt(bkg+0.75)
e_bkg_n=np.zeros_like(bkg)
e_bkg_n[bkg >= 0.25] = np.sqrt(bkg[bkg >= 0.25]-0.25)	

# This is the probability of having EXACTLY cts counts given bkg, not of having
# AT LEAST cts counts given bkg.
prob=scipy.stats.distributions.poisson.pmf(tot,bkg)

net=tot-bkg
e_net_p=np.sqrt(e_cts_p**2+e_bkg_p**2) # Propagate the errors
e_net_n=np.sqrt(e_cts_n**2+e_bkg_n**2)

cr=net*1.1/exp
e_cr_p=e_net_p*1.1/exp # Propagate the errors
e_cr_n=e_net_n*1.1/exp
flux=cr*ecf
e_flux_p=e_cr_p*ecf
e_flux_n=e_cr_n*ecf

#write catalog
cat=Table([ra,dec,prob,r90,tot,bkg,net,e_net_p,e_net_n,exp,cr,e_cr_p,e_cr_n,flux,e_flux_p,e_flux_n],names=('RA','DEC','PROB','AV_R90','TOT','BKG','NET','E_NET_+','E_NET_-','EXP','CR','E_CR_+','E_CR_-','FLUX','E_FLUX_+','E_FLUX_-'))
cat.write(wd+'cdwfs_'+band+'_cat0_'+date+'.fits',format='fits',overwrite=True)
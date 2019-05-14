# Positional errors for CDWFS, try again with zero-order approach - now included in create_catalog.py
##### OBSOLETE SCRIPT
import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from ciao_contrib.region.check_fov import FOVFiles
import subprocess as s 
import os
import time
#from scipy.optimize import curve_fit

wd = '/Users/alberto/Desktop/XBOOTES/'

# take the cdwfs catalog
cat=fits.open(wd+'cdwfs_merged_cat1.fits')
ra=cat[1].data['RA']
dec=cat[1].data['DEC']
probf=cat[1].data['PROB_F']
probs=cat[1].data['PROB_S']
probh=cat[1].data['PROB_H']
r90f=cat[1].data['R90_F']
r90s=cat[1].data['R90_S']
r90h=cat[1].data['R90_H']

r90=r90f
#print(len(r90[r90==0.0]))
r90[probs<probf]=r90s[probs<probf]
#print(len(r90[r90==0.0]))
r90[(probh<probs) & (probh<probf)]=r90h[(probh<probs) & (probh<probf)]
#print(len(r90[r90==0.0]))
#sys.exit()
cat.close()

r50=r90*5./9. # zero-order approx for r50
imp,used_r90=0,0
poserr=[]
tin=time.time()
for i in range(len(ra)):
	# try to extract total and bkg counts from mosaics
	imagemap=wd+'mosaic_broad_4rebinned.fits'
	backmap=wd+'cdwfs_broad_bkgmap_4rebinned.fits'

	s.call('dmextract infile="'+imagemap+'[bin pos=circle('+str(ra[i])+'d,'+str(dec[i])+'d,'+str(r50[i]*0.000277778)+'d)]" mode=h bkg="'+backmap+'[bin pos=circle('+str(ra[i])+'d,'+str(dec[i])+'d,'+str(r50[i]*0.000277778)+'d)]" outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
	cts=s.check_output('dmlist "counts.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
	bkg=s.check_output('dmlist "counts.fits[cols BG_COUNTS]" data,clean | grep -v BG_COUNTS',shell=True)
	cts,bkg=float(cts),float(bkg)
	net=cts-bkg
	if net > 0.:
		poserr.append(r50[i]/np.sqrt(net))
	else:
		used_r90=used_r90+1
		print('#'*20)
		print('r50 not suitable here, try with r90')
		print('#'*20)
		
		s.call('dmextract infile="'+imagemap+'[bin pos=circle('+str(ra[i])+'d,'+str(dec[i])+'d,'+str(r90[i]*0.000277778)+'d)]" mode=h bkg="'+backmap+'[bin pos=circle('+str(ra[i])+'d,'+str(dec[i])+'d,'+str(r90[i]*0.000277778)+'d)]" outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
		cts=s.check_output('dmlist "counts.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
		bkg=s.check_output('dmlist "counts.fits[cols BG_COUNTS]" data,clean | grep -v BG_COUNTS',shell=True)
		cts,bkg=float(cts),float(bkg)
		net=cts-bkg
		if net > 0.:
			poserr.append(r90[i]/np.sqrt(net))
		else:
			imp=imp+1
			poserr.append(0.1)
			print('#'*20)
			print('Not even r90 is suitable here',ra[i],dec[i],r90[i],cts,bkg,net)
			print('#'*20)

print((time.time()-tin)/60., 'minutes')
print(used_r90,'times used R90')
print(imp,'times problematic situations')
bins=np.logspace(np.log10(0.1),np.log10(10),50)

plt.figure()
plt.hist(poserr,bins=bins)
plt.xscale('log')
plt.show()
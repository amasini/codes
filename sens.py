# Script to compute sensitivity from bkg map for CDWFS and XBOOTES
import numpy as np
import subprocess as s
import matplotlib.pyplot as plt
from astropy import wcs
from astropy.wcs import WCS
from astropy.io import fits
import sys
from scipy.special import gammainc
from matplotlib.colors import LogNorm
import time
from scipy.stats import poisson
from scipy import interpolate
from collections import OrderedDict

wd='/Users/alberto/Desktop/XBOOTES/'

########### PARAMETERS ############
what='cdwfs'
band='hard'
rpsf='r90'
logcut=-4.2
pthresh=10**(logcut)

rebin_factor=4.
scale=(0.492/3600.)*rebin_factor #pixel size in deg
arcsec2pix=scale*3600.
###################################

print('Doing',what,'in the',band,'band, using '+rpsf+' and '+str(round(logcut,1)))


### Take expo map (exposure has average exposure in 4x4 pixels in s)
expmap=fits.open(wd+'new_mosaics_detection/'+what+'_'+band+'_expomap_4reb.fits')
exp=expmap[0].data
exp[np.isnan(exp)]=0.0 #put nans to 0
#exp[exp > 9e4] = 0.0
expmap.close()

### Take average psfmap squared (in arcsec)
psfmap=fits.open(wd+'new_mosaics_detection/'+what+'_'+band+'_r90sq_4reb.fits')
psf=psfmap[0].data
psfmap.close()

# No need to convert to arcsec if psfmap in arsec is used
if rpsf == 'r50':
	sourcearea=np.pi*psf*(5./9.)**2 # Use a zeroth-order approx for r50
	fpsf=0.5
elif rpsf == 'r70':
	sourcearea=np.pi*psf*(7./9.)**2 # Use a zeroth-order approx for r70
	fpsf=0.7
else:
	sourcearea=np.pi*psf # Use r90
	fpsf=0.9

### Take bkg map (bkgmap has SUMMED bkg in 4x4 pixels, like data)
bkgmap=fits.open(wd+'new_mosaics_detection/'+what+'_'+band+'_bkgmap_4reb.fits')
bkg=bkgmap[0].data
backheader=bkgmap[0].header
bkg[np.isnan(bkg)]=0.0 #put nans to 0
bkgmap.close()

### Take energy conversion factors map (weighted average of countrate to flux factor)
ecfmap=fits.open(wd+'new_mosaics_detection/'+what+'_'+band+'_ecfmap_4reb.fits')
ecf=ecfmap[0].data
ecfmap.close()
'''
ecf=np.zeros_like(exp,dtype=float)
ecf[exp==0.0]=0.0
if band == 'broad':
	#ecf[(exp<=12000.) & (exp>0.)]= 9.325E-12 #Cy3 Gamma=1.4
	#ecf[(exp<=12000.) & (exp>0.)]= 1.016E-11 #Cy8 Gamma=1.4
	ecf[exp>0.]=1.411E-11 #Cy18 Gamma=1.4
	pthresh=8e-5
elif band == 'soft':
	#ecf[(exp<=12000.) & (exp>0.)]=4.481E-12 #Cy3 Gamma=1.4
	ecf[exp>0.]=8.707E-12 #Cy18 Gamma=1.4
	pthresh=6e-4
else:
	#ecf[(exp<=12000.) & (exp>0.)]=1.973E-11 #Cy3 Gamma=1.4
	ecf[exp>0.]=2.022E-11 #Cy18 Gamma=1.4
	pthresh=4e-5
'''

### Bkg is in cts/pixel; bkg2 is in cts [(cts/arcsec2)*arcsec2]
pixarea=arcsec2pix*arcsec2pix
bkg2=bkg*(sourcearea/pixarea)

pcts=np.zeros_like(bkg2,dtype=float)

# GAMMA INCOMPLETE METHOD
# Given background and threshold probability, find how many counts are needed to have P < threhsold
tin=time.time()
tot=np.arange(1,101)
f=np.logspace(np.log10(1e-17),np.log10(2e-10),101)
totprob=np.zeros_like(f)
for i in range(bkg2.shape[0]):
	for j in range(bkg2.shape[1]):
		if exp[i][j] != 0.0:
			trial=gammainc(tot,bkg2[i][j])

			pcts[i][j]=np.min(tot[trial < pthresh])

			T=bkg2[i][j]+f*exp[i][j]/ecf[i][j]*fpsf

			prob=gammainc(pcts[i][j],T)
			totprob=totprob+prob
		else:
			pcts[i][j]=0
				
print((time.time()-tin)/60.,'minutes for the georkakais map.')

totprob=totprob*2.988e-7 # convert to area (1 4x4 pix is 2.988e-7 deg2)

### Write out result
w=open(wd+what+'_'+band+'_sens_'+str(round(logcut,1))+'_geo.dat','w')
for j in range(len(f)):
	w.write(str(f[j])+' \t '+str(totprob[j])+'\n')
w.close()

#totproblog=np.log10(totprob)
#plt.figure()
#plt.plot(f,totproblog,'k-')
#plt.xscale('log')
#plt.show()

# FOLLOWING CIVANO, EASIER METHOD
# Given background and threshold probability, find how many counts are needed to have P < threhsold
tin=time.time()
tot=np.arange(1,101)
f=np.logspace(np.log10(1e-17),np.log10(2e-10),101)
fcen = list((f[i+1]+f[i])/2. for i in range(len(f)-1))
fcen = np.array(fcen)
#totprob=np.zeros_like(f)
for i in range(bkg2.shape[0]):
	for j in range(bkg2.shape[1]):
		if exp[i][j] != 0.0:
			trial=gammainc(tot,bkg2[i][j])
			
			pct = np.min(tot[trial < pthresh])
			
			pcts[i][j] = ((pct - bkg2[i][j])*ecf[i][j])/(exp[i][j]*fpsf)
			
		else:
			pcts[i][j]=0
			
print((time.time()-tin)/60.,'minutes for the normal map.')

a,b = np.histogram(pcts, bins = f)
totprob = np.cumsum(a)

totprob=totprob*2.988e-7 # convert to area (1 4x4 pix is 2.988e-7 deg2)

### Write out result
w=open(wd+what+'_'+band+'_sens_'+str(round(logcut,1))+'_civ.dat','w')
for j in range(len(fcen)):
	w.write(str(fcen[j])+' \t '+str(totprob[j])+'\n')
w.close()

totproblog=np.log10(totprob)
plt.figure()
plt.plot(fcen,totproblog,'k-')
plt.xscale('log')
plt.show()
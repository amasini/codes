# Script to compute sensitivity from bkg map for XBOOTES
import numpy as np
import subprocess as s
import matplotlib.pyplot as plt
from astropy import wcs
from astropy.wcs import WCS
from astropy.io import fits
import sys
from scipy.special import gammainc
from matplotlib.colors import LogNorm
from ciao_contrib.region.check_fov import FOVFiles
import time
from scipy.stats import poisson

wd='/Users/alberto/Desktop/XBOOTES/'

#################################
what='xbootes'
band='broad'
gamma='1.4'
cy='Cy3' # XBOOTES was done in Cy3, CDWFS in Cy18
r50=True
'''
ok=False
if band == 'broad':
	if cy == 'Cy3':
		if gamma == '1.4':
			ecf=9.325E-12 # gamma=1.4 for Chandra Cy3, broad band
			ok=True
		elif gamma == '1.8':
			ecf= 7.796E-12 # gamma=1.8 for Chandra Cy3, broad band
			ok=True
	elif cy == 'Cy18':
		if gamma == '1.4':
			ecf=1.411E-11 # gamma=1.4 for Chandra Cy18, broad band
			ok=True
		elif gamma == '1.8':
			ecf=1.307E-11 # gamma=1.8 for Chandra Cy18, broad band
			ok=True
elif band == 'soft':
	if cy == 'Cy18':
		if gamma == '1.4':
			ecf=8.707E-12 # gamma=1.4 for Chandra Cy18, soft band
			ok=True
if ok == False:
	print('Something wrong with input parameters. Check.')
	sys.exit()
'''
rebin_factor=4.
scale=(0.492/3600.)*rebin_factor #pixel size in deg
arcsec2pix=scale*3600
print(arcsec2pix)
#################################

#take expo map (exposure has average exposure in 4x4 pixels in s)
#expmap=fits.open(wd+'new_mosaics_detection/'+what+'_'+band+'_expomap_4reb.fits')
#exp=expmap[0].data
expmap=fits.open(wd+'data/3596/repro_new_asol/out/acisf03596_broad_expomap_4reb.fits')
exp=expmap[0].data/16.0
exp[np.isnan(exp)]=0.0 #put nans to 0

ecf=np.zeros_like(exp,dtype=float)
ecf[exp==0.0]=0.0
ecf[(exp<=12000.) & (exp>0.)]=9.81E-12
ecf[exp>12000.]=1.329E-11

#plt.figure()
#plt.imshow(ecf,origin='lower')
#plt.show()
#sys.exit()

#take psfmap (in arcsec)
#psfmap=fits.open(wd+'psfmaps/'+what+'_'+band+'_r90sq_4reb.fits')
psfmap=fits.open(wd+'psfmaps/03596_r90sq_4reb.fits')
psf=psfmap[0].data

#take bkg map (bkgmap has SUMMED bkg in 4x4 pixels, like data)
#bkgmap=fits.open(wd+'new_mosaics_detection/'+what+'_'+band+'_bkgmap_4reb.fits')
bkgmap=fits.open(wd+'data/3596/repro_new_asol/out/acisf03596_broad_bkgmap_4reb.fits')
bkg=bkgmap[0].data
backheader=bkgmap[0].header
bkg[np.isnan(bkg)]=0.0 #put nans to 0

#no need to convert to arcsec if psfmap in arsec is used
if r50 == True:
	sourcearea=np.pi*psf*(5./9.)**2 # this uses a zeroth-order approx for r50
elif r50 == False:
	sourcearea=np.pi*psf # this uses r90

### bkg is in cts/pixel; bkg2 is in cts [(cts/arcsec2)*arcsec2]
pixarea=arcsec2pix*arcsec2pix
bkg2=bkg*(sourcearea/pixarea)

#bkg3=bkg2[bkg2>0]
#bins=np.logspace(np.log10(min(bkg3)),np.log10(max(bkg3)),100)
#plt.figure()
#plt.hist(bkg3,bins=bins)
#plt.xscale('log')
#plt.show()

#plim = poisson.isf(6e-5,bkg2) # this gives minimum counts to have prob equal to 6e-5
nS=1.35
plim=(nS**2+np.sqrt(nS**4+4.*nS**2*bkg2))/2. #counts needed to have a nS sigma detection
plim2=plim[plim>0]
plt.figure()
plt.hist(plim2,bins=30)
plt.yscale('log')
plt.show()

#pcts = poisson.isf(0.5,plim)
#poiss=np.random.poisson(plim-bkg2)
#pcts = poisson.isf(0.5,poiss)  # this gives counts that 50% of times make my source detected -- CHECK THIS

#plim2=plim[plim>0]
#pcts2=pcts[pcts>0]
#plt.figure()
#plt.plot(plim,pcts,'b.')
#plt.show()



flim=plim/exp*ecf
#flim=plim/4500.*ecf
plt.figure()
plt.imshow(flim)
plt.show()
flim[exp<1.]=0. #put fluxes where expo < 1 s to zero

#flim2=plim2/exp*ecf
#flim2[exp<1.]=0. #put fluxes where expo < 1 s to zero
'''
#writes output sensitivity map
hdu1 = fits.PrimaryHDU(flim,header=backheader)
hdu1.writeto(wd+'flim_'+what+'_'+cy+'.fits',overwrite=True)

flim=flim[flim>0]
bins=np.logspace(np.log10(1e-16),np.log10(1e-13),100)

plt.figure()
(n,binedges)=np.histogram(flim,bins=bins)
binc=list((binedges[i]+binedges[i+1])/2. for i in range(len(binedges)-1))
nsum=np.cumsum(n)/float(np.max(np.sum(n)))*9.3
plt.plot(binc,nsum,'k-')
plt.xscale('log')
plt.show()

sys.exit()

#thresh=5e-3 #this looks good for XBootes (corresponds to ~5e-5 wavdetect?) / also 1e-2
if band == 'broad':
	thresh=6e-5 #(99% reliabiity threshold)
elif band == 'soft':
	thresh=1.4e-4 #(99% reliabiity threshold)
else:
	thresh=6e-5 #(99% reliabiity threshold)
	
cts=np.empty_like(bkg2)
tot=np.linspace(1.0,600.,5991)

#given background and threshold probability, find how many counts are needed to have P < threhsold
for j in range(len(bkg2)):
	print(j)
	for k in range(len(bkg2[0])):
		if bkg2[j][k] > 0.0:
			i=0
			x=tot[i]
			trial=gammainc(tot[0],bkg2[j][k])
			while trial > thresh:
				i=i+1
				x=tot[i]
				trial=gammainc(x,bkg2[j][k])
			cts[j][k]=x-bkg2[j][k]
		else:
			cts[j][k]=0

cts3=cts[cts>0]
bins=np.logspace(np.log10(min(cts3)),np.log10(max(cts3)),100)
plt.figure()
plt.hist(cts3,bins=bins)
plt.xscale('log')
plt.show()

cts=cts/exp*ecf #convert cts to fluxes
cts[exp<1.0]=0. #put fluxes where expo < 1 to zero

#writes output sensitivity map
hdu1 = fits.PrimaryHDU(cts,header=backheader)
hdu1.writeto(wd+what+'_'+band+'_sens_'+gamma+cy+'.fits',overwrite=True)

sys.exit()
'''
#sensb=fits.open(wd+'cdwfs_broad_sens_1.4Cy18.fits')
#sens2b=sensb[0].data
#sens2b=sens2b[sens2b>0]

#binsb=np.logspace(np.log10(min(sens2b)),np.log10(1e-13),100)
#ab,bb=np.histogram(sens2b,bins=binsb)
#areab=ab*scale**2
#centersb=list((binsb[i+1]+binsb[i])/2. for i in range(len(binsb)-1))
#cumb=np.cumsum(areab)

#sensc=fits.open(wd+'xbootes_broad_sens_1.4Cy3.fits')
#sens2c=sensc[0].data
#sens2c=sens2c[sens2c>0]

#binsc=np.logspace(np.log10(min(sens2c)),np.log10(1e-13),100)
#ac,bc=np.histogram(sens2c,bins=binsc)
#areac=ac*scale**2
#centersc=list((binsc[i+1]+binsc[i])/2. for i in range(len(binsc)-1))
#cumc=np.cumsum(areac)

#sensd=fits.open(wd+'xbootes_broad_sens_1.4Cy18.fits')
#sens2d=sensd[0].data
#sens2d=sens2d[sens2d>0]

#binsd=np.logspace(np.log10(min(sens2d)),np.log10(1e-13),100)
#ad,bd=np.histogram(sens2d,bins=binsd)
#aread=ad*scale**2
#centersd=list((binsd[i+1]+binsd[i])/2. for i in range(len(binsd)-1))
#cumd=np.cumsum(aread)

#cumb=np.array(cumb)
#centersb=np.array(centersb)
#cumc=np.array(cumc)
#centersc=np.array(centersc)

kenterflux,kentersky=np.genfromtxt(wd+'kenter05_sens_05-7keV.dat',unpack=True)
#simflux,simsky=np.genfromtxt(wd+'skycov.dat',unpack=True)

#hdu1 = fits.PrimaryHDU(flim,header=backheader)
#hdu1.writeto(wd+'flim_'+cy+'.fits',overwrite=True)

#file=fits.open(wd+'flim_Cy3.fits')
#flim=file[0].data

flim=flim[(flim>1e-16) & (np.isfinite(flim))]
bins=np.logspace(np.log10(1e-16),np.log10(1e-13),100)
(n,binedges)=np.histogram(flim,bins=bins)
binc=list((binedges[i]+binedges[i+1])/2. for i in range(len(binedges)-1))
nsum=np.cumsum(n)/float(np.max(np.sum(n)))*9.3

#flim2=flim2[(flim2>1e-16) & (np.isfinite(flim2))]
#bins2=np.logspace(np.log10(1e-16),np.log10(1e-13),100)
#(n2,binedges2)=np.histogram(flim2,bins=bins2)
#binc2=list((binedges2[i]+binedges2[i+1])/2. for i in range(len(binedges2)-1))
#nsum2=np.cumsum(n2)/float(np.max(np.sum(n2)))*9.3

#file2=fits.open(wd+'flim_Cy18.fits')
#flim2=file2[0].data

#flim2=flim2[(flim2>1e-16) & (np.isfinite(flim2))]

#(n2,binedges2)=np.histogram(flim2,bins=bins)
#binc2=list((binedges2[i]+binedges2[i+1])/2. for i in range(len(binedges2)-1))
#nsum2=np.cumsum(n2)/float(np.max(np.sum(n2)))*9.3

#file3=fits.open(wd+'flim_cdwfs_Cy18.fits')
#flim3=file3[0].data

#flim3=flim3[(flim3>1e-16) & (np.isfinite(flim3))]

#(n3,binedges3)=np.histogram(flim3,bins=bins)
#binc3=list((binedges3[i]+binedges3[i+1])/2. for i in range(len(binedges3)-1))
#nsum3=np.cumsum(n3)/float(np.max(np.sum(n3)))*9.3

#pcts = poisson.isf(0.5,2)-bkg2

#flim4=pcts/exp*ecf
#flim=pcts/4500.*ecf
#flim4[exp<1000.]=0. #put fluxes where expo < 1 to zero
#flim4=flim4[(flim4>1e-16) & (np.isfinite(flim4))]
#(n4,binedges4)=np.histogram(flim4,bins=bins)
#binc4=list((binedges4[i]+binedges4[i+1])/2. for i in range(len(binedges4)-1))
#nsum4=np.cumsum(n4)/float(np.max(np.sum(n4)))*9.3

f,ax=plt.subplots(1)
#ax.plot(centersb,cumb,'b-',label='CDWFS')
#ax.plot(centersc,cumc,'k--',label='Cy3 XBootes')
#ax.plot(centersd,cumd,'r--',label='Cy18 XBootes')
ax.plot(kenterflux,kentersky,'k.-',label='Kenter')
ax.plot(binc,nsum,'g--',label='My XB')
#ax.plot(binc2,nsum2,'r--',label='XB Cy18')
#ax.plot(binc3,nsum3,'b--',label='CDWFS Cy18')
#ax.plot(binc4,nsum4,'y-.',label='CDWFS 6cts')
#ax.plot(simflux,simsky,'g.-',label='Sim')
ax.set_xlabel(r'0.5-7 keV flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=13)
ax.set_ylabel(r'Area [deg$^2$]',fontsize=13)
ax.set_xscale('log')
#ax.set_yscale('log')
ax.axis([5e-16,1e-13,0,11])
ax.legend()
ax.tick_params(labelsize=14)
plt.tight_layout()
plt.show()
#plt.savefig(wd+'cdwfs_full_sensitivity.pdf',format='pdf',dpi=1000)
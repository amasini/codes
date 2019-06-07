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
#from ciao_contrib.region.check_fov import FOVFiles
import time
from scipy.stats import poisson
from scipy import interpolate
from collections import OrderedDict

wd='/Users/alberto/Desktop/XBOOTES/'

#################################
what='cdwfs'
band='broad'
gamma='1.4'
#cy='Cy3' # XBOOTES was done in Cy3, CDWFS in Cy18/Cymix
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
'''
x=np.arange(0.,11.)
y=poisson.sf(x,2.1)
a,b,c=np.polyfit(x, np.log10(y), 2)
xinterp = np.linspace(0,10,100)
log_yfit = a*xinterp**2+b*xinterp + c
yinterp=10**log_yfit
diff=abs(yinterp-6e-5)
myx=xinterp[diff==min(diff)]
print(myx)
plt.figure()
plt.plot(x,y,'ro')
plt.plot(xinterp,yinterp,'k.')
plt.axhline(y=6e-5)
plt.yscale('log')
plt.show()

sys.exit()
'''
'''
print('Doing',what,'in the',band,'band, with Gamma',gamma)
#take expo map (exposure has average exposure in 4x4 pixels in s)
expmap=fits.open(wd+'new_mosaics_detection/'+what+'_'+band+'_expomap_4reb.fits')
exp=expmap[0].data
#expmap=fits.open(wd+'data/3596/repro_new_asol/out/acisf03596_broad_expomap_4reb.fits')
#exp=expmap[0].data/16.0
exp[np.isnan(exp)]=0.0 #put nans to 0

ecf=np.zeros_like(exp,dtype=float)
ecf[exp==0.0]=0.0
if band == 'broad':
	ecf[(exp<=12000.) & (exp>0.)]= 9.325E-12 #Cy3 Gamma=1.4
	ecf[exp>12000.]=1.411E-11 #Cy18 Gamma=1.4
	pthresh=8e-5
elif band == 'soft':
	ecf[(exp<=12000.) & (exp>0.)]=4.481E-12 #Cy3 Gamma=1.4
	ecf[exp>12000.]=8.707E-12 #Cy18 Gamma=1.4
	pthresh=6e-4
else:
	ecf[(exp<=12000.) & (exp>0.)]=1.973E-11 #Cy3 Gamma=1.4
	ecf[exp>12000.]=2.022E-11 #Cy18 Gamma=1.4
	pthresh=4e-5
#plt.figure()
#plt.imshow(ecf,origin='lower')
#plt.show()
#sys.exit()

#take psfmap (in arcsec)
psfmap=fits.open(wd+'psfmaps/'+what+'_'+band+'_r90sq_4reb.fits')
#psfmap=fits.open(wd+'psfmaps/03596_r90sq_4reb.fits')
psf=psfmap[0].data

#take bkg map (bkgmap has SUMMED bkg in 4x4 pixels, like data)
bkgmap=fits.open(wd+'new_mosaics_detection/'+what+'_'+band+'_bkgmap_4reb.fits')
#bkgmap=fits.open(wd+'data/3596/repro_new_asol/out/acisf03596_broad_bkgmap_4reb.fits')
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

#print(np.max(bkg2))
#sys.exit()

maxval=20
x=np.arange(0.,maxval)
xinterp = np.linspace(0,maxval-1,100)

pcts=np.zeros_like(bkg2,dtype=float)
for i in range(bkg2.shape[0]):
	for j in range(bkg2.shape[1]):
		if bkg2[i][j] != 0.0:
			#maxval=float(int(bkg2[i][j])+2)
			#x=np.arange(0.,maxval)
			#xinterp = np.linspace(0,maxval-1,100)
			y=poisson.sf(x,bkg2[i][j])
			a,b,c=np.polyfit(x, np.log10(y), 2)

			log_yfit = a*xinterp**2+b*xinterp + c
			yinterp=10**log_yfit
			
			#plt.figure()
			#plt.plot(x,y,'ro')
			#plt.plot(xinterp,yinterp,'k.')
			#plt.axhline(y=6e-5)
			#plt.yscale('log')
			#plt.show()
			
			diff=abs(yinterp-pthresh)
			#print(xinterp[diff==min(diff)])
			#sys.exit()
			pcts[i][j]=xinterp[diff==min(diff)]-bkg2[i][j]
		else:
			pcts[i][j]=0.0

flim=pcts/exp*ecf
#flim=plim/4500.*ecf
#plt.figure()
#plt.imshow(flim)
#plt.show()
flim[exp<1.]=0. #put fluxes where expo < 1 s to zero

#writes output sensitivity map
hdu1 = fits.PrimaryHDU(flim,header=backheader)
hdu1.writeto(wd+'new_mosaics_detection/'+what+'_'+band+'_sensmap_4reb.fits',overwrite=True)

sys.exit()
'''
'''
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

#sensb=fits.open(wd+'new_mosaics_detection/sens_xbootes_soft_Cy3.fits')
#sens2b=sensb[0].data
#sens2b=sens2b*16.0

#binsb=np.logspace(np.log10(1e-16),np.log10(1e-13),100)
#ab,bb=np.histogram(sens2b,bins=binsb)
#areab=ab*scale**2
#centersb=list((binsb[i+1]+binsb[i])/2. for i in range(len(binsb)-1))
#cumb=np.cumsum(areab)

'''#####
#take expo map (exposure has average exposure in 4x4 pixels in s)
expmap=fits.open(wd+'new_mosaics_detection/cdwfs_broad_expomap_4reb.fits')
exp=expmap[0].data
exp[np.isnan(exp)]=0.0 #put nans to 0

ecf=np.zeros_like(exp,dtype=float)
ecf[exp==0.0]=0.0
ecf[(exp<=12000.) & (exp>0.)]=9.81E-12
ecf[exp>12000.]=1.329E-11

cr=sens2b/ecf # this is the countrate map

ecf[exp>12000.]=9.81E-12
sensd=cr*ecf
sens2d=sensd[sensd>0]

binsd=np.logspace(np.log10(min(sens2d)),np.log10(1e-13),100)
ad,bd=np.histogram(sens2d,bins=binsd)
aread=ad*scale**2
centersd=list((binsd[i+1]+binsd[i])/2. for i in range(len(binsd)-1))
cumd=np.cumsum(aread)
######'''



f,ax=plt.subplots(3,figsize=[6,6],sharex=True)

for field in ['xbootes','cdwfs']:
	i=0
	for band in ['broad','soft','hard']:
		
		print(field,band)
		
		sensc=fits.open(wd+'new_mosaics_detection/'+field+'_'+band+'_sensmap_4reb.fits')
		sens2c=sensc[0].data
		sens2c=sens2c[sens2c>0]
		#if band == 'soft':
		#	sens2c=sens2c*16.0
		sensc.close()
		binsc=np.logspace(np.log10(min(sens2c)),np.log10(1e-13),100)
		ac,bc=np.histogram(sens2c,bins=binsc)
		areac=ac*scale**2
		centersc=list((binsc[i+1]+binsc[i])/2. for i in range(len(binsc)-1))
		cumc=np.cumsum(areac)
	
		if field == 'xbootes':
			color ='k'
		else:
			color ='b'
			
		if band == 'broad':
			style='-'
		elif band == 'soft':
			style='--'
		else:
			style='dotted'
	
		ax[i].plot(centersc,cumc,color=color,linestyle=style,label=field)
		ax[i].set_ylabel(r'Area [deg$^2$]',fontsize=13)
		ax[i].set_xscale('log')
		ax[i].axis([1e-16,5e-14,0,11])
		ax[i].tick_params(direction='inout',which='major',top=True,right=True,length=8, labelsize=14)
		ax[i].tick_params(direction='inout',which='minor',top=True,right=True,length=4, labelsize=14)

		i=i+1

kenterflux,kentersky=np.genfromtxt(wd+'kenter05_sens_05-7keV.dat',unpack=True)
ax[0].plot(kenterflux,kentersky,'k.-',label='Kenter')
#if band == 'broad':
#	ax.set_xlabel(r'0.5-7 keV flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=13)
#elif band == 'soft':
#	ax.set_xlabel(r'0.5-2 keV flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=13)
#else:
#	ax.set_xlabel(r'2-7 keV flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=13)
ax[2].set_xlabel(r'Flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=13)
handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys())
f.subplots_adjust(hspace=0)
#plt.tight_layout()
plt.show()


'''
sensc=fits.open(wd+'new_mosaics_detection/sens_xbootes_broad_Cy3.fits')
sens2c=sensc[0].data
sens2c=sens2c[sens2c>0]

binsc=np.logspace(np.log10(min(sens2c)),np.log10(1e-13),100)
ac,bc=np.histogram(sens2c,bins=binsc)
areac=ac*scale**2
centersc=list((binsc[i+1]+binsc[i])/2. for i in range(len(binsc)-1))
cumc=np.cumsum(areac)
'''

#simflux,simsky=np.genfromtxt(wd+'skycov.dat',unpack=True)

#ax.plot(centersb,cumb,'b-',label='CDWFS')

#ax.plot(centersd,cumd,'r--',label='Cy3 CDWFS')

#ax.plot(simflux,simsky,'g.-',label='Sim')

#plt.savefig(wd+'cdwfs_full_sensitivity.pdf',format='pdf',dpi=1000)
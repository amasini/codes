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

wd='/Users/alberto/Desktop/XBOOTES/'

rebin_factor=4
ecf=1.825E-11 # gamma=1.8 
#ecf=1.847E-11 # gamma=1.4
scale=(0.492/3600.)*rebin_factor #pixel size in deg
arcsec2pix=scale*3600
print(arcsec2pix)
########!!!use psfmap.py code
'''
#take expo map
expmap=fits.open(wd+'mosaic_broad_expomap_'+str(int(rebin_factor))+'rebinned.fits')
exp=expmap[0].data
exp[np.isnan(exp)]=0.0 #put nans to 0

#take WCS system to convert from pixel position to ra,dec
w=WCS(wd+'mosaic_broad_expomap_'+str(int(rebin_factor))+'rebinned.fits')

#load all FOV files
my_obs = FOVFiles('@'+wd+'fov.lis')

psf=np.zeros_like(exp)
for i in range(exp.shape[0]):
	print(i)
	tstart=time.time()
	for j in range(exp.shape[1]):
		
		(ra_aux,dec_aux) =w.all_pix2world(i,j,1) #this task starts from 1 (origin), and converts pixel x,y to ra,dec
		
		#for each ra,dec now define in which obsids it is found
		myobs = my_obs.inside(float(ra_aux),float(dec_aux))
		
		#if the pixel is covered by the expomap (expo>0), define average r90
		if len(myobs) != 0:
			r90,expo=[],[]
			#for each obsids in which the pixel is contained,
			for kk in range(len(myobs)):
				obs=myobs[kk][36:-30]
				if len(obs) == 4:
					stem='0'+obs
				elif len(obs) == 3:
					stem='00'+obs
				elif len(obs) == 5:
					stem=obs
				#extract exposure from vignetting-corrected expomap
				expomap=wd+'data/'+obs+'/repro_new_asol/out/acisf'+stem+'_broad_expomap.fits'
				s.call('dmextract infile="'+expomap+'[bin pos=circle('+str(ra_aux)+'d,'+str(dec_aux)+'d,'+str(2*0.000277778)+'d)]" mode=h outfile=expo.fits opt=generic mode=h clobber=yes',shell=True)
				s.call('dmlist "expo.fits[cols COUNTS, AREA]" data,clean | grep -v COUNTS > expo.dat',shell=True)
				(totexpo,npix)=np.genfromtxt('expo.dat',unpack=True)
				exp_i=totexpo/npix
				expo.append(exp_i)
				
			if np.sum(expo) != 0:
			
				for k in range(len(myobs)):
					#get theta to compute r90
					obs=myobs[k][36:-30]
					if len(obs) == 4:
						stem='0'+obs
					elif len(obs) == 3:
						stem='00'+obs
					elif len(obs) == 5:
						stem=obs
					path=myobs[k][:-14]+'out/acisf'+stem+'_broad_expomap.fits'
					s.call('punlearn dmcoords',shell=True)
					s.call('dmcoords '+path+' asol=none opt=cel celfmt=deg ra='+str(ra_aux)+' dec='+str(dec_aux)+'',shell=True)
					res=s.check_output('pget dmcoords theta',shell=True)
					theta=res.splitlines()
					theta=float(theta[0])
					r90_aux=1+10*(theta/10.)**2 #WARNING, this is in arcsec
			
					r90.append(r90_aux)
		
				#at the end of the loop on the obsids, get weighted average r90 using exp time as weight
				r90_av=np.average(r90,weights=expo)
				psf[i][j]=r90_av
			else: #if the pixel has total exposure = 0
				psf[i][j]=0.0
		else: #if the pixel is not covered by any obsid
			psf[i][j]=0.0

	print((time.time()-tstart)*(exp.shape[0]-(i+1))/3600.,'hours to the end.')
#print((time.time()-tstart)/60.,'minutes for the psfsize map creation.')
hdu1 = fits.PrimaryHDU(cts,header=backheader)
hdu1.writeto(wd+'/psfsize_map_'+str(int(rebin_factor))+'rebinned.fits',overwrite=True)
sys.exit()
'''
###################
'''
#from psfmap, expomap and bkgmap compute sensitivity map

#take expo map
expmap=fits.open(wd+'mosaic_broad_expomap_4rebinned.fits')
exp=expmap[0].data
exp[np.isnan(exp)]=0.0 #put nans to 0
exp=exp/16.0

#take psfmap (in arcsec)
psfmap=fits.open(wd+'average_r90squaredmap.fits')
psf=psfmap[0].data
#psf=4.

#take bkg map
bkgmap=fits.open(wd+'mosaic_broad_bkgmap_4rebinned.fits')
bkg=bkgmap[0].data
backheader=bkgmap[0].header
bkg[np.isnan(bkg)]=0.0 #put nans to 0
bkg=bkg/16.0

#sourcearea=np.pi*(psf*arcsec2pix)**2
#bkg2=bkg*sourcearea/(arcsec2pix**2)

#no need to convert to arcsec if psfmap in arsec is used
sourcearea=np.pi*psf
### bkg is in cts/pixel; bkg2 is in cts [(cts/arcsec2)*arcsec2]
pixarea=arcsec2pix*arcsec2pix
bkg2=bkg*(sourcearea/pixarea)

#thresh=5e-3 #this looks good for XBootes (corresponds to ~5e-5 wavdetect?) / also 1e-2
thresh=1e-2
cts=np.empty_like(bkg2)
tot=np.linspace(1.0,500.,3000)

#given background and threshold probability, find how many counts are needed to have P < threhsold
for j in range(len(bkg2)):
	for k in range(len(bkg2[0])):
		print(j)
		if bkg2[j][k] > 0.0:
			i=0
			x=tot[i]
			trial=gammainc(tot[0],bkg2[j][k])
			while trial > thresh:
				i=i+1
				x=tot[i]
				trial=gammainc(x,bkg2[j][k])
			cts[j][k]=x
		else:
			cts[j][k]=0

cts=cts/exp*ecf #convert cts to fluxes
cts[exp<1.0]=0. #put fluxes where expo < 1 to zero

#writes output sensitivity map
hdu1 = fits.PrimaryHDU(cts,header=backheader)
hdu1.writeto(wd+'/new_sens.fits',overwrite=True)
'''
#plt.imshow(cts,origin='lower',norm=LogNorm())
#plt.show()

#sens=fits.open(wd+'broad_sens.fits')
#sens2=sens[0].data
#sens2=sens2[sens2>0]

sensb=fits.open(wd+'new_sens.fits')
sens2b=sensb[0].data
sens2b=sens2b[sens2b>0]
print('Here.')

(kenterflux,kentersky)=np.genfromtxt('/Users/alberto/Desktop/XBOOTES/kenter05_sens_05-7keV.dat',unpack=True)
#### MY XBOOTES
sensc=fits.open(wd+'murray_sens/new_sens.fits')
sens2c=sensc[0].data
sens2c=sens2c[sens2c>0]

binsc=np.logspace(np.log10(min(sens2c)),np.log10(1e-13),100)
ac,bc=np.histogram(sens2c,bins=binsc)
areac=ac*scale**2
centersc=list((binsc[i+1]+binsc[i])/2. for i in range(len(binsc)-1))
cumc=np.cumsum(areac)
#######
'''
bins=np.logspace(np.log10(min(sens2)),np.log10(1e-13),100)
a,b=np.histogram(sens2,bins=bins)
area=a*(0.492*4/3600.)**2
centers=list((bins[i+1]+bins[i])/2. for i in range(len(bins)-1))
cum=np.cumsum(area)
'''
binsb=np.logspace(np.log10(min(sens2b)),np.log10(1e-13),100)
ab,bb=np.histogram(sens2b,bins=binsb)
areab=ab*scale**2
centersb=list((binsb[i+1]+binsb[i])/2. for i in range(len(binsb)-1))
cumb=np.cumsum(areab)



f,ax=plt.subplots(1)
#plt.hist(sens2,bins=bins)
#ax.plot(centers,cum,'k-',label='CDWFS')
ax.plot(centersb,cumb,'b-',label='CDWFS')
ax.plot(centersc,cumc,'k--',label='XBootes')
#ax.plot(kenterflux,kentersky,'k--',label='XBootes')
ax.set_xlabel(r'0.5-7 keV flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=13)
ax.set_ylabel(r'Area [deg$^2$]',fontsize=13)
ax.set_xscale('log')
#ax.set_yscale('log')
ax.axis([5e-16,1e-13,0,11])
ax.legend()
ax.tick_params(labelsize=14)
plt.tight_layout()
#plt.show()
plt.savefig(wd+'cdwfs_full_sensitivity.pdf',format='pdf',dpi=1000)
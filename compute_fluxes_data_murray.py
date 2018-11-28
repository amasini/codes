# this code extracts fluxes for detected real sources, using r90.
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
import scipy
import scipy.stats.distributions
import subprocess as s
import sys
import time
from ciao_contrib.region.check_fov import FOVFiles

def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180.*np.pi)
    return np.sqrt(((pointa[0]-pointb[0])*xx*3600.)**2 +((pointa[1]-pointb[1])*3600.)**2)

wd='/Users/alberto/Desktop/XBOOTES/'
band='broad'
band2='05to7'
if band=='broad':
	cf=1.825E-11
elif band=='soft':
	cf=1.632E-11
elif band=='hard':
	cf=2.016E-11
t_in=time.time()

#take catalog of detected sources
dat=fits.open(wd+'data/3596/repro_new_asol/out_cat0.fits')
src_ra=dat[1].data['RA']
src_dec=dat[1].data['DEC']
src_r90=dat[1].data['AV_R90']
#src_sign=dat[1].data['SRC_SIGNIFICANCE']
#src_net=dat[1].data['NET_COUNTS']
src_net=dat[1].data['NET']
dat.close()

w=open(wd+'out_cat0.reg','w')
for i in range(len(src_ra)):
	w.write('circle('+str(src_ra[i])+'d, '+str(src_dec[i])+'d, '+str(src_r90[i])+'\") #color=yellow \n')
#	#w.write('circle('+str(ra[i])+'d, '+str(dec[i])+'d, 3\") #color=cyan \n')
w.close()
sys.exit()


#ignore detected sources with NET_CTS<cut
cut=2.0
src_ra=src_ra[src_net>=cut]
src_dec=src_dec[src_net>=cut]
#src_sign=src_sign[src_net>=cut]

print('There are '+str(len(src_ra))+' detected sources with NET_CTS >='+str(cut))
sys.exit()

out_r90,out_det,out_tot,out_back,out_cts,out_exp,out_cr,out_flux=[],[],[],[],[],[],[],[]
#for each detected source compute flux from data and bkg maps
for i in range(len(src_ra)):
	print(i+1,len(src_ra))
	#see in which obsids this source is contained 
	my_obs = FOVFiles('@'+wd+'fov.lis')
	myobs = my_obs.inside(src_ra[i], src_dec[i])
	cts,bkg,exp,r=[],[],[],[]
	for j in range(len(myobs)):
		#get theta to compute r90
		obs=myobs[j][36:-30]
		#print(obs)
		if len(obs) == 4:
			stem='0'+obs
		elif len(obs) == 3:
			stem='00'+obs
		elif len(obs) == 5:
			stem=obs
		#here, I could just use the r90 maps I created before...
		path=myobs[j][:-14]+'out/acisf'+stem+'_'+band+'_expomap.fits'
		s.call('punlearn dmcoords',shell=True)
		s.call('dmcoords '+path+' asol=none opt=cel celfmt=deg ra='+str(src_ra[i])+' dec='+str(src_dec[i])+'',shell=True)
		res=s.check_output('pget dmcoords theta',shell=True)
		theta=res.splitlines()
		theta=float(theta[0])
		if band == 'hard':
			r90=1.8+10*(theta/10.)**2
		else:
			r90=1+10*(theta/10.)**2
		
		#extract counts and background using r90
		imagemap=wd+'data/'+obs+'/repro_new_asol/acisf'+stem+'_repro_'+band2+'keV.img'
		#imagemap=wd+'sim_full_bitpix-32/acisf'+stem+'_sim_poiss.fits'
		
		backmap=wd+'data/'+obs+'/repro_new_asol/out/acisf'+stem+'_'+band+'_bkgmap.fits'
		s.call('pset dmextract infile="'+imagemap+'[bin pos=circle('+str(src_ra[i])+'d,'+str(src_dec[i])+'d,'+str(r90*0.000277778)+'d)]" mode=h bkg="'+backmap+'[bin pos=circle('+str(src_ra[i])+'d,'+str(src_dec[i])+'d,'+str(r90*0.000277778)+'d)]" outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
		s.call('dmextract',shell=True)
		s.call('dmlist "counts.fits[cols COUNTS, BG_COUNTS]" data,clean | grep -v COUNTS > counts.dat', shell=True)
		(cts_i,bkg_i)=np.genfromtxt('counts.dat',unpack=True)
		cts.append(cts_i)
		bkg.append(bkg_i)
		r.append(r90)
		
		#extract exposure from vignetting-corrected expomap
		expomap=wd+'data/'+obs+'/repro_new_asol/out/acisf'+stem+'_'+band+'_expomap.fits'
		s.call('dmextract infile="'+expomap+'[bin pos=circle('+str(src_ra[i])+'d,'+str(src_dec[i])+'d,'+str(r90*0.000277778)+'d)]" mode=h outfile=expo.fits opt=generic mode=h clobber=yes',shell=True)
		s.call('dmlist "expo.fits[cols COUNTS, AREA]" data,clean | grep -v COUNTS > expo.dat',shell=True)
		(totexpo,npix)=np.genfromtxt('expo.dat',unpack=True)
		exp_i=totexpo/npix
		exp.append(exp_i)
	
	av_r90=np.sum(np.array(r)*np.array(exp))/np.sum(np.array(exp))
	prob=scipy.stats.distributions.poisson.pmf(np.sum(np.array(cts)),np.sum(np.array(bkg)))
	det=-(np.log(prob))
	if prob == 0.0:
		det=800
	net=np.sum(np.array(cts))-np.sum(np.array(bkg))
	count_rate=net*1.1/np.sum(np.array(exp))
	flux=count_rate*cf #Gamma=1.8
	#define output stuff	
	out_det.append(det)
	out_r90.append(av_r90)
	out_tot.append(np.sum(np.array(cts)))
	out_back.append(np.sum(np.array(bkg)))
	out_cts.append(net*1.1)
	out_exp.append(np.sum(np.array(exp)))
	out_cr.append(count_rate)
	out_flux.append(flux)

#write catalog
cat=Table([src_ra,src_dec,out_det,out_r90,out_tot,out_back,out_cts,out_exp,out_cr,out_flux],names=('RA','DEC','DETML','AV_R90','TOT','BKG','NET','EXP','CR','FLUX'))
cat.write(wd+'data/3596/repro_new_asol/out_cat0.fits',format='fits',overwrite=True) 

t_out=time.time()
print((t_out-t_in)/60.,'Minutes for the loop')

'''
bins=np.logspace(np.log10(1e-16),np.log10(1e-12),30)
plt.figure()
plt.hist(out_flux,bins=bins,histtype='step',color='red')
plt.xscale('log')
plt.xlabel(r'$F_{0.5-7}$ [erg cm$^{-2}$ s$^{-1}$]')
plt.tight_layout()
plt.show()

bins=np.linspace(1,301,60)
plt.figure()
plt.hist(out_cts,bins=bins,histtype='step',color='blue')
plt.xscale('log')
plt.xlabel(r'Net$_{0.5-7}$ [erg cm$^{-2}$ s$^{-1}$]')
plt.tight_layout()
plt.show()
'''

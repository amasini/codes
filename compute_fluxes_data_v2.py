# this code extracts counts, exposure, fluxes for detected real sources, using r90.
# r90 is computed with the weighted average of the r90s of all the obsids in which the
# source is contained - version 2 with fewer calls to dmextract and dmlist
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

wd="/Users/alberto/Desktop/XBOOTES/"
band='broad'
bandbkg='broad' # this is just for the soft band, use 'broad' and 'hard' for the other bands

if band=='broad':
	cf=1.825E-11
elif band=='soft':
	cf=1.632E-11
elif band=='hard':
	cf=2.016E-11
t_in=time.time()

#take catalog of detected sources - output from wavdetect
dat=fits.open(wd+'cdwfs_'+band+'_src.fits')
src_ra=dat[1].data['RA']
src_dec=dat[1].data['DEC']
src_sign=dat[1].data['SRC_SIGNIFICANCE']
src_net=dat[1].data['NET_COUNTS']
src_r=dat[1].data['R']

print(len(src_ra))
'''
mask=[]
#select only sources with physical ellipse size
for i in range(len(src_ra)):
	if src_r[i][1] != 0.:
		if src_r[i][1] >= 1.5e-6:
			mask.append(1)
			#if src_r[i][0] == src_r[i][1]:
			#	print(src_r[i])
		else:
			mask.append(0)
	else:
		mask.append(0)
dat.close()
mask=np.array(mask)
src_ra=src_ra[mask==1]
print(len(src_ra))
sys.exit()

src_dec=src_dec[mask==1]
src_sign=src_sign[mask==1]
src_net=src_net[mask==1]
src_r=src_r[mask==1]
'''
#w=open(wd+'cdwfs_broad_output_wavdetect.reg','w')
#for i in range(len(ra)):
#	w.write('circle('+str(src_ra[i])+'d, '+str(src_dec[i])+'d, 5\") #color=yellow \n')
#	#w.write('circle('+str(ra[i])+'d, '+str(dec[i])+'d, 3\") #color=cyan \n')
#w.close()
#sys.exit()

#ignore detected sources with NET_CTS<cut
#cut=0.0
#src_ra=src_ra[src_net>=cut]
#src_dec=src_dec[src_net>=cut]
#src_sign=src_sign[src_net>=cut]

#print('There are '+str(len(src_ra))+' detected sources with NET_CTS >='+str(cut))

out_r90,out_prob,out_tot,out_back,out_net,out_enetp,out_enetn,out_exp,out_cr,out_ecrp,out_ecrn,out_flux,out_efluxp,out_efluxn=[],[],[],[],[],[],[],[],[],[],[],[],[],[]
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
		path=myobs[j][:-14]+'out/acisf'+stem+'_'+bandbkg+'_expomap.fits'
		s.call('punlearn dmcoords',shell=True)
		s.call('dmcoords '+path+' asol=none opt=cel celfmt=deg ra='+str(src_ra[i])+' dec='+str(src_dec[i])+'',shell=True)
		res=s.check_output('pget dmcoords theta',shell=True)
		theta=res.splitlines()
		theta=float(theta[0])
		if band == 'hard':
			r90=1.8+10*(theta/10.)**2
		else:
			r90=1+10*(theta/10.)**2
		r.append(r90)

        
        #extract exposure from vignetting-corrected expomap
		expomap=wd+'data/'+obs+'/repro_new_asol/out/acisf'+stem+'_'+bandbkg+'_expomap.fits'
		s.call('dmextract infile="'+expomap+'[bin pos=circle('+str(src_ra[i])+'d,'+str(src_dec[i])+'d,'+str(r90*0.000277778)+'d)]" mode=h outfile=expo.fits opt=generic mode=h clobber=yes',shell=True)
		s.call('dmlist "expo.fits[cols COUNTS, AREA]" data,clean | grep -v COUNTS > expo.dat',shell=True)
		(totexpo,npix)=np.genfromtxt('expo.dat',unpack=True)
		exp_i=totexpo/npix
		exp.append(exp_i)

	av_r90=np.sum(np.array(r)*np.array(exp))/np.sum(np.array(exp))
	totexp=np.sum(np.array(exp))
	
	#extract counts and background using r90
	imagemap=wd+'mosaic_'+band+'_4rebinned.fits'
	#imagemap=wd+'sim_'+band+'/acisf'+stem+'_sim_poiss_bitpix-64.fits'
		
	backmap=wd+'mosaic_'+band+'_bkgmap_4rebinned.fits'
	s.call('dmextract infile="'+imagemap+'[bin pos=circle('+str(src_ra[i])+'d,'+str(src_dec[i])+'d,'+str(av_r90*0.000277778)+'d)]" mode=h bkg="'+backmap+'[bin pos=circle('+str(src_ra[i])+'d,'+str(src_dec[i])+'d,'+str(r90*0.000277778)+'d)]" outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
	cts=s.check_output('dmlist "counts.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
	bkg=s.check_output('dmlist "counts.fits[cols BG_COUNTS]" data,clean | grep -v BG_COUNTS',shell=True)
	#(cts_i,bkg_i)=np.genfromtxt('counts.dat',unpack=True)
	#cts,bkg=res.splitlines()
	cts,bkg=float(cts),float(bkg)
	e_cts_p=1+np.sqrt(cts+0.75) # Gehrels+86 1sigma errors
	e_cts_n=np.sqrt(cts-0.25)
	e_bkg_p=1+np.sqrt(bkg+0.75)
	e_bkg_n=np.sqrt(bkg-0.25)

	print(src_ra[i], src_dec[i], cts,bkg)
	sys.exit()
	
	prob=scipy.stats.distributions.poisson.pmf(cts,bkg)
	#det=-(np.log(prob))
	#if prob == 0.0:
	#	det=800
	net=cts-bkg
	e_net_p=np.sqrt(e_cts_p**2+e_bkg_p**2) # Propagate the errors
	e_net_n=np.sqrt(e_cts_n**2+e_bkg_n**2)

	cr=net*1.1/totexp
	e_cr_p=e_net_p*1.1/totexp # Propagate the errors
	e_cr_n=e_net_n*1.1/totexp
	flux=cr*cf #Gamma=1.8
	e_flux_p=e_cr_p*cf
	e_flux_n=e_cr_n*cf

	#define output stuff	
	out_prob.append(prob)
	out_r90.append(av_r90)
	out_tot.append(cts)
	out_back.append(bkg)
	out_net.append(net)
	out_enetp.append(e_net_p)
	out_enetn.append(e_net_n)
	out_exp.append(totexp)
	out_cr.append(cr)
	out_ecrp.append(e_cr_p)
	out_ecrn.append(e_cr_n)
	out_flux.append(flux)
	out_efluxp.append(e_flux_p)
	out_efluxn.append(e_flux_n)

#write catalog
cat=Table([src_ra,src_dec,out_prob,out_r90,out_tot,out_back,out_net,out_enetp,out_enetn,out_exp,out_cr,out_ecrp,out_ecrn,out_flux,out_efluxp,out_efluxn],names=('RA','DEC','PROB','AV_R90','TOT','BKG','NET','E_NET','e_NET','EXP','CR','E_CR','e_CR','FLUX','E_FLUX','e_FLUX'))
cat.write(wd+'cdwfs_'+band+'_cat0.fits',format='fits',overwrite=True)

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

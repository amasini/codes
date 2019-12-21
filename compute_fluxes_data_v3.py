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

t_in=time.time()
for band in ['hard']:

	#take catalog of detected sources - output from wavdetect
	dat=fits.open(wd+'new_mosaics_detection/cdwfs_'+band+'_src_exp-psf2.fits')
	src_ra=dat[1].data['RA']
	src_dec=dat[1].data['DEC']
	src_sign=dat[1].data['SRC_SIGNIFICANCE']
	src_net=dat[1].data['NET_COUNTS']
	src_r=dat[1].data['R']

	print('wavdetect found '+str(len(src_ra))+' sources in the '+band+' band.')

	#w=open(wd+'cdwfs_broad_output_wavdetect.reg','w')
	#for i in range(len(ra)):
	#	w.write('circle('+str(src_ra[i])+'d, '+str(src_dec[i])+'d, 5\") #color=yellow \n')
	#	#w.write('circle('+str(ra[i])+'d, '+str(dec[i])+'d, 3\") #color=cyan \n')
	#w.close()
	#sys.exit()
	path=wd+'new_mosaics_detection/cdwfs_'+band+'_r90_4reb.fits'
	ima=fits.open(path)
	im=ima[0].data
	ima.close()

	path=wd+'new_mosaics_detection/cdwfs_'+band+'_ecfmap_4reb.fits'
	ima=fits.open(path)
	im2=ima[0].data
	ima.close()

	out_r90,out_prob,out_tot,out_back,out_net,out_enetp,out_enetn,out_exp,out_cr,out_ecrp,out_ecrn,out_flux,out_efluxp,out_efluxn=[],[],[],[],[],[],[],[],[],[],[],[],[],[]
	#for each detected source compute flux from data and bkg maps
	for i in range(len(src_ra)):
		print(str(round(float(i+1)*100./float(len(src_ra)),2))+'%')
	
		s.call('dmcoords '+path+' asol=none opt=cel celfmt=deg ra='+str(src_ra[i])+' dec='+str(src_dec[i])+'',shell=True)
		res=s.check_output('pget dmcoords logicalx logicaly',shell=True)
		logicalx,logicaly=res.splitlines()
		logicalx,logicaly=float(logicalx),float(logicaly)
	
		# Define average R90 and ECF from the maps
		av_r90=im[int(round(logicaly)-1),int(round(logicalx)-1)]
		av_ecf=im2[int(round(logicaly)-1),int(round(logicalx)-1)]

		#extract exposure,src and bkg counts
		imagemap=wd+'new_mosaics_detection/cdwfs_'+band+'_4reb.fits'
		expomap=wd+'new_mosaics_detection/cdwfs_'+band+'_expomap_4reb.fits'	
		backmap=wd+'new_mosaics_detection/cdwfs_'+band+'_bkgmap_4reb.fits'
	
		s.call('dmextract infile="'+expomap+'[bin pos=circle('+str(src_ra[i])+'d,'+str(src_dec[i])+'d,'+str(av_r90*0.000277778)+'d)]" mode=h outfile=expo2.fits opt=generic mode=h clobber=yes',shell=True)
		s.call('dmlist "expo2.fits[cols COUNTS, AREA]" data,clean | grep -v COUNTS > expo2.dat',shell=True)
		(totexpo,npix)=np.genfromtxt('expo2.dat',unpack=True)
		av_exp=16.0*totexpo/npix # npix is the number of NATIVE CHANDRA pixels, so need to divide it by 16!
	
	
		s.call('dmextract infile="'+imagemap+'[bin pos=circle('+str(src_ra[i])+'d,'+str(src_dec[i])+'d,'+str(av_r90*0.000277778)+'d)]" mode=h bkg="'+backmap+'[bin pos=circle('+str(src_ra[i])+'d,'+str(src_dec[i])+'d,'+str(av_r90*0.000277778)+'d)]" outfile=counts2.fits opt=generic mode=h clobber=yes',shell=True)
		cts=s.check_output('dmlist "counts2.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
		bkg=s.check_output('dmlist "counts2.fits[cols BG_COUNTS]" data,clean | grep -v BG_COUNTS',shell=True)
		cts,bkg=float(cts),float(bkg)
		e_cts_p=1+np.sqrt(cts+0.75) # Gehrels+86 1sigma errors
		e_cts_n=np.sqrt(cts-0.25)
		e_bkg_p=1+np.sqrt(bkg+0.75)
		if bkg >= 0.25:
			e_bkg_n=np.sqrt(bkg-0.25)
		else:
			e_bkg_n=0	

		# This is the probability of having EXACTLY cts counts given bkg, not of having
		# AT LEAST cts counts given bkg.
		prob=scipy.stats.distributions.poisson.pmf(cts,bkg)

		net=cts-bkg
		e_net_p=np.sqrt(e_cts_p**2+e_bkg_p**2) # Propagate the errors
		e_net_n=np.sqrt(e_cts_n**2+e_bkg_n**2)

		cr=net*1.1/av_exp
		e_cr_p=e_net_p*1.1/av_exp # Propagate the errors
		e_cr_n=e_net_n*1.1/av_exp
		flux=cr*av_ecf
		e_flux_p=e_cr_p*av_ecf
		e_flux_n=e_cr_n*av_ecf

		#define output stuff	
		out_prob.append(prob)
		out_r90.append(av_r90)
		out_tot.append(cts)
		out_back.append(bkg)
		out_net.append(net)
		out_enetp.append(e_net_p)
		out_enetn.append(e_net_n)
		out_exp.append(av_exp)
		out_cr.append(cr)
		out_ecrp.append(e_cr_p)
		out_ecrn.append(e_cr_n)
		out_flux.append(flux)
		out_efluxp.append(e_flux_p)
		out_efluxn.append(e_flux_n)

	#write catalog
	cat=Table([src_ra,src_dec,out_prob,out_r90,out_tot,out_back,out_net,out_enetp,out_enetn,out_exp,out_cr,out_ecrp,out_ecrn,out_flux,out_efluxp,out_efluxn],names=('RA','DEC','PROB','AV_R90','TOT','BKG','NET','E_NET_+','E_NET_-','EXP','CR','E_CR_+','E_CR_-','FLUX','E_FLUX_+','E_FLUX_-'))
	cat.write(wd+'new_mosaics_detection/cdwfs_'+band+'_cat0_exp-psf2.fits',format='fits',overwrite=True)

t_out=time.time()
print((t_out-t_in)/60.,'Minutes for the loop')

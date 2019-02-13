import numpy as np
import sys
import subprocess as s
from astropy.io import fits
import os
import matplotlib.pyplot as plt
import random

def gauss(x,mu,sigma):
	g=np.exp(-(x-mu)**2/(2*sigma**2))
	return g
	
wd='/Users/alberto/Desktop/XBOOTES/'

band='broad'
band2='broad'
band3='05to7'

obs=np.genfromtxt(wd+'data_counts.dat',unpack=True, usecols=1,dtype='str')
diff1,diff2,exp=[],[],[]
for i in range(len(obs)):
	if len(obs[i]) == 4:
		stem='0'+obs[i]
	elif len(obs[i]) == 3:
		stem='00'+obs[i]
	elif len(obs[i]) == 5:
		stem=obs[i]
	print(obs[i])

	# Extract diffuse bkg counts from data-detected sources in F band+clusters
	s.call('dmextract infile="'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_'+band3+'keV_4rebinned.img[bin sky=region('+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_broad_src-bkg.reg)]" outfile=out2.fits opt=generic mode=h clobber=yes',shell=True)
	cts=s.check_output('dmlist "out2.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
	area=s.check_output('dmlist "out2.fits[cols AREA]" data,clean | grep -v AREA',shell=True)
	cts,area=float(cts),float(area)
	exposure=s.check_output('dmkeypar '+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img LIVETIME echo=yes',shell=True)
	exp.append(float(exposure))
	
	# Take instr bkg map
	#s.call('dmcopy \"'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_9to12keV_cl.fits\" '+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_9to12keV.img opt=image clobber=yes',shell=True)
	#image=fits.open(wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_9to12keV.img') # real 9-12 keV image
	#img=image[0].data
	#header=image[0].header
	
	#hh=fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_expomap.fits')
	#header=hh[0].header
	
	''''
	if band == 'soft':
		if os.path.isfile(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_0.5-2_bkgmap_total.fits') == True:
			image = fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_0.5-2_bkgmap_total.fits')
			out=image[0].data
			image.close()
			
			hdu = fits.PrimaryHDU(out,header=header)
			hdu.writeto(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_total.fits',overwrite=True)
			
			image = fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_0.5-2_bkgmap_instr.fits')
			out=image[0].data
			image.close()
			
			hdu = fits.PrimaryHDU(out,header=header)
			hdu.writeto(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_instr.fits',overwrite=True)
			
		else:
			image = fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_0.5-2_bkgmap_instr.fits')
			out=image[0].data
			image.close()
			
			hdu = fits.PrimaryHDU(out,header=header)
			hdu.writeto(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_instr.fits',overwrite=True)
	elif band=='hard':
		if os.path.isfile(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_hard_bkgmap_total.fits') == True:
			image = fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_hard_bkgmap_total.fits')
			out=image[0].data
			image.close()
			
			hdu = fits.PrimaryHDU(out,header=header)
			hdu.writeto(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_total.fits',overwrite=True)
			
			image = fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_hard_bkgmap_instr.fits')
			out=image[0].data
			image.close()
			
			hdu = fits.PrimaryHDU(out,header=header)
			hdu.writeto(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_instr.fits',overwrite=True)
		else:
			image = fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_hard_bkgmap_instr.fits')
			out=image[0].data
			image.close()
			
			hdu = fits.PrimaryHDU(out,header=header)
			hdu.writeto(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_instr.fits',overwrite=True)
	else:
		image=fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_total.fits')
		out=image[0].data
		image,close()
		
		hdu = fits.PrimaryHDU(out,header=header)
		hdu.writeto(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_total.fits',overwrite=True)
	'''
	#out=image[0].data
	
	#hdu = fits.PrimaryHDU(out,header=header)
	#hdu.writeto(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_total.fits',overwrite=True)
	
	
	# Rescale the counts for the total number of pixels in Chandra's FOV (16.9 arcmin^2)
	cts2=cts*(4247721./area) # This is the total bkg estimated from data
	
	# Rescale the 9-12 backmap following Hickox & Markevitch06
	if band == 'soft':
		F=0.4
	elif band == 'broad':
		F=1.52
	else:
		F=1.12
	#out=F*img
	
	'''
	#################### Can comment this
	outvalue=np.sum(out)/4247721. # cts/pixel on the instr bkgmap
	
	# Take the NON-VIGNETTING CORRECTED expomap
	dat=fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/'+band2+'_thresh.expmap')
	expo=dat[0].data
	dat.close()
	
	# Rescale expomap
	new=(expo*outvalue)/np.max(expo)
	new2=new*np.sum(out)/np.sum(new) # Renormalize in order to have the same counts as the unvignetted one
	
	# Save instr bkgmap
	hdu = fits.PrimaryHDU(new2,header=header)
	hdu.writeto(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_instr.fits',overwrite=True)
	#####################
	'''
	#print(cts2,np.sum(out))
	
	#diff=cts2-np.sum(out) # Difference between total bkg and instrumental one
	#e_diff=np.sqrt(cts2+(F*np.sum(img))) # 1sigma unc on difference
	#e_diff=np.sqrt(cts2+np.sum(out)) # 1sigma unc on difference
	
	#diff1.append(diff/e_diff)
	
	
	if band == 'soft':
		if diff > 0.:
			
			#cxb=bkg2-bkg # CXB is difference between bkg from data and just instrumental one (both not corrected for vignetting)
			mu=diff
			sigma=e_diff
			xx=np.linspace(mu-4*sigma,mu+4*sigma,1000)
			p=[]
			for ii in range(len(xx)):
				p.append(gauss(xx[ii],mu,sigma))
			p=np.array(p)
			new=p/np.sum(p)
			#plt.figure()
			#plt.plot(xx,new,'k-')
			#plt.axvline(mu)
			#plt.axvline(mu-3*sigma)
			#plt.axvline(mu+3*sigma)
			#plt.show()
			totcxb=np.random.choice(xx, p=new)
			diff2.append((diff-totcxb)/sigma)
		
			dat=fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_expomap.fits') # Take expomap VIGNETTED
			expo=dat[0].data
			dat.close()
			cxb2=diff*expo/np.max(expo) # Vignet the CXB map
			cxb3=cxb2*np.sum(totcxb)/np.sum(cxb2) # Renormalize in order to have the same counts as the unvignetted one
			hdu0 = fits.PrimaryHDU(cxb3,header=header)
			hdu0.writeto(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_cxb.fits',overwrite=True)
		
			#dat=fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_0.5-2_bkgmap_cxb.fits')
			#cxb3=dat[0].data
			#dat.close()
		
			# Sum the CXB map with the instrumental one
			totbkg=cxb3+new2
			hdu0 = fits.PrimaryHDU(totbkg,header=header)
			hdu0.writeto(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_total.fits',overwrite=True)
			
			dat=fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_total.fits')
			img=dat[0].data
			dat.close()
			diff1.append((cts2-np.sum(img))/e_diff)
		else:
			diff1.append(diff/e_diff)
	elif band == 'hard':
		if diff < 0:
			mu=diff
			sigma=e_diff
			xx=np.linspace(mu-4*sigma,mu+4*sigma,1000)
			p=[]
			for ii in range(len(xx)):
				p.append(gauss(xx[ii],mu,sigma))
			p=np.array(p)
			new=p/np.sum(p)
			#plt.figure()
			#plt.plot(xx,new,'k-')
			#plt.axvline(mu)
			#plt.axvline(mu-3*sigma)
			#plt.axvline(mu+3*sigma)
			#plt.show()
			#sys.exit()
			randomdiff=np.random.choice(xx, p=new)
			new_out=randomdiff+np.sum(out)
			e_diff2=np.sqrt(cts2+np.sum(new_out))
			diff1.append((cts2-new_out)/e_diff2)
			
			outvalue=np.sum(new_out)/4247721. # cts/pixel on the instr bkgmap
			
			# Take the NON-VIGNETTING CORRECTED expomap
			dat=fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/'+band+'_thresh.expmap')
			expo=dat[0].data
			dat.close()
	
			# Rescale expomap
			new=(expo*outvalue)/np.max(expo)
			new2=new*np.sum(new_out)/np.sum(new) # Renormalize in order to have the same counts as the unvignetted one
	
			# Save instr bkgmap
			hdu = fits.PrimaryHDU(new2,header=header)
			hdu.writeto(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_instr.fits',overwrite=True)
			
		else:
			mu=diff
			sigma=e_diff
			xx=np.linspace(mu-4*sigma,mu+4*sigma,1000)
			p=[]
			for ii in range(len(xx)):
				p.append(gauss(xx[ii],mu,sigma))
			p=np.array(p)
			new=p/np.sum(p)
			#plt.figure()
			#plt.plot(xx,new,'k-')
			#plt.axvline(mu)
			#plt.axvline(mu-3*sigma)
			#plt.axvline(mu+3*sigma)
			#plt.show()
			totcxb=np.random.choice(xx, p=new)
			diff1.append((diff-totcxb)/sigma)
		
			dat=fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_expomap.fits') # Take expomap VIGNETTED
			expo=dat[0].data
			dat.close()
			cxb2=diff*expo/np.max(expo) # Vignet the CXB map
			cxb3=cxb2*np.sum(totcxb)/np.sum(cxb2) # Renormalize in order to have the same counts as the unvignetted one
			hdu0 = fits.PrimaryHDU(cxb3,header=header)
			hdu0.writeto(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_cxb.fits',overwrite=True)
		
			#dat=fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_0.5-2_bkgmap_cxb.fits')
			#cxb3=dat[0].data
			#dat.close()
		
			# Sum the CXB map with the instrumental one
			totbkg=cxb3+new2
			hdu0 = fits.PrimaryHDU(totbkg,header=header)
			hdu0.writeto(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_total.fits',overwrite=True)
			dat.close()		
	else:
		#if os.path.isfile(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_0.5-2_bkgmap_total.fits') == True:
		#	softbkgmap = fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_0.5-2_bkgmap_total.fits')
		#else:
		softbkgmap = fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_0.5-2_bkgmap_instr.fits')
		
		softbkg=softbkgmap[0].data
		softbkghead=softbkgmap[0].header
		
		#if os.path.isfile(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_hard_bkgmap_total.fits') == True:
		#	hardbkgmap = fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_hard_bkgmap_total.fits')
		#else:
		hardbkgmap = fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_hard_bkgmap_instr.fits')
		
		hardbkg=hardbkgmap[0].data
		
		fullbkg = softbkg+hardbkg
		
		hdu = fits.PrimaryHDU(fullbkg,header=softbkghead)
		hdu.writeto(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_instr.fits',overwrite=True)
		
		diff=cts2-np.sum(fullbkg) # Difference between total bkg and instrumental one
		e_diff=np.sqrt(cts2+np.sum(fullbkg)) # 1sigma unc on difference
		diff1.append(diff/e_diff)
	
#sys.exit()
#print(np.min(diff1),np.max(diff1))
#plt.figure()
#plt.hist(factors,bins=5)
#plt.xlabel('Scaling factor from 9-12 keV to 2-7 keV')
#plt.savefig(wd+'cdwfs_bkg_'+band+'_scalingfactor.pdf',format='pdf')
#plt.show()

#mu=0
#sigma=1.5
#xx=np.linspace(mu-4*sigma,mu+4*sigma,1000)
#p=[]
#for ii in range(len(xx)):
#	p.append(gauss(xx[ii],mu,sigma))
#p=np.array(p)
#new=33*p

#print(np.min(diff2),np.max(diff2))
plt.figure()
plt.hist(diff1,bins=20)
#plt.plot(exp,diff1,'k.')
#plt.xscale('log')
plt.xlabel('Sigma deviation (Bkg_fromdata - Bkg_instr)')
plt.ylabel('N')
plt.savefig(wd+'cdwfs_bkg_'+band+'_instr.pdf',format='pdf')
#plt.show()
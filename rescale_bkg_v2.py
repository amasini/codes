import numpy as np
import sys
import subprocess as s
from astropy.io import fits
import matplotlib.pyplot as plt
import os
import random
import scipy.stats

def gauss(x,mu,sigma):
	g=np.exp(-(x-mu)**2/(2*sigma**2))
	return g
	
wd='/Users/alberto/Desktop/XBOOTES/'

band=['broad','soft','hard']
band2=['broad','0.5-2','hard']
band3=['05to7','05to2','2to7']

obs=np.genfromtxt(wd+'data_counts.dat',unpack=True, usecols=1,dtype='str')

f, ax = plt.subplots(nrows = 3, sharex = True, figsize = (6,9))
for j in range(len(band)):
	'''
	w = open(wd+'rescale_bkg_'+band[j]+'.dat','w')
	diff1,diff2=[],[]
	for i in range(len(obs)):
		if len(obs[i]) == 4:
			stem='0'+obs[i]
		elif len(obs[i]) == 3:
			stem='00'+obs[i]
		elif len(obs[i]) == 5:
			stem=obs[i]
		print(band[j],obs[i])

		# Extract diffuse bkg counts from data-detected sources in F band+clusters
		s.call('dmextract infile="'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_'+band3[j]+'keV_4rebinned.img[bin sky=region('+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_broad_src-bkg.reg)]" outfile=out2.fits opt=generic mode=h clobber=yes',shell=True)
		cts=s.check_output('dmlist "out2.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
		area=s.check_output('dmlist "out2.fits[cols AREA]" data,clean | grep -v AREA',shell=True)
		cts,area=float(cts),float(area)

		if os.path.isfile(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2[j]+'_bkgmap_total.fits') == True:
			image = fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2[j]+'_bkgmap_total.fits')
			out=image[0].data
			image.close()
		
			image = fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2[j]+'_bkgmap_instr.fits')
			out2=image[0].data
			image.close()
		
		else:
			image = fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2[j]+'_bkgmap_instr.fits')
			out=image[0].data
			image.close()
	
			out2 = out
	
		# Rescale the counts for the total number of pixels in Chandra's FOV (16.9 arcmin^2)
		cts2=cts*(4247721./area) # This is the total bkg estimated from data

		diff=cts2-np.sum(out) # Difference between data bkg and total one
		e_diff=np.sqrt(cts2+np.sum(out)) # 1sigma unc on difference
	
		diff0=cts2-np.sum(out2) # Difference between data bkg and instrumental one
		e_diff0=np.sqrt(cts2+np.sum(out2)) # 1sigma unc on difference
	
		diff1.append(diff/e_diff)
		diff2.append(diff0/e_diff0)
	
	for k in range(len(diff1)):
		w.write(str(diff1[k])+' \t '+str(diff2[k])+'\n')
	w.close()
	'''
	###
	(diff1,diff2) = np.genfromtxt(wd+'rescale_bkg_'+band[j]+'.dat',unpack=True)
	###
	
	(mu, sigma) = scipy.stats.norm.fit(diff1)
	(mu2, sigma2) = scipy.stats.norm.fit(diff2)
	
	print(band[j],mu,sigma)
	
	bins = np.linspace(-5,15,31)
	
	n, bins0, patches = ax[j].hist(diff1,bins=bins, density=True, histtype='stepfilled',linewidth=2, edgecolor='k', alpha = 0.3, label=r'($B_{\rm Data} - B_{\rm Tot}$)')
	ax[j].hist(diff2,bins=bins,density=True, histtype='stepfilled', linestyle='--', linewidth=2, edgecolor='k', alpha = 0.3, label=r'($B_{\rm Data} - B_{\rm Instr}$)')
	
	x = np.linspace(-5,15,101)
	y = scipy.stats.norm.pdf(x, mu, sigma)
	y2 = scipy.stats.norm.pdf(x, mu2, sigma2)
	ax[j].plot(x, y, 'b--', linewidth=2, label ='Gaussian Fit')
	ax[j].plot(x, y2, linestyle='--', color='orange',linewidth=2)
	ax[j].annotate(band[j].capitalize(), xy=(10,0.3))
	ax[j].set_ylabel('Fraction')
	ax[j].tick_params(direction='in', top =True, right= True)
	ax[j].axis([-5,15,0,0.55])
	ax[j].legend()
		
ax[j].set_xlabel('Sigma deviation')
#plt.legend()
plt.subplots_adjust(hspace = 0)
#plt.show()
plt.savefig(wd+'cdwfs_bkg_comparison.pdf',format='pdf')
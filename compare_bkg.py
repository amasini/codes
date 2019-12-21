# this script is old, the actual bkgmaps were done in rescale_bkg_from_blank_sky-files.py

import numpy as np
import sys
import subprocess as s
from astropy.io import fits
from ciao_contrib.region.check_fov import FOVFiles
import matplotlib.pyplot as plt
import os.path

def distance(pointa, pointb):
	xx = np.cos(pointa[1]/180*3.141592)
	return np.sqrt(((pointa[0]-pointb[0])*xx*3600)**2 +((pointa[1]-pointb[1])*3600)**2)

def gauss(x,mu,sigma):
	g=np.exp(-(x-mu)**2/(2*sigma**2))
	return g

wd="/Users/alberto/Desktop/XBOOTES/"

print('this script is old, the actual bkgmaps were done in rescale_bkg_from_blank_sky-files.py')
sys.exit()

band='soft'
band2='0.5-2' # 0.5-2 for soft
band3 = '05to2'
before = False

obs=np.genfromtxt(wd+'data_counts.dat',unpack=True, usecols=1,dtype='str')
diff1,exp=[],[]
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
	
	# Rescale the counts for the total number of pixels in Chandra's FOV (16.9 arcmin^2)
	cts2=cts*(4247721./area) # This is the total bkg estimated from data

	if band == 'soft':
		if before == False:
			if os.path.isfile(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_0.5-2_bkgmap_total.fits') == True:
				softbkgmap = fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_0.5-2_bkgmap_total.fits')
			else:
				softbkgmap = fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_0.5-2_bkgmap_instr.fits')
		else:
				softbkgmap = fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_0.5-2_bkgmap_instr.fits')
		
		bkg=softbkgmap[0].data
	elif band == 'hard':
		if before == False:
			if os.path.isfile(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_hard_bkgmap_total.fits') == True:
				hardbkgmap = fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_hard_bkgmap_total.fits')
			else:
				hardbkgmap = fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_hard_bkgmap_instr.fits')
		else:
			hardbkgmap = fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_hard_bkgmap_instr.fits')
		
		bkg=hardbkgmap[0].data
	else:
		fullbkgmap = fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_broad_bkgmap_total.fits')
		bkg=fullbkgmap[0].data
		
	diff=cts2-np.sum(bkg) # Difference between total bkg and instrumental one
	e_diff=np.sqrt(cts2+np.sum(bkg)) # 1sigma unc on difference
	diff1.append(diff/e_diff)


plt.figure()
plt.hist(diff1,bins=20)
#plt.plot(exp,diff1,'k.')
#plt.xscale('log')
plt.xlabel('Sigma deviation (Bkg_fromdata - Bkg_instr)')
plt.ylabel('N')
#plt.savefig(wd+'cdwfs_bkg_'+band+'.pdf',format='pdf')
plt.show()
sys.exit()
############
############








if band=='broad':
	R=0.13+0.27+1.12
elif band=='soft':
	R=0.13+0.27
elif band=='hard':
	R=1.12

obs=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=1,dtype='str') # List of observations
#(cts_f,cts_s,cts_h)=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=[2,3,4],dtype='float') # Total 0.5-7 keV counts
#obs=['3130']
apec_sb=np.genfromtxt(wd+'apec_cr_soft.dat',unpack=True, skip_header=1,usecols=4)

'''
cutf,cuts,cuth=1.4e-2,1e-2,3.5e-3 # These are the probability cuts in F,S,H bands at 97% rel -> 9240 srcs
# Take the catalog of detected sources (cut at 97% reliability)
cat=fits.open('/Users/alberto/Desktop/prova_cdwfs_merged_cat0.fits')
raf=cat[1].data['RA_F']
decf=cat[1].data['DEC_F']
r90f=cat[1].data['R90_F']
probf=cat[1].data['PROB_F']

ras=cat[1].data['RA_S']
decs=cat[1].data['DEC_S']
r90s=cat[1].data['R90_S']
probs=cat[1].data['PROB_S']

rah=cat[1].data['RA_H']
dech=cat[1].data['DEC_H']
r90h=cat[1].data['R90_H']
probh=cat[1].data['PROB_H']

probf[raf==0.0]=9999
probs[ras==0.0]=9999
probh[rah==0.0]=9999

# Take Kenter+05 XBOOTES Clusters
cat0=fits.open(wd+'kenter_clusters.fits')
ra_cl=cat0[1].data['_RAJ2000']
dec_cl=cat0[1].data['_DEJ2000']
r90_cl=cat0[1].data['Size']

ra_src,dec_src,r90_src=[],[],[]
for j in range(len(probf)):
	#if (probf[j] <= cutf or probs[j] <= cuts or probh[j] <= cuth):
	if band=='broad':
		if probf[j] <= cutf: # Select only srcs detected above threshold in F band
			ra_src.append(raf[j])
			dec_src.append(decf[j])
			r90_src.append(r90f[j])
	elif band=='soft':
		if probs[j] <= cuts: # Select only srcs detected above threshold in F band
			ra_src.append(ras[j])
			dec_src.append(decs[j])
			r90_src.append(r90s[j])
	elif band=='hard':
		if probh[j] <= cuth: # Select only srcs detected above threshold in F band
			ra_src.append(rah[j])
			dec_src.append(dech[j])
			r90_src.append(r90h[j])

# Add clusters to total list of sources 
for j in range(len(ra_cl)):
	ra_src.append(ra_cl[j])
	dec_src.append(dec_cl[j])
	r90_src.append(r90_cl[j])

ra_src=np.array(ra_src)
dec_src=np.array(dec_src)
r90_src=np.array(r90_src)
#'''
dev,dev_0,livetime=[],[],[]
for i in range(len(obs)):
	if len(obs[i]) == 5:
		stem=obs[i]
	elif len(obs[i]) == 4:
		stem='0'+obs[i]
	elif len(obs[i]) == 3:
		stem='00'+obs[i]
	
	sources_ra,sources_dec,sources_r90,distances=[],[],[],[]
	ra_aim=s.check_output('dmkeypar '+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img RA_PNT echo=yes',shell=True)
	dec_aim=s.check_output('dmkeypar '+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img DEC_PNT echo=yes',shell=True)
	exp=s.check_output('dmkeypar '+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img LIVETIME echo=yes',shell=True)
	ra_aim=float(ra_aim)
	dec_aim=float(dec_aim)
	
	if obs[i]=='18460':
		ra_aim=218.2957863
		dec_aim=33.043566
	
	'''#
	# Create region file with center of field of radius 6'
	w=open(wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_'+band+'_src.reg','w')
	w.write('circle('+str(ra_aim)+'d,'+str(dec_aim)+'d,360\")\n')
	w.close()
	
	my_obs = FOVFiles(wd+'data/'+obs[i]+'/repro_new_asol/fov_acisI.fits')
	
	for k in range(len(ra_src)): # We can substitute this to the direct dmextract to account for extended sources outside of FoV
		myobs = my_obs.inside(ra_src[k], dec_src[k])
		if myobs !=[]:
			sources_ra.append(ra_src[k])
			sources_dec.append(dec_src[k])
			sources_r90.append(r90_src[k])
	
	print('In obsid '+obs[i]+' there are '+str(len(sources_ra))+' sources')
	
	# Create region file with list of sources
	w=open(wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_'+band+'_bkg.reg','w')
	for k in range(len(sources_ra)):
		src=[sources_ra[k],sources_dec[k]]
		aim=[ra_aim,dec_aim]
		if distance(aim,src) <= 360.:
			radius=2.0*sources_r90[k]
			w.write('circle('+str(sources_ra[k])+'d,'+str(sources_dec[k])+'d,'+str(radius)+'\")\n')
	w.close()
	#
	# Create now appropriate regionfile
	w=open(wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_'+band+'_src-bkg.reg','w')
	file1=open(wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_'+band+'_src.reg','r')
	for line in file1:
		w.write(line)
	file1.close()
	file2=open(wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_'+band+'_bkg.reg','r')
	for line in file2:
		w.write('-'+line)
	file2.close()
	w.close()
	# DMCOPY
	#s.call('dmcopy infile="'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img[exclude sky=region('+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_'+band+'_bkg.reg)]" outfile=excised.fits clobber=yes',shell=True)
	#sys.exit()
	
	s.call('punlearn dmextract',shell=True)
	if band=='broad':
		s.call('dmextract infile="'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img[bin sky=region('+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_'+band+'_src-bkg.reg)]" outfile=out2.fits opt=generic mode=h clobber=yes',shell=True)
	elif band=='soft':
		s.call('dmextract infile="'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to2keV_4rebinned.img[bin sky=region('+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_'+band+'_src-bkg.reg)]" outfile=out2.fits opt=generic mode=h clobber=yes',shell=True)
	elif band=='hard':
		s.call('dmextract infile="'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_2to7keV_4rebinned.img[bin sky=region('+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_'+band+'_src-bkg.reg)]" outfile=out2.fits opt=generic mode=h clobber=yes',shell=True)

	cts_new=s.check_output('dmlist "out2.fits[cols SUR_BRI]" data,clean | grep -v SUR_BRI',shell=True)	
	val=float(cts_new)/float(exp) # Background counts per pixel per second from data-sources
	
	dat=fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/'+band2+'_thresh.expmap') # Open non-vignetttng corrected exposure map
	expo=dat[0].data
	header=dat[0].header
	new=expo*val
	#s.call('rm -f '+wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_new.fits',shell=True)
	hdu0 = fits.PrimaryHDU(new,header=header)
	hdu0.writeto(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_fromdata.fits',overwrite=True)
	dat.close()
	
	dat=fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_expomap.fits') # Take expomap VIGNETTED
	expo=dat[0].data
	totapec=np.sum(apec_sb[i]*expo)
	new2=apec_sb[i]*expo/np.max(expo) # Vignet the APEC map
	new3=new2*np.sum(totapec)/np.sum(new2) # Renormalize in order to have the same counts as the unvignetted one
	hdu0 = fits.PrimaryHDU(new3,header=header)
	hdu0.writeto(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_apec.fits',overwrite=True)
	dat.close()
	
	#s.call('mv '+wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap.fits '+wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_instr.fits',shell=True)
	
	# For broad and hard bands, _2.fits because the R factor was calibrated to bring the hard band into agreement with the bkg from data
	#if os.path.isfile(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_total.fits') == True:
	#	bkgmap=fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_total.fits')
	#else:
	bkgmap=fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_instr.fits')
	bkg=bkgmap[0].data
	#oldback=np.sum(bkg)
	
	apecmap=fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_apec.fits')
	apec=apecmap[0].data
	instr_plus_apec=bkg+apec
	hdu0 = fits.PrimaryHDU(instr_plus_apec,header=header)
	hdu0.writeto(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_instr+apec.fits',overwrite=True)
	#'''
	bkgmap=fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_instr+apec.fits')
	instr_plus_apec=bkgmap[0].data
	oldback=np.sum(instr_plus_apec)
	
	bkgmap2=fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_fromdata.fits')
	bkg2=bkgmap2[0].data
	bkgheader=bkgmap2[0].header
	bkgmap2.close()
	newback=np.sum(bkg2)
	# Percentage deviation
	x=(newback-oldback)/oldback
	dev_0.append(x*100)
	#print(obs[i],oldback,newback,x*100)
	
	# Sigma deviation
	x=(newback-oldback)/np.sqrt(newback+(R*oldback))
	dev.append(x)
	livetime.append(float(exp)/1e3)
	print(obs[i],oldback,newback,x)
	
	'''
	if x > 3: # Deviation is more than +3sigma
		cxb=bkg2-bkg # CXB is difference between bkg from data and just instrumental one (both not corrected for vignetting)
		#totcxb=np.random.uniform(np.sum(cxb)-3*np.sqrt(newback+(R*oldback)),np.sum(cxb)+3*np.sqrt(newback+(R*oldback)))

		mu=np.sum(cxb)
		sigma=np.sqrt(newback+(R*oldback))
		xx=np.linspace(mu-4*sigma,mu+4*sigma,1000)
		p=[]
		for ii in range(len(xx)):
			p.append(gauss(xx[ii],mu,sigma))
		p=np.array(p)
		new=p/np.sum(p)
		totcxb=np.random.choice(xx, p=new)

		#print(np.sum(cxb)-3*np.sqrt(newback),np.sum(cxb)+3*np.sqrt(newback),totcxb)
		dat=fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_expomap.fits') # Take expomap VIGNETTED
		expo=dat[0].data
		cxb2=cxb*expo/np.max(expo) # Vignet the CXB map
		cxb3=cxb2*np.sum(totcxb)/np.sum(cxb2) # Renormalize in order to have the same counts as the unvignetted one
		hdu0 = fits.PrimaryHDU(cxb3,header=bkgheader)
		hdu0.writeto(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_cxb.fits',overwrite=True)
		
		# Sum the CXB map with the instrumental one
		totbkg=cxb3+bkg
		hdu0 = fits.PrimaryHDU(totbkg,header=bkgheader)
		hdu0.writeto(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_total.fits',overwrite=True)
		dat.close()
	'''
#sys.exit()
print(np.mean(dev),np.median(dev))
plt.figure()
plt.hist(dev,bins=20)
plt.show()

dev=np.array(dev)
livetime=np.array(livetime)
xbootes=dev[livetime<6.0]
print('xbootes',np.mean(xbootes),np.median(xbootes))
plt.figure()
plt.hist(xbootes,bins=10)
plt.show()

midblock=dev[(livetime<20.0) & (livetime>=6.0)]
print('midblock',np.mean(midblock),np.median(midblock))
cdwfs=dev[livetime>20.0]
print('cdwfs',np.mean(cdwfs),np.median(cdwfs))
plt.figure()
plt.hist(cdwfs,bins=10)
plt.show()

f,(ax1,ax2)=plt.subplots(2,1,sharex=True)
ax1.plot(livetime,dev,'k.')
#plt.plot([4.5,8,30],[8.48,3.33,1.01],'go',ms=15)
ax1.axhline(y=0)
ax1.fill_between([3,200],y1=-3,y2=3)
#ax1.xlabel('Exposure time [ks]')
ax1.set_ylabel('New bkg deviation [sigma]')
ax1.set_xscale('log')

ax2.plot(livetime,dev_0,'k.')
#plt.plot([4.5,8,30],[8.48,3.33,1.01],'go',ms=15)
ax2.axhline(y=0)
ax2.set_xlabel('Exposure time [ks]')
ax2.set_ylabel('New bkg deviation [%]')
ax2.set_xscale('log')

plt.tight_layout()
plt.subplots_adjust(hspace=0.2)
plt.show()
#plt.savefig(wd+'tot_instr_bkg_dev_'+band+'.png',format='png')

'''
for k in range(len(ra_src)):
    # Compute distance array for all sources
    
    distances.append(distance(src,aim))

# Mask sources within 15' from the aimpoint to loop on
distances=np.array(distances)
ra_src2=ra_src[distances <= 900.]
dec_src2=dec_src[distances <= 900.]
r90_src2=r90_src[distances <= 900.]
    
for k in range(len(ra_src2)):
        	#radius=2.0*r90_src2[k]
        	#s.call('dmextract infile="'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img[bin pos=circle('+str(ra_src2[k])+'d,'+str(dec_src2[k])+'d,'+str(radius*0.000277778)+'d)]" outfile=out2.fits opt=generic mode=h clobber=yes',shell=True)
        	#src_cts=s.check_output('dmlist "out2.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
        	#if float(src_cts) > 0.0:
'''

# Code to extract background spectra from data and stowed background in CDWFS
import numpy as np
import sys
import subprocess as s
from astropy.io import fits
import os
import matplotlib.pyplot as plt
from ciao_contrib.region.check_fov import FOVFiles

def distance(pointa, pointb):
	xx = np.cos(pointa[1]/180*3.141592)
	return np.sqrt(((pointa[0]-pointb[0])*xx*3600)**2 +((pointa[1]-pointb[1])*3600)**2)
	
wd='/Users/alberto/Desktop/XBOOTES/'

obsid=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=1,dtype='str')

#obs='3130'
#stem='03130'

##########################################################################################
'''/// PART 1: USE CIAO TOOLS TO EXTRACT DATA AND STOWED SPECTRA ///'''
for obs in obsid:
	if len(obs) == 4:
		stem='0'+obs
	elif len(obs) == 3:
		stem='00'+obs
	elif len(obs) == 5:
		stem=obs
	print(obs)
	if obs != '3130':
		### REPROJECT THE STOWED BACKGROUND
		s.call('reproject_events infile='+wd+'stowed_2005.fits outfile='+wd+'data/'+obs+'/repro_new_asol/'+stem+'_stowed.fits aspect='+wd+'data/'+obs+'/repro_new_asol/acisf'+stem+'_new_asol1.fits match='+wd+'data/'+obs+'/repro_new_asol/acisf'+stem+'_repro_evt2.fits random=0 clobber=yes',shell=True)

		### CHECK CALIBRATION BETWEEN DATA AND STOWED BACK FILES
		gaina=s.check_output('dmkeypar '+wd+'data/'+obs+'/repro_new_asol/acisf'+stem+'_repro_evt2.fits GAINFILE echo+',shell=True)
		gainb=s.check_output('dmkeypar '+wd+'data/'+obs+'/repro_new_asol/'+stem+'_stowed.fits GAINFILE echo+',shell=True)
		if gaina[:-10] != gainb[:-10]:
			print('Check GAINFILE calibration for obsid',obs)
			sys.exit()

		### FILTERING APPROPRIATE FOR VFAINT MODE
		s.call('dmcopy "'+wd+'data/'+obs+'/repro_new_asol/'+stem+'_stowed.fits[status=0]" '+wd+'data/'+obs+'/repro_new_asol/'+stem+'_stowed_VFcl.fits clobber=yes',shell=True)

		#----------------------------------------------------------------------------------------#
		### CREATE APPROPRIATE REGIONFILES
		# Take Kenter+05 XBOOTES Clusters
		cat0=fits.open(wd+'kenter_clusters.fits')
		ra_cl=cat0[1].data['_RAJ2000']
		dec_cl=cat0[1].data['_DEJ2000']
		r90_cl=cat0[1].data['Size']
		cat0.close()

		cdwfs=fits.open(wd+'cdwfs_merged_cat1.fits')
		ra=cdwfs[1].data['RA']
		dec=cdwfs[1].data['DEC']
		r90f=cdwfs[1].data['R90_F']
		r90s=cdwfs[1].data['R90_S']
		r90h=cdwfs[1].data['R90_H']
		cdwfs.close()

		r90=r90f
		r90[r90==0.0]=r90s[r90==0.0]
		r90[r90==0.0]=r90h[r90==0.0]

		ra_src,dec_src,r90_src=[],[],[]
		for j in range(len(ra)):
			ra_src.append(ra[j])
			dec_src.append(dec[j])
			r90_src.append(r90[j])

		# Add clusters to total list of sources 
		for j in range(len(ra_cl)):
			ra_src.append(ra_cl[j])
			dec_src.append(dec_cl[j])
			r90_src.append(r90_cl[j])

		ra_src=np.array(ra_src)
		dec_src=np.array(dec_src)
		r90_src=np.array(r90_src)

		ra_aim=s.check_output('dmkeypar '+wd+'data/'+obs+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img RA_PNT echo=yes',shell=True)
		dec_aim=s.check_output('dmkeypar '+wd+'data/'+obs+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img DEC_PNT echo=yes',shell=True)
		ra_aim=float(ra_aim)
		dec_aim=float(dec_aim)
	
		if obs == '18460':
			ra_aim=218.2957863
			dec_aim=33.043566

		sources_ra,sources_dec,sources_r90=[],[],[]
		my_obs = FOVFiles(wd+'data/'+obs+'/repro_new_asol/fov_acisI.fits')
		for k in range(len(ra_src)):
			myobs = my_obs.inside(ra_src[k], dec_src[k])
			if myobs !=[]:
				sources_ra.append(ra_src[k])
				sources_dec.append(dec_src[k])
				sources_r90.append(2.0*r90_src[k])
		print('In obsid '+obs+' there are '+str(len(sources_ra))+' sources')

		w=open(wd+'data/'+obs+'/repro_new_asol/acisf'+stem+'_broad_bkg.reg','w')
		for l in range(len(sources_ra)):
			src=[sources_ra[l],sources_dec[l]]
			aim=[ra_aim,dec_aim]
			if distance(aim,src) <= 360.:
				w.write('circle('+str(sources_ra[l])+'d,'+str(sources_dec[l])+'d,'+str(sources_r90[l])+'\")\n')
		w.close()

		w=open(wd+'data/'+obs+'/repro_new_asol/acisf'+stem+'_broad_src-bkg.reg','w')
		file1=open(wd+'data/'+obs+'/repro_new_asol/acisf'+stem+'_broad_src.reg','r')
		for line in file1:
			w.write(line)
		file1.close()
		file2=open(wd+'data/'+obs+'/repro_new_asol/acisf'+stem+'_broad_bkg.reg','r')
		for line in file2:
			w.write('-'+line)
		file2.close()
		w.close()
		#----------------------------------------------------------------------------------------#

		### CREATE 0.3-12 keV FILE FILTERED FOR BACKGROUND FLARES
		s.call('dmcopy "'+wd+'/data/'+obs+'/repro_new_asol/acisf'+stem+'_repro_evt2.fits[@'+wd+'/data/'+obs+'/repro_new_asol/acisf'+stem+'_lc_05to7keV.gti][ccd_id<4][energy=300:12000]" '+wd+'/data/'+obs+'/repro_new_asol/acisf'+stem+'_repro_03to12keV_cl.fits clobber=yes',shell=True)

		### EXTRACT DATA SPECTRUM WITH SPECEXTRACT
		s.call('specextract "'+wd+'data/'+obs+'/repro_new_asol/acisf'+stem+'_repro_03to12keV_cl.fits[sky=region('+wd+'data/'+obs+'/repro_new_asol/acisf'+stem+'_broad_src-bkg.reg)]" '+wd+'data/'+obs+'/repro_new_asol/extract/spec_broad energy=0.3:12.0:0.01 grouptype=NONE binspec=NONE clobber=yes',shell=True)

		### EXTRACT STOWED SPECTRUM
		s.call('dmextract "'+wd+'data/'+obs+'/repro_new_asol/'+stem+'_stowed_VFcl.fits[sky=region('+wd+'data/'+obs+'/repro_new_asol/acisf'+stem+'_broad_src-bkg.reg)][bin pi=1:1024:1]" '+wd+'data/'+obs+'/repro_new_asol/extract/'+stem+'_stowed.pha wmap="det=8" clobber=yes',shell=True)

		### MAKE THE RMF FOR THE STOWED SPECTRUM
		s.call('mkacisrmf infile=CALDB outfile='+wd+'data/'+obs+'/repro_new_asol/extract/'+stem+'_stowed.rmf energy=0.3:12.0:0.01 channel=1:1024:1 chantype=PI wmap='+wd+'data/'+obs+'/repro_new_asol/extract/'+stem+'_stowed.pha gain=CALDB asolfile='+wd+'data/'+obs+'/repro_new_asol/acisf'+stem+'_new_asol1.fits clobber=yes',shell=True)
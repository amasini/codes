# This script is deprecated from 05/23/18. It has been implemented into repro_new_obsid.py instead. 
import numpy as np
import subprocess as s
import sys
import os

wd='/Users/alberto/Desktop/XBOOTES/'

obs=np.genfromtxt(wd+'data_counts.dat',unpack=True, usecols=1,dtype='str')

for i in range(len(obs)):
	if len(obs[i]) == 4:
		stem='0'+obs[i]
	elif len(obs[i]) == 3:
		stem='00'+obs[i]
	elif len(obs[i]) == 5:
		stem=obs[i]
	print(obs[i])
	
	#s.call('punlearn dmcopy',shell=True)
	#s.call('dmcopy \"'+wd+'data/'+obs[i]+'/repro_new_asol/expo_broad_flux.img[ccd_id=0:3][bin x=::4,y=::4]\" '+wd+'data/'+obs[i]+'/repro_new_asol/expo_broad_flux_4rebinned.img opt=image clobber=yes',shell=True)
	#s.call('dmcopy \"'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_evt2.fits[ccd_id=0:3][energy=500:7000][bin x=::4,y=::4]\" '+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img opt=image clobber=yes',shell=True)
	
	#if os.path.isfile(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_hard_bkgmap_total.fits') == True:
	#	s.call('dmcopy \"'+wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_hard_bkgmap_total.fits[bin (x,y)=::4]\" outfile='+wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_hard_bkgmap_4reb.fits clobber=yes',shell=True)
	#else:
	#	s.call('dmcopy \"'+wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_hard_bkgmap_instr.fits[bin (x,y)=::4]\" outfile='+wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_hard_bkgmap_4reb.fits clobber=yes',shell=True)
	
	s.call('dmcopy \"'+wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_0.5-2_bkgmap_instr.fits[bin (x,y)=::4]\" outfile='+wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_soft_bkgmap_instr_4reb.fits clobber=yes',shell=True)
	
	#s.call('dmcopy \"'+wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_hard_expomap.fits[bin (x,y)=::4]\" outfile='+wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_hard_expomap_4reb.fits clobber=yes',shell=True)
	#s.call('dmcopy \"'+wd+'/psfmaps/'+stem+'_hard_r90squared-x-expo.fits[bin (x,y)=::4]\" outfile='+wd+'/psfmaps/'+stem+'_hard_r90sq-x-exp_4reb.fits clobber=yes',shell=True)
	
	#s.call('dmcopy \"'+wd+'/ecfmaps/'+stem+'_broad_ecf-x-expo.fits[bin (x,y)=::4]\" outfile='+wd+'/ecfmaps/'+stem+'_broad_ecf-x-exp_4reb.fits clobber=yes',shell=True)
	#s.call('dmcopy \"'+wd+'/ecfmaps/'+stem+'_soft_ecf-x-expo.fits[bin (x,y)=::4]\" outfile='+wd+'/ecfmaps/'+stem+'_soft_ecf-x-exp_4reb.fits clobber=yes',shell=True)
	#s.call('dmcopy \"'+wd+'/ecfmaps/'+stem+'_hard_ecf-x-expo.fits[bin (x,y)=::4]\" outfile='+wd+'/ecfmaps/'+stem+'_hard_ecf-x-exp_4reb.fits clobber=yes',shell=True)
	
	#for band in ['soft','hard']:
	#	s.call('dmcopy \"'+wd+'/sim_all_FINAL/acisf'+stem+'_'+band+'_sim.fits[bin (x,y)=::4]\" outfile='+wd+'/sim_all_FINAL/acisf'+stem+'_'+band+'_sim_4reb.fits clobber=yes',shell=True)
	
	#s.call('dmcopy \"'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_cl.fits[energy=500:2000]\" '+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to2keV_cl.fits clobber=yes',shell=True)
	#s.call('dmcopy \"'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_cl.fits[energy=500:2000][bin x=::4,y=::4]\" '+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to2keV_4rebinned.img opt=image clobber=yes',shell=True)
	
	#s.call('dmcopy \"'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_cl.fits[energy=2000:7000]\" '+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_2to7keV_cl.fits clobber=yes',shell=True)
	#s.call('dmcopy \"'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_cl.fits[energy=2000:7000][bin x=::4,y=::4]\" '+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_2to7keV_4rebinned.img opt=image clobber=yes',shell=True)
# THIS SCRIPT IS SIMPLY A WRAPPER TO WRITE XIMAGE COMMAND FILES AND CALL XIMAGE
import numpy as np
import sys
import os
import subprocess as s
import time

start_time=time.time()
wd='/Users/alberto/Desktop/XBOOTES/'

#(dirs,obs)=np.genfromtxt(wd+'murray_sens/murray_data_counts.dat',unpack=True,usecols=[0,1],dtype=str)
(dirs,obs)=np.genfromtxt(wd+'/data_counts.dat',unpack=True,usecols=[0,1],dtype=str)

'''
#obs=obs[dirs=='murray05']
w=open(wd+'/xbootes_regfile.reg','w')
for i in range(len(obs)):
    if len(obs[i]) == 4:
        stem='0'+obs[i]
    elif len(obs[i]) == 3:
        stem='00'+obs[i]
    elif len(obs[i]) == 5:
        stem=obs[i]
    ra=s.check_output('dmkeypar '+wd+dirs[i]+'/'+obs[i]+'/repro/acisf'+stem+'_repro_evt2.fits RA_PNT echo=yes',shell=True)
    dec=s.check_output('dmkeypar '+wd+dirs[i]+'/'+obs[i]+'/repro/acisf'+stem+'_repro_evt2.fits DEC_PNT echo=yes',shell=True)    
    ra,dec=float(ra),float(dec)
    w.write('circle('+str(ra)+'d,'+str(dec)+'d,3\") #text={'+obs[i]+'}\n')
w.close()
'''
for band in ['05to7','05to2','2to7']:
	if band =='05to7':
		band2='broad'
	elif band=='05to2':
		band2='0.5-2'
	else:
		band2='hard'
	#start from 4221 in a 7200x7200 square
	#w=open(wd+'commands_r90.xco','w')
	w2=open(wd+'commands_bkg_instr.xco','w')
	#w3=open(wd+'commands_exp.xco','w')
	
	#w.write('read/fits/size=7200 '+wd+'data/4221/repro_new_asol/acisf04221_repro_'+band+'keV_4rebinned.img\n')
	#w.write('read/fits/size=7200/rebin=4 '+wd+'sim_all/acisf04221_'+band+'_sim_poiss.fits\n')
	#w.write('read/fits/size=7200/rebin=4 '+wd+'psfmaps/04221_r90-x-expo.fits\n')
	'''
	if band2 == 'broad':
		w2.write('read/fits/size=7200/rebin=4 '+wd+'data/4221/repro_new_asol/out/acisf04221_'+band2+'_bkgmap_total.fits\n')
	else:
		if os.path.isfile(wd+'data/4221/repro_new_asol/out/acisf04221_'+band2+'_bkgmap_total.fits') == True:
			w2.write('read/fits/size=7200/rebin=4 '+wd+'data/4221/repro_new_asol/out/acisf04221_'+band2+'_bkgmap_total.fits\n')
		else:
			w2.write('read/fits/size=7200/rebin=4 '+wd+'data/4221/repro_new_asol/out/acisf04221_'+band2+'_bkgmap_instr.fits\n')
	'''
	w2.write('read/fits/size=7200/rebin=4 '+wd+'data/4221/repro_new_asol/out/acisf04221_'+band2+'_bkgmap_instr.fits\n')
	#w3.write('read/fits/size=7200/rebin=4 '+wd+'data/4221/repro_new_asol/out/acisf04221_'+band2+'_expomap.fits\n')
	#w3.write('rebin/mode=0 4 \n')
	#w.write('read/fits/size=7200 '+wd+'data/4221/repro_new_asol/expo_0.5-2.0_flux.img\n')
	#w2.write('read/fits/size=7200 '+wd+'data/4221/repro_new_asol/expo_2.0-4.5_flux.img\n')
	#w3.write('read/fits/size=7200 '+wd+'data/4221/repro_new_asol/expo_4.5-7.0_flux.img\n')
	
	#w.write('save_image\n')
	w2.write('save_image\n')
	#w3.write('save_image\n')
	
	for i in range(len(obs)):
		if len(obs[i]) == 4:
			stem='0'+obs[i]
		elif len(obs[i]) == 3:
			stem='00'+obs[i]
		elif len(obs[i]) == 5:
			stem=obs[i]
		if obs[i]!='4221':
			#w.write('read/fits '+wd+dirs[i]+'/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_'+band+'keV_4rebinned.img\n')			
			#w.write('read/fits/rebin=4 '+wd+'sim_all/acisf'+stem+'_'+band+'_sim_poiss.fits\n')
			#w.write('read/fits/rebin=4 '+wd+'psfmaps/'+stem+'_r90-x-expo.fits\n')
			'''
			if band2 == 'broad':
				w2.write('read/fits/rebin=4 '+wd+dirs[i]+'/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_total.fits\n')
			else:
				if os.path.isfile(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_total.fits') == True:
					w2.write('read/fits/size=7200/rebin=4 '+wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_total.fits\n')
				else:
					w2.write('read/fits/size=7200/rebin=4 '+wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_instr.fits\n')
			'''
			w2.write('read/fits/rebin=4 '+wd+dirs[i]+'/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_instr.fits\n')
			#w3.write('read/fits/rebin=4 '+wd+dirs[i]+'/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_expomap.fits\n')
			#w3.write('rebin/mode=0 4 \n')
			#w.write('read/fits '+wd+dirs[i]+'/'+obs[i]+'/repro_new_asol/expo_0.5-2.0_flux.img\n')
			#w2.write('read/fits '+wd+dirs[i]+'/'+obs[i]+'/repro_new_asol/expo_2.0-4.5_flux.img\n')
			#w3.write('read/fits '+wd+dirs[i]+'/'+obs[i]+'/repro_new_asol/expo_4.5-7.0_flux.img\n')
			
			#w.write('sum_image\n')
			w2.write('sum_image\n')
			#w3.write('sum_image\n')
			
			#w.write('save_image\n')
			w2.write('save_image\n')
			#w3.write('save_image\n')
			

	if band =='05to2':
		#w.write('write/fits '+wd+'mosaic_soft_4rebinned.fits\n')
		w2.write('write/fits '+wd+'mosaic_soft_bkgmap_4rebinned_instr.fits\n')
		#w3.write('write/fits '+wd+'mosaic_soft_expomap_4rebinned.fits\n')
	elif band =='2to7':
		#pass
		#w.write('write/fits '+wd+'mosaic_hard_4rebinned.fits\n')
		w2.write('write/fits '+wd+'mosaic_hard_bkgmap_4rebinned_instr.fits\n')
		#w3.write('write/fits '+wd+'mosaic_hard_expomap_4rebinned.fits\n')
	else:
		#pass
		w2.write('write/fits '+wd+'mosaic_broad_bkgmap_4rebinned_instr.fits\n')
		#w.write('write/fits '+wd+'mosaic_broad_4rebinned.fits\n')
		#w.write('write/fits '+wd+'mosaic_broad_sim_poiss_4rebinned.fits\n')
		#w2.write('write/fits '+wd+'/murray_sens/mosaic_broad_bkgmap_4rebinned_murray.fits\n')
		#w3.write('write/fits '+wd+'/mosaic_broad_expomap_4rebinned.fits\n')
		
		#w.write('write/fits '+wd+'/cdwfs_soft_4rebin.fits\n')
		#w2.write('write/fits '+wd+'/cdwfs_medium_4rebin.fits\n')
		#w3.write('write/fits '+wd+'/cdwfs_hard_4rebin.fits\n')
	
	#w.write('exit\n')
	w2.write('exit\n')
	#w3.write('exit\n')
	
	#w.close()
	w2.close()
	#w3.close()
	#s.call('ximage @'+wd+'commands_r90.xco',shell=True)
	s.call('ximage @'+wd+'commands_bkg_instr.xco',shell=True)
	#s.call('ximage @'+wd+'commands_exp.xco',shell=True)

	elapsed_time=time.time()-start_time
	print(float(elapsed_time)/60.)

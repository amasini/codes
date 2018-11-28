import numpy as np
import sys
import os
import subprocess as s
import time

wd='/Users/alberto/Desktop/XBOOTES/'
#################
dirs='data'    #
obsid=['19780','19781']   #
#################
'''
#add exposure corrected img to mosaic
w=open(wd+'commands.xco','w')
w.write('read/fits/size=7200 '+wd+'mosaic_broad_4rebinned.fits\n')
w.write('save_image\n')
for obs in obsid:
	if len(obs) == 4:
		stem='0'+obs
	elif len(obs) == 3:
		stem='00'+obs
	elif len(obs) == 5:
		stem=obs
	w.write('read/fits/rebin=4 '+wd+dirs+'/'+obs+'/repro_new_asol/expo_broad_flux.img\n')
	w.write('sum_image\n')
	w.write('save_image\n')
w.write('write/fits '+wd+'mosaic_broad_4rebinned2.fits\n')
w.write('exit\n')
w.close()
s.call('ximage @'+wd+'commands.xco',shell=True)

#add NOT exposure corrected img to mosaic
w=open(wd+'commands.xco','w')
w.write('read/fits/size=7200 '+wd+'mosaic_broad_4rebinned_not-expo-corr.fits\n')
w.write('save_image\n')
for obs in obsid:
	if len(obs) == 4:
		stem='0'+obs
	elif len(obs) == 3:
		stem='00'+obs
	elif len(obs) == 5:
		stem=obs
	w.write('read/fits/size=1000 '+wd+dirs+'/'+obs+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img\n')
	w.write('sum_image\n')
	w.write('save_image\n')
w.write('write/fits '+wd+'mosaic_broad_4rebinned_not-expo-corr2.fits\n')
w.write('exit\n')
w.close()
s.call('ximage @'+wd+'commands.xco',shell=True)
'''
#add sim to sim mosaic
w=open(wd+'commands.xco','w')
w.write('read/fits/size=7200 '+wd+'mosaic_broad_sim_poiss_4rebinned.fits\n')
w.write('save_image\n')
for obs in obsid:
	if len(obs) == 4:
		stem='0'+obs
	elif len(obs) == 3:
		stem='00'+obs
	elif len(obs) == 5:
		stem=obs
	w.write('read/fits/rebin=4 '+wd+'/sim_full/acisf'+stem+'_sim_poiss.fits\n')
	w.write('sum_image\n')
	w.write('save_image\n')
w.write('write/fits '+wd+'mosaic_broad_sim_poiss_4rebinned2.fits\n')
w.write('exit\n')
w.close()
s.call('ximage @'+wd+'commands.xco',shell=True)
'''
#add expomap to expomap mosaic
w=open(wd+'commands.xco','w')
w.write('read/fits/size=7200 '+wd+'mosaic_broad_expomap_4rebinned.fits\n')
w.write('save_image\n')
for obs in obsid:
	if len(obs) == 4:
		stem='0'+obs
	elif len(obs) == 3:
		stem='00'+obs
	elif len(obs) == 5:
		stem=obs
	w.write('read/fits/rebin=4 '+wd+dirs+'/'+obs+'/repro_new_asol/out/acisf'+stem+'_broad_expomap.fits\n')
	w.write('sum_image\n')
	w.write('save_image\n')
w.write('write/fits '+wd+'mosaic_broad_expomap_4rebinned2.fits\n')
w.write('exit\n')
w.close()
s.call('ximage @'+wd+'commands.xco',shell=True)
'''

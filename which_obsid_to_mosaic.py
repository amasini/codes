#define a center for the chunk of mosaic over which run wavdetect
import numpy as np
import subprocess as s
import sys
import time
import os

def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*xx)**2 +((pointa[1]-pointb[1]))**2)

wd='/Users/alberto/Desktop/XBOOTES/'

ra_c=[218.9132491,217.1703769,218.6973158,217.0839299,218.0404204,216.7921841]
dec_c=[35.20442382,35.25246411,33.93039948,33.92450157,32.83570264,32.82283117]
radius=[0.85,0.85,0.85,0.85,0.6,0.6]

(dirs,obsid)=np.genfromtxt(wd+'murray_sens/murray_data_counts.dat',unpack=True,usecols=[0,1],dtype='str')
'''
w4=open(wd+'/xbootes_regfile.reg','w')
for i in range(len(obsid)):
    if len(obsid[i]) == 4:
        stem='0'+obsid[i]
    elif len(obsid[i]) == 3:
        stem='00'+obsid[i]
    elif len(obsid[i]) == 5:
        stem=obsid[i]
    ra=s.check_output('dmkeypar '+wd+dirs[i]+'/'+obsid[i]+'/repro/acisf'+stem+'_repro_evt2.fits RA_PNT echo=yes',shell=True)
    dec=s.check_output('dmkeypar '+wd+dirs[i]+'/'+obsid[i]+'/repro/acisf'+stem+'_repro_evt2.fits DEC_PNT echo=yes',shell=True)    
    ra,dec=float(ra),float(dec)
    w4.write('circle('+str(ra)+'d,'+str(dec)+'d,3\") #text={'+obsid[i]+'}\n')
w4.close()
'''
#band=['05to2','2to7']
band=['05to7']
band2=['broad']
#band2=['soft','hard']
for chunk_number in range(0,len(ra_c)):
#for chunk_number in range(0,6):
	
	list_to_mosaic=[]
	center=[ra_c[chunk_number],dec_c[chunk_number]]
	'''
	#take the coordinates
	for i in range(len(obsid)):
		if len(obsid[i]) == 4:
			stem='0'+obsid[i]
		elif len(obsid[i]) == 3:
			stem='00'+obsid[i]
		elif len(obsid[i]) == 5:
			stem=obsid[i]
		ra=s.check_output('dmkeypar '+wd+dirs[i]+'/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_evt2.fits RA_PNT echo=yes',shell=True)
		dec=s.check_output('dmkeypar '+wd+dirs[i]+'/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_evt2.fits DEC_PNT echo=yes',shell=True)
		pnt=[float(ra),float(dec)]
		dist=distance(center,pnt)
		if dist < radius[chunk_number]:
			list_to_mosaic.append((dirs[i],obsid[i],dist))
	
	#sort with criterion of minimum distance from center		
	newlist=sorted(list_to_mosaic, key=lambda cr: cr[2])
	'''
	for j in range(len(band)):
		newlist=np.genfromtxt(wd+'murray_sens/mosaic_chunk_'+str(chunk_number)+'_obsids.dat',unpack=True,dtype='str')
		print('Now mosaic '+str(len(newlist))+' obsids for chunk '+str(chunk_number)+'.')
		'''
		w0=open(wd+'murray_sens/mosaic_chunk_'+str(chunk_number)+'_obsids.dat','w')
		#start from closest to center in a 16000x16000 square
		#w=open(wd+'commands.xco','w')
		w2=open(wd+'commands_bkg.xco','w')
		w3=open(wd+'commands_dat.xco','w')
		for i in range(0,len(newlist)):
			if len(newlist[i][1]) == 4:
				stem='0'+newlist[i][1]
			elif len(newlist[i][1]) == 3:
				stem='00'+newlist[i][1]
			elif len(newlist[i][1]) == 5:
				stem=newlist[i][1]
			w0.write(newlist[i][1]+'\n')
		'''
		
		#start from closest to center in a 16000x16000 square
		#w=open(wd+'commands.xco','w')
		w2=open(wd+'commands_bkg.xco','w')
		w3=open(wd+'commands_dat.xco','w')
		for i in range(0,len(newlist)):
			if len(newlist[i]) == 4:
				stem='0'+newlist[i]
			elif len(newlist[i]) == 3:
				stem='00'+newlist[i]
			elif len(newlist[i]) == 5:
				stem=newlist[i]
			if i==0:
				#w.write('read/fits/size=16000 '+wd+'sim_full/acisf'+stem+'_sim_poiss.fits\n')
				w2.write('read/fits/size=4000/rebin=4 '+wd+'/data/'+newlist[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2[j]+'_bkgmap.fits\n')
				#w2.write('read/fits/size=16000 '+wd+'/data/'+newlist[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2[j]+'_expomap.fits\n')
				w3.write('read/fits/size=4000 '+wd+'/data/'+newlist[i]+'/repro_new_asol/acisf'+stem+'_repro_'+band[j]+'keV_4rebinned.img\n')
				#w3.write('read/fits/size=16000 '+wd+'/psfmaps/'+stem+'_r90-x-expo.fits\n')
				
				#w.write('save_image\n')
				w2.write('save_image\n')
				w3.write('save_image\n')
			else:
				#w.write('read/fits/size=3000 '+wd+'sim_full/acisf'+stem+'_sim_poiss.fits\n')
				w2.write('read/fits/rebin=4 '+wd+'/data/'+newlist[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2[j]+'_bkgmap.fits\n')
				#w2.write('read/fits/size=3000 '+wd+'/data/'+newlist[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2[j]+'_expomap.fits\n')
				w3.write('read/fits '+wd+'/data/'+newlist[i]+'/repro_new_asol/acisf'+stem+'_repro_'+band[j]+'keV_4rebinned.img\n')
				#w3.write('read/fits/size=3000 '+wd+'/psfmaps/'+stem+'_r90-x-expo.fits\n')
				
				#w.write('sum_image\n')
				w2.write('sum_image\n')
				w3.write('sum_image\n')
				#w.write('save_image\n')
				w2.write('save_image\n')
				w3.write('save_image\n')
		#w.write('write/fits '+wd+'mosaic_broad_chunk_'+str(chunk_number)+'_normal.fits\n')
		w2.write('write/fits '+wd+'murray_sens/mosaic_'+band2[j]+'_chunk_'+str(chunk_number)+'_bkg_4rebinned.fits\n')
		#w2.write('write/fits '+wd+'mosaic_'+band2[j]+'_chunk_'+str(chunk_number)+'_exp.fits\n')
		#w3.write('write/fits '+wd+'mosaic_'+band2[j]+'_chunk_'+str(chunk_number)+'_r90.fits\n')
		w3.write('write/fits '+wd+'murray_sens/mosaic_'+band2[j]+'_chunk_'+str(chunk_number)+'_dat_4rebinned.fits\n')
		#w.write('exit\n')
		w2.write('exit\n')
		w3.write('exit\n')
		
		#w.close()
		w2.close()
		w3.close()
		#w0.close()

		start_time=time.time()
		if os.path.isfile(wd+'murray_sens/mosaic_'+band2[j]+'_chunk_'+str(chunk_number)+'_normal.fits'):
			pass
		else:
			pass
			#s.call('ximage @'+wd+'commands.xco',shell=True)
		if os.path.isfile(wd+'murray_sens/mosaic_'+band2[j]+'_chunk_'+str(chunk_number)+'_bkg_4rebinned.fits'):
			pass
		else:
			s.call('ximage @'+wd+'commands_bkg.xco',shell=True)
			pass
		if os.path.isfile(wd+'murray_sens/mosaic_'+band2[j]+'_chunk_'+str(chunk_number)+'_dat_4rebinned.fits'):
			pass
		else:
			pass
			s.call('ximage @'+wd+'commands_dat.xco',shell=True)
	
		elapsed_time=time.time()-start_time
		print(float(elapsed_time)/3600.,' hours for the mosaics of chunk '+str(chunk_number)+'.')
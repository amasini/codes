#Use CIAO tool reproject_obs to reprject set of event file sto a common tangent point and create a mosaic.
import numpy as np
import subprocess as s
import sys

def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*xx)**2 +((pointa[1]-pointb[1]))**2)

wd='/Users/alberto/Desktop/XBOOTES/'

ra_c=[218.9132491,217.1703769,218.6973158,217.0839299,218.0404204,216.7921841]
dec_c=[35.20442382,35.25246411,33.93039948,33.92450157,32.83570264,32.82283117]
radius=[0.85,0.85,0.85,0.85,0.6,0.6]

for j in range(1,6):
	list_to_mosaic=[]
	center=[ra_c[j],dec_c[j]]
	(dirs,obsid)=np.genfromtxt(wd+'chunks_of_mosaics_fullres/mosaic_chunk_'+str(j)+'_obsids.dat',unpack=True,dtype='str')
	#w=open(wd+'stack.lis','w')
	for i in range(len(obsid)):
		if len(obsid[i]) == 4:
			stem='0'+obsid[i]
		elif len(obsid[i]) == 3:
			stem='00'+obsid[i]
		elif len(obsid[i]) == 5:
			stem=obsid[i]
		
		#filename=dirs[i]+'/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV.fits'
		#w.write(filename+'\n')
		
		ra=s.check_output('dmkeypar '+wd+dirs[i]+'/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_evt2.fits RA_PNT echo=yes',shell=True)
		dec=s.check_output('dmkeypar '+wd+dirs[i]+'/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_evt2.fits DEC_PNT echo=yes',shell=True)
		pnt=[float(ra),float(dec)]
		dist=distance(center,pnt)
		if dist < radius[j]:
			list_to_mosaic.append((dirs[i],obsid[i],dist))
		'''
		if i==0:
			filename=dirs[i]+'/'+obsid[i]+'/repro/acisf'+stem+'_repro_05to7keV.fits'
			#filename=dirs[i]+'/'+obsid[i]
		else:
			filename=filename+','+dirs[i]+'/'+obsid[i]+'/repro/acisf'+stem+'_repro_05to7keV.fits'
			#filename=filename+','+dirs[i]+'/'+obsid[i]
		'''
	#w.close()
	#s.call('merge_obs infiles='+filename+ outroot=test_merge/ band=broad binsize=1 parallel=yes nproc=4 clobber=yes',shell=True)
	#s.call('merge_obs infiles=@stack.lis[ccd_id=0:3] outroot=test_merge/'+str(j)+' band=broad binsize=1 parallel=yes nproc=4 clobber=yes',shell=True)
	#s.call('reproject_obs infiles=@stack.lis[ccd_id=0:3] outroot=test_merge_back/'+str(j)+' parallel=yes nproc=4 clobber=yes',shell=True)
	
	
	#sort with criterion of minimum distance from center		
	newlist=sorted(list_to_mosaic, key=lambda cr: cr[2])
	print('Now mosaic '+str(len(newlist))+' obsids for chunk '+str(j)+'.')
	#start from closest to center in a 16000x16000 square
	w2=open(wd+'commands_dat.xco','w')
	for i in range(0,len(newlist)):
		if len(newlist[i][1]) == 4:
			stem='0'+newlist[i][1]
		elif len(newlist[i][1]) == 3:
			stem='00'+newlist[i][1]
		elif len(newlist[i][1]) == 5:
			stem=newlist[i][1]
		#creates full band evt file binned to native pixel scale
		s.call('rm -f '+wd+'data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV.img',shell=True)
		#s.call('dmcopy \"'+wd+newlist[i][0]+'/'+newlist[i][1]+'/repro_new_asol/acisf'+stem+'_repro_evt2.fits[events][ccd_id=0:3][energy=500:7000]\" \"'+wd+newlist[i][0]+'/'+newlist[i][1]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV.img\" opt=image clobber=yes',shell=True)
		if i==0:
			w2.write('read/fits/size=16000 '+wd+newlist[i][0]+'/'+newlist[i][1]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV.img\n')
			w2.write('save_image\n')
		else:
			w2.write('read/fits/size=3000 '+wd+newlist[i][0]+'/'+newlist[i][1]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV.img\n')
			w2.write('sum_image\n')
			w2.write('save_image\n')
	w2.write('write/fits '+wd+'mosaic_broad_chunk_'+str(j)+'_dat.fits\n')
	w2.write('exit\n')
	w2.close()

	s.call('ximage @'+wd+'commands_dat.xco',shell=True)
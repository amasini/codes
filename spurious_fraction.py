# THIS SCRIPT MAKES A POISSONIAN REALIZATION OF THE BKGMAPS, IN ORDER TO DETECT SPURIOUS 
# SOURCES - THE PROBLEM IS THAT I WAS FINDING NO SOURCES WITH WAVDETECT - TRY AGAIN WITH
# 4X4 MOCK IMAGES?
#make Poissonian realization of background maps
import numpy as np
import subprocess as s
import sys
from astropy.io import fits

wd='/Users/alberto/Desktop/XBOOTES/'

obsid=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=1,dtype='str')

for band in ['broad']:
    #create Poissonian realization of backmaps (they are instrumental only...)    
    for obs in obsid:
        if len(obs) == 4:
        	stem='0'+obs
        elif len(obs) == 3:
        	stem='00'+obs
        elif len(obs) == 5:
        	stem=obs
        
        print(obs)
        #take bkgmap
        bkgmap=fits.open(wd+'data/'+obs+'/repro_new_asol/out/acisf'+stem+'_'+band+'_bkgmap.fits')
        back=bkgmap[0].data
        backheader=bkgmap[0].header

        #make Poissonian and save file
        poiss_back=np.random.poisson(back)
        hdu1 = fits.PrimaryHDU(poiss_back,header=backheader)
        hdu1.writeto(wd+'data/'+obs+'/repro_new_asol/out/acisf'+stem+'_'+band+'_bkgmap_poiss.fits',overwrite=True)

    #create full resolution back mosaics
    ra_c=[218.9132491,217.1703769,218.6973158,217.0839299,218.0404204,216.7921841]
    dec_c=[35.20442382,35.25246411,33.93039948,33.92450157,32.83570264,32.82283117]
    radius=[0.85,0.85,0.85,0.85,0.6,0.6]

    for chunk_number in range(0,len(ra_c)):

        newlist=np.genfromtxt(wd+'mosaic_chunk_'+str(chunk_number)+'_obsids.dat',unpack=True,dtype='str')
        print('Now mosaic '+str(len(newlist))+' obsids for chunk '+str(chunk_number)+'.')
        w3=open(wd+'commands_dat.xco','w')
        for i in range(0,len(newlist)):
            if len(newlist[i]) == 4:
                stem='0'+newlist[i]
            elif len(newlist[i]) == 3:
		        	stem='00'+newlist[i]
            elif len(newlist[i]) == 5:
			    stem=newlist[i]
            if i==0:
			    w3.write('read/fits/size=16000 '+wd+'data/'+newlist[i]+'/repro_new_asol/out/acisf'+stem+'_'+band+'_bkgmap.fits\n')
			    w3.write('save_image\n')
            else:
			    w3.write('read/fits/size=3000 '+wd+'data/'+newlist[i]+'/repro_new_asol/out/acisf'+stem+'_'+band+'_bkgmap.fits\n')
			    w3.write('sum_image\n')
			    w3.write('save_image\n')
        w3.write('write/fits '+wd+'chunks_of_mosaics_fullres/mosaic_'+band+'_chunk_'+str(chunk_number)+'_bkg_poiss.fits\n')
        w3.write('exit\n')
		
        w3.close()

        s.call('ximage @'+wd+'commands_dat.xco',shell=True)


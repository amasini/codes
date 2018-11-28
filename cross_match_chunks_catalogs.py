# THIS CODE CROSS-MATCHES CATALOGS OF CHUNKS OF MOSAIC, MERGING THEM IN COUPLES AND
# USING FTMERGE TO APPEND ONE CATALOG TO THE OTHER
import numpy as np
import sys
from astropy.io import fits
import subprocess as s

def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*xx*3600)**2 +((pointa[1]-pointb[1])*3600)**2)

wd='/Users/alberto/Desktop/XBOOTES/murray_sens/'
for band in ['broad']:
	for k in [0,2,4]:

		#w=open(wd+str(k)+'-'+str(k+1)+'.dat','w')

		#open first catalog
		cat0=fits.open(wd+'mosaic_'+band+'_'+str(k)+'_src_4reb.fits')
		#cat0=fits.open(wd+'mosaic_'+band+'_'+str(k)+'_src_new.fits')
		dat0=cat0[1].data
		ra0=cat0[1].data['RA']
		dec0=cat0[1].data['DEC']
		cat0.close()
	
		#open second catalog
		cat1=fits.open(wd+'mosaic_'+band+'_'+str(k+1)+'_src_4reb.fits')
		#cat1=fits.open(wd+'mosaic_'+band+'_'+str(k+1)+'_src_new.fits')
		dat1=cat1[1].data
		ra1=cat1[1].data['RA']
		dec1=cat1[1].data['DEC']
		cat1.close()
	
		goodvalues_cat0=[]
		for i in range(len(ra0)):
			p0=[ra0[i],dec0[i]]
			double=0
		
			for j in range(len(ra1)):
				p1=[ra1[j],dec1[j]]
				#write to the output file all the sources of cat1
				if i==0:
					pass
					#w.write(str(ra1[j])+'\t'+str(dec1[j])+'\n')
				if distance(p0,p1) < 0.5:
					double=1
		
			#add to the output files only the unique sources of cat0
			if double==0:
				goodvalues_cat0.append(i)
				#w.write(str(ra0[i])+'\t'+str(dec0[i])+'\n')
	
		newdata=dat0[goodvalues_cat0]
		hdu = fits.BinTableHDU(data=newdata)
		
		hdu.writeto(wd+'unique_'+str(k)+'_dat_src.fits',overwrite=True)
		#hdu.writeto(wd+'unique_'+str(k)+'_src.fits',overwrite=True)
	
		s.call('ftmerge '+wd+'unique_'+str(k)+'_dat_src.fits,'+wd+'mosaic_'+band+'_'+str(k+1)+'_src_4reb.fits '+wd+'out_'+str(k)+'-'+str(k+1)+'_dat.fits clobber=yes',shell=True)
		#s.call('ftmerge '+wd+'unique_'+str(k)+'_src.fits,'+wd+'mosaic_'+str(k+1)+'_src_new_bothscales.fits out_'+str(k)+'-'+str(k+1)+'.fits clobber=yes',shell=True)
		#w.close()


	#open first catalog
	cat0=fits.open(wd+'out_0-1_dat.fits')
	#cat0=fits.open(wd+'out_0-1.fits')
	dat0=cat0[1].data
	ra0=cat0[1].data['RA']
	dec0=cat0[1].data['DEC']
	cat0.close()
	
	#open second catalog
	cat1=fits.open(wd+'out_2-3_dat.fits')
	#cat1=fits.open(wd+'out_2-3.fits')
	dat1=cat1[1].data
	ra1=cat1[1].data['RA']
	dec1=cat1[1].data['DEC']
	cat1.close()
		
	goodvalues_cat0=[]
	for i in range(len(ra0)):
		p0=[ra0[i],dec0[i]]
		double=0	
		
		for j in range(len(ra1)):
			p1=[ra1[j],dec1[j]]
			#write to the output file all the sources of cat1
			if distance(p0,p1) < 0.5:
				double=1
		
		#add to the output files only the unique sources of cat0
		if double==0:
			goodvalues_cat0.append(i)
	
	newdata=dat0[goodvalues_cat0]

	hdu = fits.BinTableHDU(data=newdata)

	hdu.writeto(wd+'unique_0-1_dat_src.fits',overwrite=True)
	#hdu.writeto(wd+'unique_0-1_src.fits',overwrite=True)

	s.call('ftmerge '+wd+'unique_0-1_dat_src.fits,'+wd+'out_2-3_dat.fits '+wd+'out_0-1-2-3_dat.fits clobber=yes',shell=True)
	#s.call('ftmerge '+wd+'unique_0-1_src.fits,'+wd+'out_2-3.fits out_0-1-2-3.fits clobber=yes',shell=True)


	#open first catalog
	cat0=fits.open(wd+'out_0-1-2-3_dat.fits')	
	#cat0=fits.open(wd+'out_0-1-2-3.fits')
	dat0=cat0[1].data
	ra0=cat0[1].data['RA']
	dec0=cat0[1].data['DEC']
	cat0.close()
	
	#open second catalog
	cat1=fits.open(wd+'out_4-5_dat.fits')
	#cat1=fits.open(wd+'out_4-5.fits')
	dat1=cat1[1].data
	ra1=cat1[1].data['RA']
	dec1=cat1[1].data['DEC']
	cat1.close()
	
	goodvalues_cat0=[]
	for i in range(len(ra0)):
		p0=[ra0[i],dec0[i]]
		double=0	
		
		for j in range(len(ra1)):
			p1=[ra1[j],dec1[j]]
			#write to the output file all the sources of cat1
			if distance(p0,p1) < 0.5:
				double=1
		
		#add to the output files only the unique sources of cat0
		if double==0:
			goodvalues_cat0.append(i)
	
	newdata=dat0[goodvalues_cat0]

	hdu = fits.BinTableHDU(data=newdata)

	hdu.writeto(wd+'unique_0-1-2-3_dat_src.fits',overwrite=True)
	#hdu.writeto(wd+'unique_0-1-2-3_src.fits',overwrite=True)

	s.call('ftmerge '+wd+'unique_0-1-2-3_dat_src.fits,'+wd+'out_4-5_dat.fits '+wd+'out_'+band+'_0-1-2-3-4-5_dat_4reb.fits clobber=yes',shell=True)
	#s.call('ftmerge '+wd+'unique_0-1-2-3_src.fits,'+wd+'out_4-5.fits out_0-1-2-3-4-5_new_bothscales.fits clobber=yes',shell=True)


sys.exit()

'''
w=open(wd+'0-1-2-3.dat','w')

#open first catalog
(ra0,dec0)=np.genfromtxt(wd+'0-1.dat',unpack=True)

for i in range(len(ra0)):
	p0=[ra0[i],dec0[i]]
	double=0
	
	#open second catalog
	(ra1,dec1)=np.genfromtxt(wd+'2-3.dat',unpack=True)
	for j in range(len(ra1)):
		p1=[ra1[j],dec1[j]]
		#write to the output file all the sources of cat1
		if i==0:
			w.write(str(ra1[j])+'\t'+str(dec1[j])+'\n')
		if distance(p0,p1) < 0.5:
			double=1
		
	#add to the output files only the unique sources of cat0
	if double==0:
		w.write(str(ra0[i])+'\t'+str(dec0[i])+'\n')

w.close()


w=open(wd+'0-1-2-3-4-5.dat','w')

#open first catalog
(ra0,dec0)=np.genfromtxt(wd+'0-1-2-3.dat',unpack=True)

for i in range(len(ra0)):
	p0=[ra0[i],dec0[i]]
	double=0
	
	#open second catalog
	(ra1,dec1)=np.genfromtxt(wd+'4-5.dat',unpack=True)
	for j in range(len(ra1)):
		p1=[ra1[j],dec1[j]]
		#write to the output file all the sources of cat1
		if i==0:
			w.write(str(ra1[j])+'\t'+str(dec1[j])+'\n')
		if distance(p0,p1) < 0.5:
			double=1
		
	#add to the output files only the unique sources of cat0
	if double==0:
		w.write(str(ra0[i])+'\t'+str(dec0[i])+'\n')

w.close()
'''
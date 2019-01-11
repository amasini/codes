#cross-correlate Xbootes and CDWFS catalogs
import numpy as np
import sys
from astropy.io import fits
import matplotlib.pyplot as plt

def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*xx*3600)**2 +((pointa[1]-pointb[1])*3600)**2)

wd='/Users/alberto/Desktop/XBOOTES/'

#take first catalog
cat1=fits.open(wd+'xbootes_kenter_cat0.fits')
#cat1=fits.open(wd+'murray_sens/xbootes_broad_cat0_3.fits')
ra_k=cat1[1].data['RA']
dec_k=cat1[1].data['DEC']
cts_full=cat1[1].data['TOT']


#w=open(wd+'murray_sens/xbootes_kenter_cat0.reg','w')
#for i in range(len(ra_k)):
#	w.write('circle('+str(ra_k[i])+'d, '+str(dec_k[i])+'d, 5\") #width=2 color=red \n')
#w.close()
print(len(ra_k))

#take second catalog
#cat2=fits.open(wd+'cdwfs_broad_cat0.fits')
#cat2=fits.open(wd+'xbootes_kenter_cat0.fits')
cat2=fits.open(wd+'murray_sens/mosaic_broad_cat0_3.fits')
#cat2=fits.open(wd+'murray_sens/out_src_fullmosaic4x4.fits')
#cat2=fits.open(wd+'chunks_of_mosaics_fullres/out_broad_0-1-2-3-4-5_dat.fits')
ra_cdwfs=cat2[1].data['RA']
dec_cdwfs=cat2[1].data['DEC']
cts_full_cdwfs=cat2[1].data['TOT']

#ra_cdwfs=ra_cdwfs[cts_full_cdwfs>=4.]
#dec_cdwfs=dec_cdwfs[cts_full_cdwfs>=4.]
#w=open(wd+'murray_sens/xbootes_broad_cat0_3.reg','w')
#for i in range(len(ra_cdwfs)):
#	w.write('circle('+str(ra_cdwfs[i])+'d, '+str(dec_cdwfs[i])+'d, 20\") #width=2 color=cyan \n')
#w.close()
print(len(ra_cdwfs))


match_rad=10.0
#w=open(wd+'murray_sens/xbootes_unmatched_'+str(int(match_rad))+'arcsec_3.reg','w')
unmatched=0
cts=[]
for i in range(len(ra_k)):
	kenter_source=[ra_k[i],dec_k[i]]
	found=0
	for j in range(len(ra_cdwfs)):
		cdwfs_source=[ra_cdwfs[j],dec_cdwfs[j]]
		d=distance(kenter_source,cdwfs_source)
		if d < match_rad:
			found=1
	
	if found==0:
		unmatched=unmatched+1
		#w.write('circle('+str(ra_k[i])+'d, '+str(dec_k[i])+'d, 20\") #width=2 color=cyan \n')
		cts.append(cts_full[i])
		if cts_full[i] > 10.0:
			print(ra_k[i],dec_k[i])

#w.close()

print(unmatched)
plt.figure()
plt.hist(cts,bins=10)
plt.show()
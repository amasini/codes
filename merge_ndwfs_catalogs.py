# Code to merge the NDWFS strips catalogs - IDL is WAY faster
import numpy as np
from astropy.io import fits
from astropy.table import Table
import sys
import time

t_start=time.time()

def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*3600*xx)**2 +((pointa[1]-pointb[1])*3600)**2)

wd='/Users/alberto/Downloads/bootes_mw_data/'
'''
# Open full band catalog (no multiple, output from clean_multiple_sources.py)
#fcat=fits.open(wd+'NDWFS_I_32_33_cat_m-cp.fits')
#filef=fcat[1].data
#cols=fcat[1].columns
#print(cols.names)
#fcat.close()
#sys.exit()

# Open soft band catalog (no multiple)
#scat=fits.open(wd+'NDWFS_I_33_34_cat_m-cp.fits')
#files=scat[1].data

outcat=[] # List containing the output file
#outcat = np.empty((2066890,101))
#j=0
for i in range(32,36):
	fcat=fits.open(wd+'NDWFS_I_'+str(i)+'_'+str(i+1)+'_cat_m-cp.fits')
	filef=fcat[1].data
	cols=fcat[1].columns
	fcat.close()
	
	names=cols.names
	for linef in range(len(filef)):
		outcat.append(filef[linef])
		#j=j+1

print(len(names),len(outcat))

outcat=np.array(outcat)

#print(len(outcat[:,0]))
print((time.time()-t_start)/60.,'minutes')

vet_col=[]
for i in range(len(names)):
	vet_col.append(outcat[:,i])

cat=Table(vet_col,names=names)
cat.write(wd+'NDWFS_I_merged.fits',format='fits',overwrite=True)

sys.exit()
'''
# Create a FITS file SDWFS (Ashby+09) catalog
'''
names=np.genfromtxt(wd+'SDWFS_ch1_stack.v34.txt',skip_header=17,skip_footer=677525,dtype='str',delimiter='|',autostrip=True)
names=names[1:-1]

a=np.genfromtxt(wd+'SDWFS_ch1_stack.v34.txt',unpack=True,skip_header=21)
vet_col=[]
for i in range(len(names)):
	vet_col.append(a[i])

cat=Table(vet_col,names=names)
cat.write(wd+'SDWFS_3.6_sel.fits',format='fits',overwrite=True)
'''
# Create a FITS file IBIS Ks (Gonzalez+10) catalog
'''
names=np.genfromtxt(wd+'bootes_dr1/BOOTES_ks_V1.0-cp.cat',unpack=True,usecols=1,skip_footer=157613,dtype='str')

a=np.genfromtxt(wd+'bootes_dr1/BOOTES_ks_V1.0-cp.cat',unpack=True,skip_header=109)
vet_col=[]
for i in range(len(names)):
	vet_col.append(a[i])

cat=Table(vet_col,names=names)
cat.write(wd+'IBIS_Ks.fits',format='fits',overwrite=True)
'''
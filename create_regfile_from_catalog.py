# Simple script to create region file from catalog

import numpy as np
from astropy.io import fits

wd='/Users/alberto/Desktop/XBOOTES/'
'''
# Define band
band='soft'

if band=='broad':
	color='green'
	cut=6e-5
elif band=='soft':
	color='red'
	cut=1.4e-4
else:
	color='cyan'
	cut=6e-5

# Open the catalog
cat=fits.open(wd+'cdwfs_'+band+'_cat1.fits')

ra=cat[1].data['RA']
dec=cat[1].data['DEC']
r90=cat[1].data['AV_R90']
prob=cat[1].data['PROB']

ra=ra[prob<=cut]
dec=dec[prob<=cut]
r90=r90[prob<=cut]

w=open(wd+'cdwfs_'+band+'_cat2.reg','w')
for i in range(len(ra)):
	w.write('circle('+str(ra[i])+'d,'+str(dec[i])+'d,'+str(r90[i])+'\") # width=2 color='+color+'\n')
w.close()

##########################
# Open the catalog
cat=fits.open(wd+'cdwfs_merged_cat0.fits')

cutf,cuts,cuth=6e-5,1.4e-4,6e-5

raf=cat[1].data['RA_F']
decf=cat[1].data['DEC_F']
r90f=cat[1].data['R90_F']
probf=cat[1].data['PROB_F']

ras=cat[1].data['RA_S']
decs=cat[1].data['DEC_S']
r90s=cat[1].data['R90_S']
probs=cat[1].data['PROB_S']

rah=cat[1].data['RA_H']
dech=cat[1].data['DEC_H']
r90h=cat[1].data['R90_H']
probh=cat[1].data['PROB_H']
'''

# Open the catalog
cat=fits.open('/Users/alberto/Desktop/XBOOTES/new_mosaics_detection/cdwfs_merged_cat1.fits')

ra=cat[1].data['RA']
dec=cat[1].data['DEC']

color='yellow'
w=open(wd+'new_mosaics_detection/cdwfs_merged_cat1.reg','w')
for i in range(len(ra)):
	w.write('circle('+str(ra[i])+'d,'+str(dec[i])+'d,4.5\") # width=2 color='+color+'\n')
w.close()
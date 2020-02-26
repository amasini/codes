# Simple script to create region file from catalog

import numpy as np
from astropy.io import fits
import sys 

wd='/Users/alberto/Downloads/'
input = 'bootes_LOFAR_lba_retana-montenegro18.fit'
output = 'bootes_LOFAR_lba_retana-montenegro18.reg'

size = '8'
color='yellow'

# Open the catalog
cat=fits.open(wd+input)
ra=cat[1].data['RAJ2000']
dec=cat[1].data['DEJ2000']
cat.close()

w=open(wd+output,'w')
for i in range(len(ra)):
	w.write('circle('+str(ra[i])+'d,'+str(dec[i])+'d,'+size+'\") # width=2 color='+color+'\n')
w.close()
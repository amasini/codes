#THIS SCRIPT IS SIMILAR TO PLOT.PY - QUICK LOOK TO DATA AND SIMULATIONS, IN PARTICULAR
# TO THEIR NUMBER OS SOURCES AT A GIVEN NUMBER OF COUNTS
import sys as s
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

wd='/Users/alberto/Desktop/XBOOTES/chunks_of_mosaics_fullres/'

dat=fits.open(wd+'out_0-1-2-3-4-5_new_bothscales.fits')
cts_sim=dat[1].data['NET_COUNTS']

dat2=fits.open(wd+'out_0-1-2-3-4-5_dat_bothscales.fits')
cts_dat=dat2[1].data['NET_COUNTS']

mincts=7
maxcts=14000
cts_dat2=cts_dat[cts_dat<maxcts]
cts_sim2=cts_sim[cts_sim<maxcts]

bins=np.logspace(np.log10(mincts),np.log10(maxcts),20)
bincenters=list((bins[i+1]+bins[i])/2 for i in range(len(bins)-1))

a,b=np.histogram(cts_dat2,bins)
c,d=np.histogram(cts_sim2,bins)

cum_dat=np.cumsum(a)
cum_sim=np.cumsum(c)

print(max(cum_dat),max(cum_sim))

plt.figure()
plt.plot(bincenters,cum_dat,'k-',label='Data')
plt.plot(bincenters,cum_sim,'b-',label='Sim')
plt.xlabel('Net counts')
plt.ylabel('# of sources')
plt.legend()
plt.tight_layout()
plt.show()

plt.figure()
plt.hist(cts_dat2,bins=bins,histtype='step',color='black',label='Data')
plt.hist(cts_sim2,bins=bins,histtype='step',color='blue',label='Sim')
plt.xscale('log')
plt.xlabel('Net counts')
plt.legend()
plt.tight_layout()
plt.show()

diff=a-c

plt.figure()
plt.plot(bincenters,diff,'k-')
plt.show()
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys

wd='/Users/alberto/Desktop/XBOOTES/'

dat=fits.open(wd+'mosaic_broad_expomap_4rebinned.fits')

ima=dat[0].data
dat.close()
ima=ima[~np.isnan(ima)]
ima=ima[ima>0]
ima=ima/(16.0*1e3)

bins=np.logspace(np.log10(min(ima)),np.log10(max(ima)),200)
a,b=np.histogram(ima,bins=bins)

cum=np.cumsum(-a)+np.max(np.cumsum(a))

cha_pixscale=0.492
pixscale=4*cha_pixscale
pixarea_deg=(pixscale/3600.)**2
cum=cum*pixarea_deg

bincenters=[]
for i in range(len(b)-1):
	bincenters.append((b[i]+b[i+1])/2.)

xticks=[1,10,100]
xlabels=['1','10','100']
yticks=[0.1,1,10]
ylabels=['0.1','1','10']

f,ax=plt.subplots(1)
ax.plot(bincenters,cum,'b-')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Exposure (ks)',fontsize=14)
ax.set_ylabel(r'Area (deg$^2$)',fontsize=14)
ax.tick_params(axis='both',which='major',direction='in',labelsize=14,length=8)
ax.tick_params(axis='both',which='minor',direction='in',labelsize=14,length=6)
ax.set_xticks(xticks)
ax.set_xticklabels(xlabels)
ax.set_yticks(yticks)
ax.set_yticklabels(ylabels)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.axis([1,120,0.1,10])
plt.tight_layout()
plt.savefig(wd+'area_expo.png',format='png')
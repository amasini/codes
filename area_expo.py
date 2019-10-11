import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys

wd='/Users/alberto/Desktop/XBOOTES/'
band='broad'

dat=fits.open(wd+'new_mosaics_detection/cdwfs_'+band+'_expomap_4reb.fits')
ima=dat[0].data
dat.close()

dat=fits.open(wd+'new_mosaics_detection/xbootes_'+band+'_expomap_4reb.fits')
ima2=dat[0].data
dat.close()

ima=ima[~np.isnan(ima)]
ima=ima[ima>0]
ima=ima/1e3

ima2=ima2[~np.isnan(ima2)]
ima2=ima2[ima2>0]
ima2=ima2/1e3

#bins=np.logspace(np.log10(min(ima)),np.log10(max(ima)),200)
bins=np.linspace(0.1,np.max(ima),401)
a,b=np.histogram(ima,bins=bins)
c,d=np.histogram(ima2,bins=bins)

cum=np.cumsum(-a)+np.max(np.cumsum(a))
cum2=np.cumsum(-c)+np.max(np.cumsum(c))

cha_pixscale=0.492
pixscale=4*cha_pixscale
pixarea_deg=(pixscale/3600.)**2

cum=cum*pixarea_deg
cum2=cum2*pixarea_deg

bincenters=list((b[i]+b[i+1])/2. for i in range(len(b)-1))

xticks=[1,10,100]
xlabels=['1','10','100']
yticks=[0.1,1,10]
ylabels=['0.1','1','10']

f,ax=plt.subplots(1)
ax.plot(bincenters,cum,'k-',linewidth=3,label='CDWFS')
ax.plot(bincenters,cum2,'k--',label='XBootes')

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
plt.legend()
plt.tight_layout()
plt.show()
#plt.savefig(wd+'area_expo.png',format='png')
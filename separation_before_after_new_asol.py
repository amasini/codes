#This script computes separation between common sources from wavdetect and brand+05,
#before and after reprojecting the wcs system.

import numpy as np
import sys
import subprocess as s
import matplotlib.pyplot as plt
from astropy.io import fits
import matplotlib.gridspec as gridspec

def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*3600*xx)**2 +((pointa[1]-pointb[1])*3600)**2)

wd='/Users/alberto/Desktop/XBOOTES/'

obsid=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=1,dtype='str')
delta_ra,delta_dec,dist=[],[],[]
delta_ra2,delta_dec2,dist2=[],[],[]

for i in range(len(obsid)):
	if len(obsid[i]) == 4:
		stem='0'+obsid[i]
	elif len(obsid[i]) == 3:
		stem='00'+obsid[i]
	elif len(obsid[i]) == 5:
		stem=obsid[i]
	
	#get X-ray sources from wavdetect
	xrayfile=fits.open(wd+'data/'+obsid[i]+'/src_to_reproject.fits')
	xray=xrayfile[1].data
	ra_x=xray['RA']
	dec_x=xray['DEC']
	xrayfile.close()
	
	#get new X-ray sources from wavdetect after new asp sol has been applied
	xrayfile2=fits.open(wd+'data/'+obsid[i]+'/src_new_asol.fits')
	xray2=xrayfile2[1].data
	ra_x2=xray2['RA']
	dec_x2=xray2['DEC']
	xrayfile2.close()
	
	#get optical source from brand+05
	optfile=fits.open(wd+'data/'+obsid[i]+'/acisf'+stem+'_brand+05_renamedcols.fits')
	opt=optfile[1].data
	ra_o=opt['RA']
	dec_o=opt['DEC']
	optfile.close()
	
	for j in range(len(ra_o)):
		src_o=[ra_o[j],dec_o[j]]
		for k in range(len(ra_x)):
			src_x=[ra_x[k],dec_x[k]]
			d=distance(src_o,src_x)
			if d < 2.0: #if found a match
				xx = np.cos(src_o[1]/180*3.141592)
				dist.append(d)
				delta_ra.append((src_o[0]-src_x[0])*3600*xx)
				delta_dec.append((src_o[1]-src_x[1])*3600)
		
		for m in range(len(ra_x2)):
			src_x2=[ra_x2[m],dec_x2[m]]
			d2=distance(src_o,src_x2)
			if d2 < 2.0: #if found a match
				xx = np.cos(src_o[1]/180*3.141592)
				dist2.append(d2)
				delta_ra2.append((src_o[0]-src_x2[0])*3600*xx)
				delta_dec2.append((src_o[1]-src_x2[1])*3600)
	
		
print('Percentiles before:',np.percentile(dist,68),np.percentile(dist,90),np.percentile(dist,95))
print('Percentiles after:',np.percentile(dist2,68),np.percentile(dist2,90),np.percentile(dist2,95))
print(np.mean(delta_ra),np.median(delta_ra),np.mean(delta_dec),np.median(delta_dec))
print(np.mean(delta_ra2),np.median(delta_ra2),np.mean(delta_dec2),np.median(delta_dec2))

plt.figure()
plt.hist(dist,bins=20,histtype='step')
plt.hist(dist2,bins=20,histtype='step')
plt.show()

circle1 = plt.Circle((0, 0), np.percentile(dist,68), color='g',linewidth=2,fill=False)
circle2 = plt.Circle((0, 0), np.percentile(dist,90), color='g',linewidth=2,fill=False)
circle3 = plt.Circle((0, 0), np.percentile(dist,95), color='g',linewidth=2,fill=False)

circle4 = plt.Circle((0, 0), np.percentile(dist2,68), color='k',linewidth=2,fill=False)
circle5 = plt.Circle((0, 0), np.percentile(dist2,90), color='k',linewidth=2,fill=False)
circle6 = plt.Circle((0, 0), np.percentile(dist2,95), color='k',linewidth=2,fill=False)

fig = plt.figure(tight_layout=True,figsize=[6,6])
gs = gridspec.GridSpec(3, 3)
gs.update(wspace=0., hspace=0.)

ax1 = plt.subplot(gs[0:1, 0:2])
ax1.hist(delta_ra,bins=20,histtype='step',color='g',density=True,label='Before')
ax1.hist(delta_ra2,bins=20,histtype='step',color='k',density=True,label='After')
ax1.axvline(x=0,color='k')
ax1.xaxis.set_visible(False)
ax1.yaxis.set_visible(False)
ax1.legend()

ax2 = plt.subplot(gs[1:3, 2:3])
ax2.hist(delta_dec,bins=20,histtype='step',color='g',density=True,orientation='horizontal')
ax2.hist(delta_dec2,bins=20,histtype='step',color='k',density=True,orientation='horizontal')
ax2.axhline(y=0,color='k')
ax2.xaxis.set_visible(False)
ax2.yaxis.set_visible(False)

ax3 = plt.subplot(gs[1:3, 0:2])
ax3.plot(delta_ra,delta_dec,'g.',ms=1)
ax3.plot(delta_ra2,delta_dec2,'k.',ms=1)
ax3.set_xlabel(r'$\Delta$ RA (arcsec)')
ax3.set_ylabel(r'$\Delta$ DEC (arcsec)')
ax3.set_yticks(np.arange(-2, 3, 1))
ax3.axhline(y=0,color='k')
ax3.axvline(x=0,color='k')
ax3.add_artist(circle1)
ax3.add_artist(circle2)
ax3.add_artist(circle3)
ax3.add_artist(circle4)
ax3.add_artist(circle5)
ax3.add_artist(circle6)
plt.savefig(wd+'astrometric_correction.pdf',format='pdf')
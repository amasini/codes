# THIS SCRIPT IS USELESS, IS JUST TO QUICKLOOK DATA AND CAN BE CHANGED EVERYTIME IS NEEDED
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from scipy.stats import poisson
import subprocess as s

wd='/Users/alberto/Desktop/XBOOTES/'
'''
obsid=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=1,dtype='str')
for i in range(len(obsid)):
	if len(obsid[i]) == 4:
		stem='0'+obsid[i]
	elif len(obsid[i]) == 3:
		stem='00'+obsid[i]
	elif len(obsid[i]) == 5:
		stem=obsid[i]
	s.call('rm -f '+wd+'data/'+obsid[i]+'/repro_new_asol/out/*soft*',shell=True)
	
print('Done')
sys.exit()
'''
(ra,dec,flux)=np.genfromtxt(wd+'cdwfs_broad_sim_poiss_matched.dat',unpack=True,skip_header=1)
bins=np.logspace(np.log10(5e-17),np.log10(9e-13),100)
a2,b2=np.histogram(flux,bins=bins)
for i in range(len(a2)):
	print(b2[i],a2[i])

ra=ra[flux>5.5e-14]
dec=dec[flux>5.5e-14]
flux=flux[flux>5.5e-14]
print(len(ra))
w=open(wd+'simul_double_matches.reg','w')
for j in range(len(ra)):
	w.write('circle('+str(ra[j])+'d, '+str(dec[j])+'d, 20\") #width=2 color=cyan \n')
	print(ra[j],dec[j],flux[j])
w.close()
sys.exit()
'''
N=np.genfromtxt(wd+'wav_xbootes_singleacis_5e-5/xbootes_Nsources_3.dat',unpack=True,usecols=0)

r = poisson.rvs(35.07, size=126)
#r = poisson.rvs(27.22, size=126)
#N=N[N<=60]
print(np.mean(N))
print(len(N[N>60]))
plt.figure()
bins=(np.max(N)-np.min(N))/4
print(np.max(N),np.min(N))
plt.hist(N,bins=int(bins),histtype='step',alpha=0.5)
plt.hist(r,bins=int(bins),histtype='step',alpha=0.5)
plt.show()
sys.exit()
'''

#cat=fits.open(wd+'cdwfs_broad_cat0.fits')
#cat=fits.open(wd+'murray_sens/xbootes_unique_2.fits')
#cat=fits.open(wd+'murray_sens/xbootes_broad_cat0.fits')
cat=fits.open(wd+'mosaic_broad_sim_poiss_cat0_3.fits') #7881 SIMULATIONS [#7211 sources] BLUE
ra=cat[1].data['RA']
dec=cat[1].data['DEC']
detml=cat[1].data['DETML']
net=cat[1].data['NET']

catb=fits.open(wd+'mosaic_broad_cat0_3.fits') #8800 sources RED
rab=catb[1].data['RA']
decb=catb[1].data['DEC']
detmlb=catb[1].data['DETML']
netb=catb[1].data['NET']
#$r=cat[1].data['AV_R90']
#cts=cat[1].data['NET']
#flux=cat[1].data['FLUX']

#inpflux=cat[1].data['CTRP_FLUX']
prob=np.e**(-detml)
probb=np.e**(-detmlb)
print(len(ra[net>200.]))
net2=net[net<200]
prob2=prob[net<200]

net3=netb[netb<200]
prob3=probb[netb<200]

plt.figure()
plt.hist(net2,bins=200,color='blue',alpha=0.5)
plt.hist(net3,bins=200,color='red',alpha=0.5)
plt.show()

bins=np.logspace(np.log10(1e-10),np.log10(1),11)
plt.figure()
plt.hist(prob2,bins=bins,color='blue',alpha=0.5)
plt.hist(prob3,bins=bins,color='red',alpha=0.5)
plt.xscale('log')
plt.show()

plt.figure()
plt.plot(prob2,net2,'b.')
plt.plot(prob3,net3,'r.')
plt.xscale('log')
#plt.axis([1e-100,1,-10,400])
plt.show()

#print(len(ra[net>2.]))
#print(len(ra[net>4.]))
print(len(ra[prob>=1e-4]))
print(len(rab[probb>=1e-4]))
sys.exit()

#w=open(wd+'xbootes_ugly_P>1e-4.reg','w')
#for i in range(len(ra2)):
#	w.write('circle('+str(ra2[i])+'d, '+str(dec2[i])+'d, 4\") #color=green\n')
#w.close()


#print(len(ra[net>=4.]))



wd='/Users/alberto/Desktop/XBOOTES/'
band='broad'
band2='05to7'
if band=='broad':
	cf=1.825E-11
elif band=='soft':
	cf=1.632E-11
elif band=='hard':
	cf=2.016E-11
	
#take catalog of detected sources
dat=fits.open(wd+'chunks_of_mosaics_fullres/out_'+band+'_0-1-2-3-4-5_dat.fits')
src_ra=dat[1].data['RA']
src_dec=dat[1].data['DEC']
src_sign=dat[1].data['SRC_SIGNIFICANCE']
src_net=dat[1].data['NET_COUNTS']
dat.close()

#ignore detected sources with NET_CTS<cut
cut=1.0
src_ra=src_ra[src_net>=cut]
src_dec=src_dec[src_net>=cut]
src_sign=src_sign[src_net>=cut]

w=open(wd+'cdwfs_broad_output_wavdetect.reg','w')
for i in range(len(src_ra)):
	w.write('circle('+str(src_ra[i])+'d, '+str(src_dec[i])+'d, 5\") #color=yellow \n')
	#w.write('circle('+str(ra[i])+'d, '+str(dec[i])+'d, 3\") #color=cyan \n')
w.close()
'''

w=open(wd+'cdwfs_broad_cat0.reg','w')
for i in range(len(ra)):
	w.write('circle('+str(ra[i])+'d, '+str(dec[i])+'d, '+str(r[i])+'\") #color=cyan \n')
	#w.write('circle('+str(ra[i])+'d, '+str(dec[i])+'d, 3\") #color=cyan \n')
w.close()
sys.exit()

(ra_h,ra_m,ra_s,dec_d,dec_m,dec_s)=np.genfromtxt(wd+'xbootes_kenter.txt',unpack=True,skip_header=45,usecols=[2,3,4,5,6,7],dtype='str')
(cts_f,cts_s,cts_h,flux_f,flux_s,flux_h)=np.genfromtxt(wd+'xbootes_kenter.txt',unpack=True,skip_header=45,usecols=[9,10,11,12,15,18])

out_ra,out_dec=[],[]
w=open(wd+'xbootes_kenter.reg','w')
for i in range(len(ra_h)):
	c = SkyCoord(ra_h[i]+' '+ra_m[i]+' '+ra_s[i]+' +'+dec_d[i]+' '+dec_m[i]+' '+dec_s[i], unit=(u.hourangle, u.deg))
	ra_k=c.ra.deg
	dec_k=c.dec.deg
	out_ra.append(ra_k)
	out_dec.append(dec_k)
	w.write('circle('+str(ra_k)+'d, '+str(dec_k)+'d, 2\") #color=magenta width=2 \n')
w.close()

#write a new version of kenter catalog in fits format
cat=Table([out_ra,out_dec,cts_f,cts_s,cts_h,flux_f*1e-15,flux_s*1e-15,flux_h*1e-15],names=('RA','DEC','CTS_FULL','CTS_SOFT','CTS_HARD','FLUX_FULL','FLUX_SOFT','FLUX_HARD'))
cat.write(wd+'xbootes_kenter_cat0.fits',format='fits',overwrite=True) 

sys.exit()
'''
bins=np.logspace(np.log10(1e-16),np.log10(1e-12),20)

plt.figure()
plt.hist(flux,bins=bins,histtype='step',linewidth=3,color='red')
plt.xscale('log')
plt.xlabel(r'0.5-7 keV flux (erg cm$^{-2}$ s$^{-1}$)',fontsize=13)
plt.ylabel('N',fontsize=13)
plt.tight_layout()
plt.savefig(wd+'cdwfs_broad_flux_distr.pdf',format='pdf',dpi=1000)
'''

plt.figure()
plt.plot(det,cts,'k.')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('DET_ML',fontsize=13)
plt.ylabel(r'Net cts',fontsize=13)
plt.axhline(y=6,color='red')
plt.axvline(x=8,color='green')
#plt.axis([1,1000,8e-18,1e-12])
#plt.axis([5e-16,1e-12,5e-16,1e-12])
plt.tight_layout()
plt.show()
#plt.savefig(wd+'simul_broad_flux_flux.pdf',format='pdf',dpi=1000)
'''
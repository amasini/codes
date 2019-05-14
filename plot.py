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
import scipy.special
from scipy.integrate import quad
import time
from scipy.interpolate import interp1d
import os
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from scipy import stats

def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*xx*3600)**2 +((pointa[1]-pointb[1])*3600)**2)
    
def integrand(x,a):
    return x**(-a+1)

def gauss(x,mu,sigma):
	g=np.exp(-(x-mu)**2/(2*sigma**2))
	return g
wd='/Users/alberto/Desktop/XBOOTES/'

#expmap=fits.open(wd+'new_mosaics_detection/cdwfs_broad_expomap_4reb.fits')
#exp=expmap[0].data

#with fits.open(wd+'psfmaps/cdwfs_broad_r90sq-x-exp_4reb-cp.fits', mode='update') as hdul:
#	# Change something in hdul.
#	hdul[0].data=hdul[0].data/(16.0*exp)
#	hdul[0].data[np.isnan(hdul[0].data)]=0.0
#	hdul.flush()  # changes are written back
#sys.exit()

cat=fits.open(wd+'new_mosaics_detection/cdwfs_merged_cat1.fits')
ctsf=cat[1].data['NET_F']
ctss=cat[1].data['NET_S']
ctsh=cat[1].data['NET_H']
poserr=cat[1].data['POS_ERR']
hr=cat[1].data['HR']
fluxf=cat[1].data['FLUX_F']
fluxs=cat[1].data['FLUX_S']
fluxh=cat[1].data['FLUX_H']
probf=cat[1].data['PROB_F']

poserr=poserr[fluxf>0]
ctsf=ctsf[ctsf>0]
ctss=ctss[ctss>0]
ctsh=ctsh[ctsh>0]
probf=probf[fluxf>0]

hrb=hr[fluxf>0]
fluxf=fluxf[fluxf>0]
fluxs=fluxs[fluxs>0]
fluxh=fluxh[fluxh>0]

fluxf2=fluxf[hrb!=-99]
hr2=hrb[hrb!=-99]

fig = plt.figure(tight_layout=True,figsize=[6,4])
gs = gridspec.GridSpec(2,3)
gs.update(wspace=0., hspace=0.)

ax1 = plt.subplot(gs[0:2, 0:2])
ax1.hexbin(fluxf2,hr2,xscale='log',bins='log',gridsize=30)
#ax1.xaxis.set_visible(False)
#ax1.yaxis.set_visible(False)
ax1.set_xscale('log')
ax1.set_xlabel(r'0.5-7 keV flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=13)
ax1.set_ylabel('HR',fontsize=13)
#ax1.legend()

# make the contour plot
#plt.figure()
#plt.plot(fluxf2,hr2,marker='.',c='gold',mec='k',linestyle='')
#plt.hexbin(fluxf2,hr2,xscale='log',bins='log',gridsize=30)
#plt.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],linewidths=3,colors='red',linestyles='solid')
#plt.xlabel(r'0.5-7 keV flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=13)
#plt.ylabel('HR',fontsize=13)
#plt.tick_params(which='major',labelsize=13)
#plt.xscale('log')
hr=hr[hr!=-99]

ax2 = plt.subplot(gs[0:2, 2:3])
ax2.hist(hr,bins=30,color='gold',edgecolor='k',orientation='horizontal')
#ax2.xaxis.set_visible(False)
ax2.set_xlabel('#')
ax2.yaxis.set_visible(False)
plt.show()

#plt.savefig(wd+'cdwfs_hr.pdf',format='pdf')
#hr=hr[hr!=-99]
#plt.figure()
#plt.hist(hr,bins=30,color='gold',edgecolor='k')
#plt.xlabel('HR',fontsize=13)
#plt.tick_params(which='major',labelsize=13)
#plt.xscale('log')
#plt.legend()
#plt.tight_layout()
#plt.show()
#plt.savefig(wd+'cdwfs_net_distr.pdf',format='pdf')

print('min(pos_err):',min(poserr),'max(pos_err):',max(poserr))
bins=np.logspace(np.log10(min(poserr)),np.log10(max(poserr)),50)
plt.figure()
plt.hist(poserr,bins=bins)
plt.xscale('log')
plt.xlabel(r'Pos err ["]',fontsize=13)
plt.show()

plt.figure()
plt.plot(fluxf,poserr,'b+')
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Pos err ["]',fontsize=13)
plt.xlabel(r'0.5-7 keV flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=13)
plt.tick_params(which='major',labelsize=13)
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(probf,fluxf,'k.')
plt.axvline(x=6e-5,linestyle='dashed')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Prob of being spurious',fontsize=13)
plt.ylabel(r'0.5-7 keV flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=13)
plt.axis([1e-30,1,5e-17,5e-13])
plt.tick_params(which='major',labelsize=13)
plt.tight_layout()
plt.show()
#plt.savefig(wd+'cdwfs_full_flux.pdf',format='pdf')

bins=np.logspace(np.log10(min(ctsf)),np.log10(max(ctsf)),30)
plt.figure()
plt.hist(ctsf,bins=bins,label='F')
plt.hist(ctss,bins=bins,label='S')
plt.hist(ctsh,bins=bins,label='H')
plt.xlabel('Net cts',fontsize=13)
plt.tick_params(which='major',labelsize=13)
plt.xscale('log')
plt.legend()
plt.tight_layout()
plt.show()
#plt.savefig(wd+'cdwfs_net_distr.pdf',format='pdf')
#sys.exit()

#with fits.open(wd+'new_mosaics_detection/cdwfs_hard_expomap_4reb.fits', mode='update') as #hdul:
#	# Change something in hdul.
#	hdul[0].data=hdul[0].data/16.0
#	hdul.flush()  # changes are written back

#sys.exit()
'''
cat=fits.open(wd+'new_mosaics_detection/xbootes_broad_cat0.fits')
ra=cat[1].data['RA']
dec=cat[1].data['DEC']
r90=cat[1].data['AV_R90']
tot=cat[1].data['TOT']

ra=ra[tot>=4]
dec=dec[tot>=4]
r90=r90[tot>=4]

ken=fits.open(wd+'xbootes_kenter+05.fits')
rak=ken[1].data['RAJ2000']
deck=ken[1].data['DEJ2000']
found=0
for i in range(len(ra)):
	my=[ra[i],dec[i]]
	at_least_one=0
	for j in range(len(rak)):
		xb=[rak[j],deck[j]]
		d=distance(my,xb)
		if d <= 1.1*r90[i]:
			if at_least_one==0:
				at_least_one=1
				found=found+1
print(len(ra),len(rak))
print(found,'matches')
			
sys.exit()
'''

obsid=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=1,dtype='str')
#obsid=['10450']
kt,apec,cxb=[],[],[]
for obs in obsid:
	kt_v=np.genfromtxt(wd+'data/'+obs+'/repro_new_asol/extract/fitparams.dat',skip_footer=4,unpack=True,usecols=0)
	cr=np.genfromtxt(wd+'data/'+obs+'/repro_new_asol/extract/fitparams.dat',skip_header=3,unpack=True,usecols=2)
	if kt_v < 72.0:
		kt.append(kt_v)
	else:
		print(kt_v,'check fit obsid',obs)
	apec.append(cr[0])
	cxb.append(cr[1])

exp=np.genfromtxt(wd+'apec_cr_soft_with-gauss-line.dat',unpack=True, skip_header=1,usecols=3)
(mjd,bkg,e_bkg)=np.genfromtxt(wd+'avg_bkg_new.dat',unpack=True,usecols=[1,2,3],skip_header=1)

#bins=np.linspace(min(mjd),max(mjd),7)
bins=[52100,52440,52750,53700,54000,54300,56000,57000,57800,58400]
nyears=int((max(bins)-min(bins))/365.)

#bin_means, bin_edges, binnumber = stats.binned_statistic(mjd,apec, statistic='median', bins=nyears)

apec=np.array(apec)
n,be=np.histogram(mjd,bins=nyears)
bc=list((be[i+1]+be[i])/2. for i in range(len(be)-1))
median,lo,up=[],[],[]
for i in range(len(be)-1):
	newapec=apec[(mjd>=be[i]) & (mjd<be[i+1])]
	if len(newapec)==0:
		median.append(0)
		lo.append(0)
		up.append(0)
	else:
		median.append(np.median(newapec))
		lo.append(np.median(newapec)-np.percentile(newapec,16))
		up.append(np.percentile(newapec,84)-np.median(newapec))

binexp=np.logspace(np.log10(min(exp)),np.log10(max(exp)),11)
n,be=np.histogram(exp,bins=binexp)
bc_exp=list((be[i+1]+be[i])/2. for i in range(len(be)-1))
median_exp,lo_exp,up_exp=[],[],[]
for i in range(len(be)-1):
	newapec=apec[(exp>=be[i]) & (exp<be[i+1])]
	if len(newapec)==0:
		median_exp.append(0)
		lo_exp.append(0)
		up_exp.append(0)
	else:
		median_exp.append(np.median(newapec))
		lo_exp.append(np.median(newapec)-np.percentile(newapec,16))
		up_exp.append(np.percentile(newapec,84)-np.median(newapec))

cxb=np.array(cxb)
n,be=np.histogram(mjd,bins=nyears)
bc_cxb=list((be[i+1]+be[i])/2. for i in range(len(be)-1))
median_cxb,lo_cxb,up_cxb=[],[],[]
for i in range(len(be)-1):
	newcxb=cxb[(mjd>=be[i]) & (mjd<be[i+1])]
	if len(newcxb)==0:
		median_cxb.append(0)
		lo_cxb.append(0)
		up_cxb.append(0)
	else:
		median_cxb.append(np.median(newcxb))
		lo_cxb.append(np.median(newcxb)-np.percentile(newcxb,16))
		up_cxb.append(np.percentile(newcxb,84)-np.median(newcxb))

n,be=np.histogram(exp,bins=binexp)
bc_cxb_exp=list((be[i+1]+be[i])/2. for i in range(len(be)-1))
median_cxb_exp,lo_cxb_exp,up_cxb_exp=[],[],[]
for i in range(len(be)-1):
	newcxb=cxb[(exp>=be[i]) & (exp<be[i+1])]
	if len(newcxb)==0:
		median_cxb_exp.append(0)
		lo_cxb_exp.append(0)
		up_cxb_exp.append(0)
	else:
		median_cxb_exp.append(np.median(newcxb))
		lo_cxb_exp.append(np.median(newcxb)-np.percentile(newcxb,16))
		up_cxb_exp.append(np.percentile(newcxb,84)-np.median(newcxb))

#tot=cxb+apec

#plt.figure()
#plt.hist(kt,bins=10)
#plt.show()

#plt.figure()
#plt.plot(exp,kt,'bx')
#plt.xscale('log')
#plt.show()
#plt.figure()
#plt.plot(mjd,kt,'rx')
#plt.show()

plt.figure()
#plt.errorbar(mjd,bkg,yerr=e_bkg,color='red',fmt='.',label='9-12 keV')
plt.plot(exp,cxb,marker='o',markerfacecolor="None",markeredgecolor='green',markeredgewidth=1,label='CXB',linestyle="None")
plt.plot(exp,apec,marker='o',markerfacecolor="None",markeredgecolor='gold',markeredgewidth=1,label='APEC',linestyle="None")
plt.errorbar(bc_exp,median_exp,yerr=[lo_exp,up_exp],fmt='s',color='red')
plt.errorbar(bc_cxb_exp,median_cxb_exp,yerr=[lo_cxb_exp,up_cxb_exp],fmt='s',color='blue')
#plt.errorbar(mjd,tot,yerr=e_bkg,color='black',fmt='+',ms=10)
#plt.errorbar(mjd,apec_sb,yerr=e_apec_sb,color='green',fmt='.',label='APEC')
#plt.errorbar(mjd,pl_sb,yerr=e_pl_sb,color='blue',fmt='.',label='PL')
#plt.errorbar(mjd,total_sb,yerr=e_total_sb,color='black',fmt='.',label='APEC+PL')
#plt.errorbar(mjd,rescaled_bkg,yerr=e_rescaled_bkg,color='gray',fmt='.',label=r'9-12keV $\rightarrow$ 0.5-1 keV')
plt.xlabel('Exp [s]')
plt.ylabel('Cts/s')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.tight_layout()

plt.figure()
#plt.errorbar(mjd,bkg,yerr=e_bkg,color='red',fmt='.',label='9-12 keV')
plt.plot(mjd,cxb,marker='o',markerfacecolor="None",markeredgecolor='green',markeredgewidth=1,label='CXB',linestyle="None")
plt.plot(mjd,apec,marker='o',markerfacecolor="None",markeredgecolor='gold',markeredgewidth=1,label='APEC',linestyle="None")
plt.errorbar(bc,median,yerr=[lo,up],fmt='s',color='red')
plt.errorbar(bc_cxb,median_cxb,yerr=[lo_cxb,up_cxb],fmt='s',color='blue')
#plt.errorbar(mjd,tot,yerr=e_bkg,color='black',fmt='+',ms=10)
#plt.errorbar(mjd,apec_sb,yerr=e_apec_sb,color='green',fmt='.',label='APEC')
#plt.errorbar(mjd,pl_sb,yerr=e_pl_sb,color='blue',fmt='.',label='PL')
#plt.errorbar(mjd,total_sb,yerr=e_total_sb,color='black',fmt='.',label='APEC+PL')
#plt.errorbar(mjd,rescaled_bkg,yerr=e_rescaled_bkg,color='gray',fmt='.',label=r'9-12keV $\rightarrow$ 0.5-1 keV')
plt.xlabel('MJD')
plt.ylabel('Cts/s')
plt.yscale('log')
plt.legend()
plt.tight_layout()
plt.show()

sys.exit()

#bins=np.logspace(np.log10(min(fluxf)),np.log10(max(fluxf)),20)
#plt.figure()
#plt.hist(fluxf,bins=bins)
#plt.hist(fluxs,bins=bins)
#plt.hist(fluxh,bins=bins)
#plt.xscale('log')
#plt.show()
cut=np.logspace(np.log10(5e-5),np.log10(1e-1),15)
sp_frac_f=np.genfromtxt(wd+'cdwfs_broad_sp-frac.dat',unpack=True,usecols=1)
sp_frac_s=np.genfromtxt(wd+'cdwfs_soft_sp-frac.dat',unpack=True,usecols=1)
sp_frac_h=np.genfromtxt(wd+'cdwfs_hard_sp-frac.dat',unpack=True,usecols=1)

plt.figure()

plt.plot(cut,sp_frac_f,'ko-',label='Full band')
plt.plot(cut,sp_frac_s,'go-',label='Soft band')
plt.plot(cut,sp_frac_h,'ro-',label='Hard band')

plt.axhline(y=1.0)
plt.axhline(y=3.0)
plt.xscale('log')
plt.xlabel('Prob cut',fontsize=13)
plt.ylabel('Sp fraction (%)',fontsize=13)
plt.tick_params(which='major',labelsize=13)
plt.legend()
plt.tight_layout()
plt.savefig(wd+'cdwfs_sp-frac.pdf',format='pdf')
#plt.show()

sys.exit()

'''
x=np.linspace(0,5,101)
sigma=np.sqrt(0.2**2+1.0**2)
y=2*np.pi*x/(np.sqrt(2*np.pi)*sigma)*gauss(x,0,sigma)
y2=1./(np.sqrt(2*np.pi)*sigma)*gauss(x,0,sigma)
normy=y/np.sum(y)
normy2=y2/np.sum(y2)
plt.figure()
plt.plot(x,normy,'r-')
plt.plot(x,normy2,'b-')
plt.show()

sys.exit()
band='I' # I, K, 3.6mu

nmI=np.genfromtxt(wd+'LR_nm_'+band+'.dat',unpack=True)
qmI=np.genfromtxt(wd+'LR_qm_'+band+'.dat',unpack=True)
binwidth=0.5 # mag

nbinsI=int((max(nmI)-min(nmI))/binwidth)
qbinsI=int((max(qmI)-min(qmI))/binwidth)
(n,b)=np.histogram(nmI,bins=nbinsI)
(q,bq)=np.histogram(qmI,bins=qbinsI)

#print(n,b) # n is the number of sources in each bin of magnitude 
r=5.0
R=30.0
A=np.pi*(R**2-r**2)
nnew=(n/A)*np.pi # number of sources per arcsec, for the whole 6963 sources, times the area of a circular 1"-radius region (A=pi*R^2)
bnew=list((b[i+1]+b[i])/2. for i in range(len(b)-1))
bqnew=list((bq[i+1]+bq[i])/2. for i in range(len(bq)-1))
#sys.exit()

# Interpolate the two distributions
xvals=np.linspace(11,30,191) # Magnitudes
qvals = np.interp(xvals, bqnew, q)
nvals = np.interp(xvals, bnew, nnew)
diff=qvals-nvals

'''
# Cubic interpolation
#fq2 = interp1d(bqnew, q, kind='cubic')
#fn2 = interp1d(bnew, nnew, kind='cubic')

#xvalsq2=np.linspace(10,20,(20-10)*10+1) # Magnitudes
#xvalsn2=np.linspace(10,20,(20-10)*10+1) # Magnitudes
#qvals2=fq2(xvalsq2)
#nvals2=fn2(xvalsn2)
#diff2=qvals2-nvals2
'''
w=open(wd+'LR_qprimom_'+band+'.dat','w')
for k in range(len(diff)):
	w.write(str(xvals[k])+' \t '+str(diff[k])+'\n')
w.close()

plt.figure()
#plt.hist(qmI,bins=qbinsI)
plt.plot(bqnew,q,'k+')
#plt.plot(xvals,qvals,'b.')
plt.plot(bnew,nnew,'r+')
#plt.plot(xvals,nvals,'g.')
plt.plot(xvals,diff,'y-')
#plt.plot(xvalsq2,diff2,'b--')
plt.xlabel('Magnitude')
plt.ylabel('# of counterparts per arcsec')
plt.show()

diff2=diff/np.max(np.sum(diff)) # this is normalized such as np.sum(diff2)=1
if band == 'I':
	diff3=diff2*(3700./6963.)
elif band == 'K':
	diff3=diff2*(1303./6963.)
elif band == '3.6mu':
	diff3=diff2*(3395./6963.)
print(diff3,np.sum(diff3))

sys.exit()

'''

'''
# Open Chung+14 catalog
cat3=fits.open(wd+'xbootes_chung+14.fits')
ra_chu=cat3[1].data['RAJ2000']
dec_chu=cat3[1].data['DEJ2000']
w=open(wd+'chung+14.reg','w')
for k in range(len(ra_chu)):
	w.write('circle('+str(ra_chu[k])+'d,'+str(dec_chu[k])+'d,2\") #color=magenta \n')
w.close()
sys.exit()

cat=fits.open(wd+'cdwfs_merged_cat1.fits')
ra=cat[1].data['RA']
dec=cat[1].data['DEC']
r90=cat[1].data['R90_S']
probf=cat[1].data['PROB_F']
probs=cat[1].data['PROB_S']
probh=cat[1].data['PROB_H']
cutf,cuts,cuth=6e-5,1.4e-4,6e-5
fsh,fs,fh,sh,f,s,h=0,0,0,0,0,0,0
for i in range(len(probf)):
	if (probf[i] != 9999.0 and probf[i] <= cutf):
		if (probs[i] != 9999.0 and probs[i] <= cuts):
			if (probh[i] != 9999.0 and probh[i] <= cuth):
				fsh=fsh+1
			else:
				fs=fs+1
		else:
			if (probh[i] != 9999.0 and probh[i] <= cuth):
				fh=fh+1
			else:
				f=f+1
	else:
		if (probs[i] != 9999.0 and probs[i] <= cuts):
			if (probh[i] != 9999.0 and probh[i] <= cuth):
				sh=sh+1
				print(ra[i],dec[i])
				sys.exit()
			else:
				s=s+1
		else:
			if (probh[i] != 9999.0 and probh[i] <= cuth):
				h=h+1
			else:
				print('something wrong',totf[i],tots[i],toth[i])
print('fsh:',fsh)
print('fs:',fs)
print('fh:',fh)
print('sh:',sh)
print('f:',f)
print('s:',s)
print('h:',h)
print(fsh+fs+fh+sh+f+s+h)
#print(len(totf),len(totf[totf>=30.]))
sys.exit()

gamma=np.arange(1.3,2.4,0.1)
cf_f=[6.243,6.355,6.456,6.544,6.617,6.674,6.712,6.731,6.731,6.709,6.668]
cf_s=[9.592,9.391,9.186,8.976,8.763,8.547,8.328,8.107,7.885,7.662,7.438]
cf_h=[4.798,4.864,4.930,4.996,5.061,5.126,5.191,5.255,5.318,5.381,5.442]
cf_f=np.array(cf_f)*1e10
cf_s=np.array(cf_s)*1e10
cf_h=np.array(cf_h)*1e10

fluxrat_s=[0.3016,0.3295,0.3587,0.3890,0.4203,0.4524,0.4849,0.5176,0.5502,0.5825,0.6141]
fluxrat_h=1-np.array(fluxrat_s)

xvals = np.linspace(1.3, 2.3, 101)
yinterp_f = np.interp(xvals, gamma, cf_f)
yinterp_s = np.interp(xvals, gamma, cf_s)
yinterp_h = np.interp(xvals, gamma, cf_h)

fluxinterp_s = np.interp(xvals, gamma, fluxrat_s)
fluxinterp_h = 1-fluxinterp_s

mu=1.8
sigma=0.2
p=[]
for ii in range(len(xvals)):
	p.append(gauss(xvals[ii],mu,sigma))
p=np.array(p)
new=p/np.sum(p)
#plt.figure()
#plt.plot(xvals,new,'k-')
#plt.axvline(mu)
#plt.axvline(mu-sigma)
#plt.axvline(mu+sigma)
#plt.show()


#(full_flux,ra,dec)=np.genfromtxt(wd+'poiss_rand_lehmerx20.dat',unpack=True, skip_header=1,usecols=[0,1,2])
#w=open(wd+'poiss_rand_lehmerx20_NEW.dat','w')
#w.write('Full flux \t RA \t DEC \t Gamma\n')
#for i in range(len(ra)):
#	random_gamma=np.random.choice(xvals, p=new)
#	w.write(str(full_flux[i])+' \t '+str(ra[i])+' \t '+str(dec[i])+' \t '+str(random_gamma)+'\n')
#w.close()
#sys.exit()

#print(random_gamma,yinterp_f[xvals==random_gamma],yinterp_s[xvals==random_gamma],yinterp_h[xvals==random_gamma])
#print(fluxinterp_s[xvals==random_gamma],1-fluxinterp_s[xvals==random_gamma])
plt.figure()
plt.plot(gamma, cf_f,'ko',label='Full')
plt.plot(xvals, yinterp_f,'k--')
plt.plot(gamma, cf_s,'go',label='Soft')
plt.plot(xvals, yinterp_s,'g--')
plt.plot(gamma, cf_h,'ro',label='Hard')
plt.plot(xvals, yinterp_h,'r--')
plt.xlabel('Gamma')
plt.xlabel('Flux to CRate CF')
plt.legend()
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(gamma, fluxrat_s,'go',label='Soft')
plt.plot(xvals, fluxinterp_s,'g--')
plt.plot(gamma, fluxrat_h,'ro',label='Hard')
plt.plot(xvals, fluxinterp_h,'r--')
plt.xlabel('Gamma')
plt.ylabel('Ratio to 0.5-7 keV flux')
plt.legend()
plt.tight_layout()
plt.show()
#sys.exit()
'''

#(obs,apec_cr,pl_cr,exp,pixarea,apec_sb,apec,pl_sb,pl)=np.genfromtxt(wd+'apec_cr_soft_with-gauss-line.dat',unpack=True, skip_header=1)
(cxb,apec)=np.genfromtxt(wd+'BACK_soft.dat',unpack=True, skip_header=1,usecols=[1,2])
(mjd,bkg,e_bkg)=np.genfromtxt(wd+'avg_bkg_new.dat',unpack=True,usecols=[1,2,3],skip_header=1)
#errs=e_bkg/bkg
#err=np.median(errs)
#e_apec_sb=err*apec_sb
#e_pl_sb=err*pl_sb

#total_sb=apec_sb+pl_sb
#e_total_sb=np.sqrt(e_apec_sb**2+e_pl_sb**2)
#rescaled_bkg=0.4*bkg
#rescaled_bkg=0.13*bkg
#e_rescaled_bkg=0.13*e_bkg

#cxb=pl-ib
#diff=(total_sb-rescaled_bkg)/total_sb

#plt.figure()
#plt.plot(exp,cxb,'b.')
#plt.xscale('log')
#plt.show()
tot=cxb+apec

plt.figure()
#plt.errorbar(mjd,bkg,yerr=e_bkg,color='red',fmt='.',label='9-12 keV')
plt.errorbar(mjd,cxb,yerr=e_bkg,color='green',fmt='.',label='CXB')
plt.errorbar(mjd,apec,yerr=e_bkg,color='gold',fmt='.',label='APEC')
plt.errorbar(mjd,tot,yerr=e_bkg,color='black',fmt='+',ms=10)
#plt.errorbar(mjd,apec_sb,yerr=e_apec_sb,color='green',fmt='.',label='APEC')
#plt.errorbar(mjd,pl_sb,yerr=e_pl_sb,color='blue',fmt='.',label='PL')
#plt.errorbar(mjd,total_sb,yerr=e_total_sb,color='black',fmt='.',label='APEC+PL')
#plt.errorbar(mjd,rescaled_bkg,yerr=e_rescaled_bkg,color='gray',fmt='.',label=r'9-12keV $\rightarrow$ 0.5-1 keV')
plt.xlabel('MJD')
plt.ylabel('Cts/s/pixel')
plt.yscale('log')
plt.legend()
plt.tight_layout()
plt.show()

#plt.figure()
#plt.plot(mjd,diff,'k.')
#plt.yscale('log')
#plt.show()

sys.exit()
'''
obsid=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=1,dtype='str')
for i in range(len(obsid)):
	if len(obsid[i]) == 4:
		stem='0'+obsid[i]
	elif len(obsid[i]) == 3:
		stem='00'+obsid[i]
	elif len(obsid[i]) == 5:
		stem=obsid[i]
		
	s.call('ds9 '+wd+'data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img -scale log -scale mode 90 -zoom 0.65 -region '+wd+'data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_broad_src-bkg.reg',shell=True)
sys.exit()
'''


#cut=np.logspace(np.log10(1e-4),np.log10(1e-2),10)
#cut2=np.logspace(np.log10(1e-5),np.log10(1e-4),10)
cut=np.logspace(np.log10(5e-5),np.log10(1e-1),15)

#sp_frac=np.array([41./3999.,58./4145.,68./4288.,83./4437.,98./4571.,113./4677.,131./4780.,146./4878.,169./4974.,191./5039.])*100.
#sp_frac_f=np.array([95./6345.,118./6526.,147./6691.,175./6854.,198./6988.,223./7100.,259./7230.,292./7351.,336/7463.,373./7547.])*100.
#sp_frac_f2=np.array([42./5579.,47./5665.,50./5741.,56./5812.,59./5879.,64./5985.,71./6073.,81./6176.,90/6257.,95./6345.])*100.

#sp_frac_f_new=np.array([16./5487.,25./5710.,34./5968.,50./6206.,68./6455.,89./6688.,118./6906.,147./7098.,171./7228.,201./7355.])*100.
#sp_frac_s_new=np.array([8./3228.,10./3420.,17./3611.,25./3775.,43./3944.,62./4102.,78./4231.,96./4332.,116./4410.,132./4466.])*100.
#sp_frac_h_new=np.array([15./3276.,22./3468.,33./3672.,43./3874.,57./4079.,69./4264.,95./4448.,119./4588.,157./4719.,192./4821.])*100.

sp_frac_f=np.genfromtxt(wd+'cdwfs_broad_sp-frac.dat',unpack=True,usecols=1)
sp_frac_s=np.genfromtxt(wd+'cdwfs_soft_sp-frac.dat',unpack=True,usecols=1)
sp_frac_h=np.genfromtxt(wd+'cdwfs_hard_sp-frac.dat',unpack=True,usecols=1)

#add last 3 points to the full band with cut4=np.logspace(np.log10(1e-2),np.log10(3e-2),3) 
#cut4=np.logspace(np.log10(1e-2),np.log10(3e-2),3)
#sp_frac_f_new2=np.array([201./7355.,234./7425.,258./7483.])*100

plt.figure()

plt.plot(cut,sp_frac_f,'k-')
plt.plot(cut,sp_frac_f,'ko',label='Full band')
plt.plot(cut,sp_frac_s,'g-')
plt.plot(cut,sp_frac_s,'go',label='Soft band')
plt.plot(cut,sp_frac_h,'r-')
plt.plot(cut,sp_frac_h,'ro',label='Hard band')

plt.axhline(y=1.0)
plt.axhline(y=3.0)
plt.xscale('log')
plt.xlabel('Prob cut')
plt.ylabel('Sp fraction (%)')
plt.legend()
plt.tight_layout()
plt.show()
sys.exit()
'''
d=np.genfromtxt(wd+'sim_all/cdwfs_broad_sim_poiss_matched.dat',unpack=True,usecols=5)
d2=np.genfromtxt(wd+'sim_all/cdwfs_soft_sim_poiss_matched.dat',unpack=True,usecols=5)
d3=np.genfromtxt(wd+'sim_all/cdwfs_hard_sim_poiss_matched.dat',unpack=True,usecols=5)
print(np.median(d),np.percentile(d, 68),np.percentile(d, 90),np.percentile(d, 99))
print(np.median(d2),np.percentile(d2, 68),np.percentile(d2, 90),np.percentile(d2, 99))
print(np.median(d3),np.percentile(d3, 68),np.percentile(d3, 90),np.percentile(d3, 99))

val,b=np.histogram(d,bins=50)
val2,b2=np.histogram(d2,bins=50)
val3,b3=np.histogram(d3,bins=50)
val=val.astype(float)
val2=val2.astype(float)
val3=val3.astype(float)

bincenters=list((b[i+1]+b[i])/2. for i in range(len(b)-1))
bincenters2=list((b2[i+1]+b2[i])/2. for i in range(len(b)-1))
bincenters3=list((b3[i+1]+b3[i])/2. for i in range(len(b)-1))

cum=np.cumsum(val)/np.sum(val)
cum2=np.cumsum(val2)/np.sum(val2)
cum3=np.cumsum(val3)/np.sum(val3)

plt.figure()
plt.hist(d,bins=50,alpha=0.7,normed=True,label='Full')
plt.hist(d2,bins=50,alpha=0.7,normed=True,label='Soft')
plt.hist(d3,bins=50,alpha=0.7,normed=True,label='Hard')
plt.xlabel('d ["]')
plt.ylabel('N')
plt.legend()
plt.axis([-0.5,6,0,0.8])
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(bincenters,cum,'k-')
plt.plot(bincenters2,cum2,'r-')
plt.plot(bincenters3,cum3,'b-')
plt.axhline(y=0.68)
plt.axhline(y=0.90)
plt.axhline(y=0.99)
plt.xlabel('d ["]')
plt.ylabel('Cumulative fraction')
plt.axis([0,6,0,1])
plt.tight_layout()
plt.show()
'''
band='hard'
if band=='broad':
	cut97=1.4e-2
	cut99=2e-4
elif band=='soft':
	cut97=1e-2
	cut99=2e-4
elif band=='hard':
	cut97=3.5e-3
	cut99=1e-4

cat=fits.open(wd+'mosaic_'+band+'_cat1_3.fits')
ra=cat[1].data['RA']
dec=cat[1].data['DEC']
prob=cat[1].data['PROB']
flux=cat[1].data['FLUX']
r90=cat[1].data['AV_R90']
cts=cat[1].data['TOT']

print(band,len(ra))
print(len(ra[prob<=cut97]))
print(len(ra[prob<=cut99]))

ra=ra[prob<=cut97]
dec=dec[prob<=cut97]
r90=r90[prob<=cut97]

w=open(wd+'cdwfs_'+band+'_cat1_3.reg','w')
for i in range(len(ra)):
	w.write('circle('+str(ra[i])+'d, '+str(dec[i])+'d, '+str(r90[i])+'\") #color=yellow \n')
w.close()

plt.figure()
plt.plot(prob,cts,'k.')
plt.xscale('log')
plt.yscale('log')
plt.show()
sys.exit()


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
band='soft'
#cat=fits.open(wd+'cdwfs_broad_cat0.fits')
#cat=fits.open(wd+'murray_sens/xbootes_unique_2.fits')
#cat=fits.open(wd+'murray_sens/xbootes_broad_cat0.fits')
cat=fits.open(wd+'sim_all/mosaic_'+band+'_sim_poiss_cat0_3.fits') #7881 SIMULATIONS [#7211 sources] BLUE
ra=cat[1].data['RA']
dec=cat[1].data['DEC']
prob=cat[1].data['PROB']
net=cat[1].data['NET']

catb=fits.open(wd+'mosaic_'+band+'_cat0_3.fits') #8800 sources RED
rab=catb[1].data['RA']
decb=catb[1].data['DEC']
probb=catb[1].data['PROB']
#detmlb=catb[1].data['DETML']
netb=catb[1].data['NET']
#$r=cat[1].data['AV_R90']
#cts=cat[1].data['NET']
#flux=cat[1].data['FLUX']

#inpflux=cat[1].data['CTRP_FLUX']
#prob=np.e**(-detml)
#probb=np.e**(-detmlb)

print(len(ra[prob<1e-4]))
print(len(rab[probb<1e-4]))

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

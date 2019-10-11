import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
from scipy.stats import kde

wd='/Users/alberto/Desktop/XBOOTES/'

# Open the matched master catalog (-cp version contains only one CDWFS source per line) 
# 7234 sources 
#cat=fits.open('/Users/alberto/Downloads/nway-master/cdwfs_I-Ks-3.6-cp.fits')
cat=fits.open(wd+'CDWFS_I-Ks-3.6_v2.fits')
data=cat[1].data
cols=cat[1].columns
names=cols.names
pany=data['p_any']
poserr0=data['CHA_POS_ERR']
xb_id0=data['CHA_XB_ID']
fluxf0=data['CHA_FLUX_F']
efluxfp0=data['CHA_E_FLUX_F_+']
fluxs0=data['CHA_FLUX_S']
efluxsp0=data['CHA_E_FLUX_S_+']
fluxh0=data['CHA_FLUX_H']
efluxhp0=data['CHA_E_FLUX_H_+']
hr0=data['CHA_HR']
ehr0p=data['CHA_E_HR_+']
ehr0n=data['CHA_E_HR_-']
imag0=data['NDWFS_MAG_AUTO']
kmag0=data['IBIS_MAG_BEST']
spmag0=data['SDWFS_ch1_ma']
sep0=data['Separation_NDWFS_CHA']
sep0k=data['Separation_IBIS_CHA']
sep0sp=data['Separation_SDWFS_CHA']
zspec0=data['zsp']
zph_g0=data['zph_G']
chi_g0=data['chi2_G']
zph_ga0=data['zph_G+A']
chi_ga0=data['chi2_G+A']
ebv0=data['E_B-V']
chi_s0=data['chi2_S']
cat.close()
print('='*30)
print('The full CDWFS catalog has',len(data),'sources.')

# X-ray Flux distributions
fluxf=fluxf0[efluxfp0!=0]
fluxs=fluxs0[efluxsp0!=0]
fluxh=fluxh0[efluxhp0!=0]
bins=np.logspace(np.log10(5e-16),np.log10(2e-13),30)
f,ax=plt.subplots(3,sharex=True)
ax[0].hist(fluxf,bins=bins,histtype='step',linewidth=3,color='k')
ax[0].tick_params(which='major',direction='inout',length=7,labelsize=13)
ax[0].tick_params(which='minor',direction='inout',length=4,labelsize=13)
ax[0].annotate('0.5-7 keV',xy=(7e-14,500),fontsize=13)

ax[1].hist(fluxs,bins=bins,histtype='step',linewidth=3,color='r')
ax[1].tick_params(which='major',direction='inout',top=True,length=7,labelsize=13)
ax[1].tick_params(which='minor',direction='inout',top=True,length=4,labelsize=13)
ax[1].annotate('0.5-2 keV',xy=(7e-14,400),fontsize=13)

ax[2].hist(fluxh,bins=bins,histtype='step',linewidth=3,color='b')
ax[2].set_xlabel(r'Flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=12)
ax[2].set_xscale('log')
ax[2].tick_params(which='major',direction='inout',top=True,length=8,labelsize=13)
ax[2].tick_params(which='minor',direction='inout',top=True,length=4,labelsize=13)
ax[2].annotate('2-7 keV',xy=(7e-14,400),fontsize=13)
plt.subplots_adjust(hspace=0)
plt.show()
#plt.savefig('/Users/alberto/Desktop/cdwfs_flux.pdf',format='pdf')


p_any_cut=0.12 # 0.54 (0.12 as of 04-Sep-19) needed to have <5% false associations

print('A total of',len(pany[pany>0]),' (p_any>0) opt-NIR associations found.')

poserr=poserr0[pany>p_any_cut]
xb_id=xb_id0[pany>p_any_cut]
fluxf=fluxf0[pany>p_any_cut]
efluxfp=efluxfp0[pany>p_any_cut]
imag=imag0[pany>p_any_cut]
kmag=kmag0[pany>p_any_cut]
spmag=spmag0[pany>p_any_cut]
sep=sep0[pany>p_any_cut]
sepk=sep0k[pany>p_any_cut]
sepsp=sep0sp[pany>p_any_cut]
zsp=zspec0[pany>p_any_cut]
zph_g=zph_g0[pany>p_any_cut]
zph_ga=zph_ga0[pany>p_any_cut]
ebv=ebv0[pany>p_any_cut]
chi_g=chi_g0[pany>p_any_cut]
chi_ga=chi_ga0[pany>p_any_cut]
chi_s=chi_s0[pany>p_any_cut]
hr1=hr0[pany>p_any_cut]
ehr1p=ehr0p[pany>p_any_cut]
ehr1n=ehr0n[pany>p_any_cut]
pany=pany[pany>p_any_cut]

print('A total of',len(imag),'ROBUST (p_any>'+str(p_any_cut)+') opt-NIR associations found.')

# Plots of STARS with HR and photometric label from Chung
what=np.full_like(chi_g,'A',dtype='str')
what[(chi_g<chi_ga) & (chi_g<chi_s)]='G'
what[(chi_s<chi_ga) & (chi_s<chi_g)]='S'
print('Based on minimum Chi^2 from Chung+14, we have',len(what[what=='A']),'AGN,',len(what[what=='G']),'GAL, and',len(what[what=='S']),'STARS.')

chi_galoragn=chi_g
chi_galoragn[chi_ga<chi_g]=chi_ga[chi_ga<chi_g]

plt.figure()
plt.scatter(chi_galoragn[(chi_galoragn!=-99) & (pany<p_any_cut)],chi_s[(chi_s!=-99) & (pany<p_any_cut)], c=hr1[(chi_s!=-99) & (pany<p_any_cut)], cmap='viridis',alpha=0.2)
plt.scatter(chi_galoragn[(chi_galoragn!=-99) & (pany>p_any_cut)],chi_s[(chi_s!=-99) & (pany>p_any_cut)], c=hr1[(chi_s!=-99) & (pany>p_any_cut)], cmap='viridis')
plt.plot(np.arange(0,400),np.arange(0,400),'r--')
plt.annotate('GAL/AGN',xy=(0.065,0.3),rotation=45)
plt.annotate('STARS',xy=(150,150),rotation=45)
#plt.axis([-5,55,-5,55])
plt.axis([5e-2,400,5e-2,400])
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'$\chi^2_{\rm STAR}$',fontsize=12)
plt.xlabel(r'$\chi^2_{\rm AGN/GAL}$',fontsize=12)
plt.xticks(ticks=[0.1,1,10,100],labels=[0.1,1,10,100])
plt.yticks(ticks=[0.1,1,10,100],labels=[0.1,1,10,100])
plt.tick_params(axis='both', which='major', labelsize=12)
cbar = plt.colorbar()
#cbar.ax.set_yticklabels(['0','1','2','>3'])
cbar.set_label('HR', rotation=270)
plt.tight_layout()
plt.show()


hr_g=hr1[what=='G']
hr_a=hr1[what=='A']
hr_s=hr1[what=='S']

#plt.figure()
#plt.hist(hr_a,bins=20,color='b',histtype='step',linewidth=3,label='AGN',density=1)
#plt.hist(hr_g,bins=20,color='g',histtype='step',linewidth=3,label='GAL',density=1)
#plt.hist(hr_s,bins=20,color='r',histtype='step',linewidth=3,label='STAR',density=1)
#zstar=zsp[what=='S']
#print(zstar[zstar!=-99])
#plt.hist(zstar[zstar!=-99],bins=20,color='r',histtype='step',linewidth=3,label='STAR')
#plt.legend()
#plt.show()


ztot,specz,photz,zflag,sp,ph=[],[],[],[],[],[]
for i in range(len(zsp)):
	if zsp[i] != -99:
		ztot.append(zsp[i])
		specz.append(zsp[i])
		zflag.append(1)
		sp.append(1)
		ph.append(0)
	elif zph_ga[i] != -99:
		if what[i] == 'A':
			ztot.append(zph_ga[i])
			photz.append(zph_ga[i])
			zflag.append(1)
			sp.append(0)
			ph.append(1)
		elif what[i] =='G':
			ztot.append(zph_g[i])
			photz.append(zph_g[i])
			zflag.append(1)
			sp.append(0)
			ph.append(1)
		else:
			zflag.append(0)
			sp.append(0)
			ph.append(0)
	else:
		zflag.append(0)
		sp.append(0)
		ph.append(0)
ztot=np.array(ztot)
photz=np.array(photz)
specz=np.array(specz)
zflag=np.array(zflag)
sp=np.array(sp)
ph=np.array(ph)
hr=hr1[zflag==1]
ehrp=ehr1p[zflag==1]
ehrn=ehr1n[zflag==1]
print(len(ztot),' redshifts:',len(specz),'spec-z,',len(photz),'photo-z in the robust subsample.')
print('Redshifts range:',min(ztot),max(ztot))

#bins=np.linspace(0,5,30)
#plt.figure()
#plt.hist(ztot,bins=bins,histtype='step',linewidth=3,color='k')
#plt.hist(specz,bins=bins,histtype='step',linewidth=2,color='lime',label='Spec-z')
#plt.hist(photz,bins=bins,histtype='step',linewidth=2,color='blueviolet',linestyle='dashed',label='hoto-z')
#plt.ylabel('#',fontsize=12)
#plt.xlabel('Redshift',fontsize=12)
#plt.yscale('log')
#plt.yticks(ticks=[1,10,100],labels=[1,10,100])
#plt.tick_params(axis='both', which='major', labelsize=12)
#plt.legend()
#plt.savefig('/Users/alberto/Desktop/cdwfs_z.pdf',format='pdf',dpi=500)
#plt.show()

################################# HR-z
fig = plt.figure(tight_layout=True,figsize=[8,8])
gs = gridspec.GridSpec(3,3)
gs.update(wspace=0., hspace=0.)

pany_z=pany[zflag==1]
pany_zsp=pany[(zflag==1) & (sp==1)]
pany_zph=pany[(zflag==1) & (ph==1)]

ax2 = plt.subplot(gs[0:1, 0:2])
bins=np.linspace(0,5,30)
n,b=np.histogram(ztot,bins=bins)
bcen=list((b[i]+b[i+1])/2. for i in range(len(b)-1))
cum=np.cumsum(n)/np.sum(n)
ax2.hist(ztot,bins=bins,histtype='step',linewidth=3,color='k')
ax2.axvline(x=np.median(ztot),linestyle='dashed',linewidth=3,color='red')
ax2.hist(ztot[pany_z>p_any_cut],bins=bins,histtype='step',linewidth=3,color='k',alpha=0.2)
ax2.hist(specz,bins=bins,histtype='step',linewidth=2,color='lime',label='Spec-z')
ax2.hist(specz[pany_zsp>p_any_cut],bins=bins,histtype='step',linewidth=2,color='lime',alpha=0.2)
ax2.hist(photz,bins=bins,histtype='step',linewidth=2,color='blueviolet',linestyle='dashed',label='Photo-z')
ax2.hist(photz[pany_zph>p_any_cut],bins=bins,histtype='step',linewidth=2,color='blueviolet',linestyle='dashed',alpha=0.2)
ax2.set_ylabel('#')
ax2.tick_params(which='major',direction='inout',labelsize=12,length=7)

ax2b = ax2.twinx()
ax2b.plot(bcen,cum,'k-',linewidth=1)
ax2b.set_ylabel('Fraction')
ax2b.set_yticks(ticks=[0.2,0.4,0.6,0.8,1.0])
ax2b.tick_params('y',direction='inout',labelsize=12,length=7)
ax2.legend()

hr_con0=hr[((hr+ehrp)!=1) & ((hr-ehrn)!=-1)]
hr_con0B=hr[((hr+ehrp)!=1) & ((hr-ehrn)!=-1) & (pany_z>p_any_cut)]

ehrp_con0=ehrp[((hr+ehrp)!=1) & ((hr-ehrn)!=-1)]
ehrn_con0=ehrn[((hr+ehrp)!=1) & ((hr-ehrn)!=-1)]

z_con0=ztot[((hr+ehrp)!=1) & ((hr-ehrn)!=-1)]
z_con0B=ztot[((hr+ehrp)!=1) & ((hr-ehrn)!=-1) & (pany_z>p_any_cut)]

hr_lol0=hr[(hr+ehrp)==1]
hr_lol0B=hr[((hr+ehrp)==1) & (pany_z>p_any_cut)]

z_lol0=ztot[(hr+ehrp)==1]
z_lol0B=ztot[((hr+ehrp)==1) & (pany_z>p_any_cut)]

hr_upl0=hr[(hr-ehrn)==-1]
hr_upl0B=hr[((hr-ehrn)==-1) & (pany_z>p_any_cut)]

z_upl0=ztot[(hr-ehrn)==-1]
z_upl0B=ztot[((hr-ehrn)==-1) & (pany_z>p_any_cut)]
#what=what[zph_ga!=-99]
what=what[zflag==1]
print('Based on SED fits from Chung+14, we have',len(what[what=='A']),'AGN,',len(what[what=='G']),'GAL, and',len(what[what=='S']),'STARS. (The stars are the ones with spec-z.)')

bins=np.linspace(-1,1,30)
ax3 = plt.subplot(gs[1:3, 2:3])
ax3.hist(hr_con0,bins=bins,color='black',histtype='step',orientation='horizontal',linewidth=3)
ax3.hist(hr_con0B,bins=bins,color='black',histtype='step',orientation='horizontal',linewidth=3,alpha=0.2)
ax3.hist(hr_upl0,bins=bins,color='red',histtype='step',orientation='horizontal',linewidth=2,label='UL')
ax3.hist(hr_upl0B,bins=bins,color='red',histtype='step',orientation='horizontal',linewidth=2,alpha=0.2)
ax3.hist(hr_lol0,bins=bins,color='blue',histtype='step',orientation='horizontal',linewidth=2,label='LL')
ax3.hist(hr_lol0B,bins=bins,color='blue',histtype='step',orientation='horizontal',linewidth=2,alpha=0.2)
ax3.axhline(y=np.median(hr_con0),color='r',linestyle='dashed',linewidth=3)
print('Median HR:',np.median(hr),'Median z:',np.median(ztot),'Number of points in the plot:',len(hr),len(ztot))
ax3.legend()
ax3.set_xlabel('#')
ax3.tick_params('major',direction='inout',labelsize=12,length=7)
ax3.yaxis.set_visible(False)

z_mo=np.arange(0,5.5,0.5)
hr_mo0=[]
ax1 = plt.subplot(gs[1:3, 0:2])
for Gamma in [1.4,1.8,2.0]:
	if Gamma == 1.8:
		# Gamma 1.8 Cycle 18
		soft_cf=1.062E+11
		hard_cf=5.221E+10
		c='deepskyblue'
		label=r'$\Gamma = 1.8$'
	elif Gamma == 1.4:
		# Gamma 1.4 Cycle 18
		soft_cf=1.149E+11
		hard_cf=4.945E+10
		c='mediumspringgreen'
		label=r'$\Gamma = 1.4$'
	elif Gamma== 2.0:
		# Gamma 2.0 Cycle 18
		soft_cf=1.015E+11
		hard_cf=5.357E+10
		c='gold'
		label=r'$\Gamma = 2.0$'
	for NH in [1e21,1e22,1e23]:
		if NH == 1e21:
			style='dotted'
		elif NH == 1e22:
			style='dashed'
		elif NH == 1e23:
			style='solid'		
		flux=np.genfromtxt('/Users/alberto/Desktop/HR_Z_MYT/hr-z-pl_'+str(Gamma)+'-'+str(int(np.log10(NH)))+'.dat',unpack=True,usecols=0)
		soft=flux[::2]*soft_cf
		hard=flux[1::2]*hard_cf
		hr_mo=(hard-soft)/(hard+soft)
		hr_mo0.append(hr_mo)
		#ax1.plot(z_mo,hr_mo,color=c,linestyle=style,linewidth=3,label=label)

ax1.fill_between(z_mo,hr_mo0[0],hr_mo0[6],color='gold',label=r'log(NH)=21')
ax1.fill_between(z_mo,hr_mo0[1],hr_mo0[7],color='deepskyblue',label=r'log(NH)=22')
ax1.fill_between(z_mo,hr_mo0[2],hr_mo0[8],color='mediumspringgreen',label=r'log(NH)=23')
#plt.scatter(ztot[what=='A'],hr[what=='A'],marker='o',color='blue',alpha=0.6)
#plt.scatter(ztot[what=='G'],hr[what=='G'],marker='s',color='cyan',alpha=0.6)

#sns.kdeplot(z_lol0,hr_lol0,cmap="Blues", shade=True, shade_lowest=False,zorder=-1)
ax1.scatter(z_lol0,hr_lol0,marker='.',color='b',alpha=0.2)
#ax1.scatter(z_lol0B,hr_lol0B,marker='.',color='b')

#sns.kdeplot(z_upl0,hr_upl0,cmap="Reds", shade=True, shade_lowest=False,zorder=-1)
ax1.scatter(z_upl0,hr_upl0,marker='.',color='r',alpha=0.2)
#ax1.scatter(z_upl0B,hr_upl0B,marker='.',color='r')

#sns.kdeplot(z_con0,hr_con0,cmap="Greys", shade=True, shade_lowest=False,zorder=-1)
ax1.scatter(z_con0,hr_con0,marker='.',color='k',alpha=0.2)
#ax1.scatter(z_con0B,hr_con0B,marker='.',color='k')
ax1.errorbar(4.5,0.45,yerr=[[np.median(ehrn_con0)],[np.median(ehrp_con0)]],color='k',marker='s')

ax1.axis([-0.5,5.5,-1.2,1.2])
ax1.set_xlabel(r'$z$',fontsize=12)
ax1.set_ylabel(r'HR=(H-S)/(H+S)',fontsize=12)
ax1.tick_params(which='major',direction='inout',top=True,right=True,labelsize=12,length=7)
handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
ax1.legend(by_label.values(), by_label.keys(),loc='upper right')
plt.show()
#plt.savefig('/Users/alberto/Desktop/cdwfs_hrz_6.pdf',format='pdf')

###################
#plt.figure()
#plt.scatter(ebv[zph_ga!=-99],hr,marker='.',color='g')
#plt.xlabel(r'$E(B-V)$')
#plt.ylabel(r'HR=(H-S)/(H+S)')
#plt.show()

# LX - z plot
efluxfp0=efluxfp[zflag==1]
efluxfp0B=efluxfp0[pany_z>p_any_cut]
goodflux0=fluxf[zflag==1]
goodflux0B=goodflux0[pany_z>p_any_cut]
xb_id1=xb_id[zflag==1]
xb_id1B=xb_id1[pany_z>p_any_cut]
uplim=np.zeros_like(goodflux0)
uplimB=np.zeros_like(goodflux0B)
uplim[efluxfp0==0]=1
uplimB[efluxfp0B==0]=1

fluxrestframe=goodflux0*(1+ztot)**(-0.2)
fluxrestframeB=goodflux0B*(1+ztot[pany_z>p_any_cut])**(-0.2)
dl=cosmo.luminosity_distance(ztot)
dlB=cosmo.luminosity_distance(ztot[pany_z>p_any_cut])
dl2=dl.value*3.086e24
dl2B=dlB.value*3.086e24
lfull=4*3.141592*fluxrestframe*dl2**2
lfullB=4*3.141592*fluxrestframeB*dl2B**2
print(len(lfull[lfull<1e42]),'objects have L<10^42 erg/s, ',len(lfull[lfull>=1e42]),'otherwise.')
print(len(lfullB[lfullB<1e42]),' ROBUST objects have L<10^42 erg/s, ',len(lfullB[lfullB>=1e42]),'otherwise.')

zz=np.logspace(np.log10(1e-4),np.log10(5),100)
flim=1e-15
gamma=1.8
flimrestframe=flim*((1+zz)**(gamma-2))
dl=cosmo.luminosity_distance(zz)
dl2=dl.value*3.086e24 # from Mpc to cm
l=flimrestframe*4*3.141592*dl2**2

(x,y)=np.genfromtxt('/Users/alberto/Desktop/aird_lstar_un.dat',unpack=True)
(x2,y2)=np.genfromtxt('/Users/alberto/Desktop/aird_lstar_ob.dat',unpack=True)

lstar_un=10**(y)
lstar_ob=10**(y2)
z1=10**(x)-1
z2=10**(x2)-1

plt.figure(figsize=[7,7])
#plt.plot(ztot[xb_id1!='0'],lfull[xb_id1!='0'],'.',color='tomato',label='XBootes')
#plt.plot(ztot[xb_id1=='0'],lfull[xb_id1=='0'],'.',color='dodgerblue',label='CDWFS',zorder=-1)
plt.plot(ztot[uplim==0],lfull[uplim==0],'.',color='dodgerblue',zorder=-1,alpha=0.2)
plt.plot(ztot[(uplim==0) & (pany_z>p_any_cut)],lfullB[uplimB==0],'.',color='dodgerblue')
plt.errorbar(ztot[uplim==1],lfull[uplim==1],marker='*',yerr=0.3*lfull[uplim==1],uplims=np.ones_like(ztot[uplim==1]),linestyle='none',color='tomato',zorder=-1,alpha=0.2)
plt.errorbar(ztot[(uplim==1) & (pany_z>p_any_cut)],lfullB[uplimB==1],marker='*',yerr=0.3*lfullB[uplimB==1],uplims=np.ones_like(ztot[(uplim==1) & (pany_z>p_any_cut)]),linestyle='none',color='tomato')
plt.plot(zz,l,'k--')
plt.plot(z1,lstar_un,color='lime',linestyle='dashed',linewidth=3,label='A15 Unobscured')
plt.plot(z2,lstar_ob,color='yellow',linestyle='dashed',linewidth=3,label='A15 Obscured')
#plt.xscale('log')
plt.yscale('log')
plt.xlabel('Redshift',fontsize=20)
plt.ylabel(r'$L_{0.5-7}$ (erg/s)',fontsize=20)
plt.axis([0.05,5,5e39,1e46])
plt.tick_params(which='major',direction='inout',length=8,labelsize=15)
plt.xticks(ticks=[0.1,1,2,3,4,5],labels=[0.1,1,2,3,4,5])
#plt.annotate(r'$N=$'+str(len(ztot)),xy=(0.02,5e44))
plt.legend(loc='lower right')
plt.tight_layout()
plt.show()
#plt.savefig('/Users/alberto/Desktop/cdwfs_lxz.pdf',format='pdf')

# PHOTOMETRY
poserr=poserr[(imag<99.0) & (imag!=-99)]
newposerr=np.sqrt(poserr**2+0.1**2)

sep=sep[(imag<99.0) & (imag!=-99)]
pany_i=pany[(imag<99.0) & (imag!=-99)]
imag1=imag[(imag<99.0) & (imag!=-99)]

pany_k=pany[(kmag<99.0) & (kmag!=-99)]
sepk=sepk[(kmag<99.0) & (kmag!=-99)]
kmag=kmag[(kmag<99.0) & (kmag!=-99)]

pany_sp=pany[(spmag<99.0) & (spmag!=-99)]
sepsp=sepsp[(spmag<99.0) & (spmag!=-99)]
spmag=spmag[(spmag<99.0) & (spmag!=-99)]

#####
from scipy.stats import gaussian_kde
x=sep
y=newposerr
k = gaussian_kde(np.vstack([sep, newposerr]))
xi, yi = np.mgrid[x.min():x.max():x.size**0.5*1j,y.min():y.max():y.size**0.5*1j]
zi = k(np.vstack([xi.flatten(), yi.flatten()]))

#set zi to 0-1 scale
zi = (zi-zi.min())/(zi.max() - zi.min())
zi =zi.reshape(xi.shape)
#plt.imshow(zi,origin='lower')
#plt.show()

#set up plot
origin = 'lower'
levels = [0,0.1,0.25,0.5,0.68,0.95,0.975,1]

CS = plt.contour(xi, yi, zi,levels = levels,
              colors=('k',),
              linewidths=(1,),
              origin=origin)

plt.clabel(CS, fmt='%.3f', colors='b', fontsize=8)
plt.scatter(sep,newposerr,color='gray',marker='.',zorder=-1,alpha=0.1)
plt.show()
####

z,xedges,yedges=np.histogram2d(sep, newposerr,bins=40)
z = z / z.sum()
xcen=list((xedges[i]+xedges[i+1])/2. for i in range(len(xedges)-1))
ycen=list((yedges[i]+yedges[i+1])/2. for i in range(len(yedges)-1))

n = 1000
t = np.linspace(0, z.max(), n)
integral = ((z >= t[:, None, None]) * z).sum(axis=(1,2))

from scipy import interpolate
from scipy.interpolate import Rbf
f = interpolate.interp1d(integral, t)
t_contours = f(np.array([0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]))

#plt.scatter(sep,newposerr,color='gray',marker='.',zorder=-1,alpha=0.2)
#plt.imshow(z.T, origin='lower', extent=[-3,3,-3,3], cmap="gray")
#plt.contour(z.T, t_contours,extent=[0,5,0,5])
#plt.show()

print('I band good photometry:',len(imag1))
print('Ks band good photometry:',len(kmag))
print('[3.6] band good photometry:',len(spmag))
print('='*30)

bins=30
plt.figure()
plt.hist(imag1,bins=bins,histtype='step',color='limegreen',linewidth=3,label='I band')
plt.hist(imag1[pany_i>p_any_cut],bins=bins,histtype='step',color='limegreen',linewidth=3,alpha=0.2)
plt.hist(kmag,bins=bins,histtype='step',color='orange',linewidth=3,label='Ks band')
plt.hist(kmag[pany_k>p_any_cut],bins=bins,histtype='step',color='orange',linewidth=3,alpha=0.2)
plt.hist(spmag,bins=bins,histtype='step',color='firebrick',linewidth=3,label='[3.6] band')
plt.hist(spmag[pany_sp>p_any_cut],bins=bins,histtype='step',color='firebrick',linewidth=3,alpha=0.2)
plt.xlabel(r'Magnitude')
plt.legend(loc='upper left')
plt.tick_params(which='major',direction='inout',length=8,labelsize=13)
plt.tight_layout()
plt.show()
#plt.savefig('/Users/alberto/Desktop/cdwfs_opt-nir_mag.pdf',format='pdf')

plt.figure()
plt.hist(sep,bins=bins,histtype='step',color='limegreen',linewidth=3,label='I band')
plt.hist(sep[pany_i>p_any_cut],bins=bins,histtype='step',color='limegreen',linewidth=3,alpha=0.2)
plt.hist(sepk,bins=bins,histtype='step',color='orange',linewidth=3,label='Ks band')
plt.hist(sepk[pany_k>p_any_cut],bins=bins,histtype='step',color='orange',linewidth=3,alpha=0.2)
plt.hist(sepsp,bins=bins,histtype='step',color='firebrick',linewidth=3,label='[3.6] band')
plt.hist(sepsp[pany_sp>p_any_cut],bins=bins,histtype='step',color='firebrick',linewidth=3,alpha=0.2)
plt.xlabel(r'Separation [arcsec]')
plt.legend()
plt.tick_params(which='major',direction='inout',length=8,labelsize=13)
plt.tight_layout()
plt.show()
#plt.savefig('/Users/alberto/Desktop/cdwfs_opt-nir_sep.pdf',format='pdf')


#fig,ax=plt.subplots()
#sns.set_style("white")
#sns.kdeplot(sep, newposerr,cmap="Blues", shade=True, shade_lowest=False,zorder=-1)
#plt.scatter(sep,newposerr,color='gray',marker='.')
#plt.plot(np.linspace(0,5),np.linspace(0,5),'k--')
#plt.xlabel(r'Separation ["]')
#plt.ylabel(r'Pos Err ["]')
#plt.tight_layout()
#plt.show()

f,ax=plt.subplots(1,3,figsize=[11,5])
ax[0].scatter(sep,imag1,color='limegreen',marker='.',alpha=0.2)
ax[0].scatter(sep[pany_i>p_any_cut],imag1[pany_i>p_any_cut],color='limegreen',marker='.')
ax[0].set_xlabel('Separation ["]')
ax[0].set_ylabel('I-band mag')

ax[1].scatter(sepk,kmag,color='orange',marker='.',alpha=0.2)
ax[1].scatter(sepk[pany_k>p_any_cut],kmag[pany_k>p_any_cut],color='orange',marker='.')
ax[1].set_xlabel('Separation ["]')
ax[1].set_ylabel('K-band mag')

ax[2].scatter(sepsp,spmag,color='firebrick',marker='.',alpha=0.2)
ax[2].scatter(sepsp[pany_sp>p_any_cut],spmag[pany_sp>p_any_cut],color='firebrick',marker='.')
ax[2].set_xlabel('Separation ["]')
ax[2].set_ylabel('[3.6]-band mag')
#plt.tick_params(which='major',direction='inout',length=8,labelsize=13)
plt.tight_layout()
plt.show()
#plt.savefig('/Users/alberto/Desktop/cdwfs_opt-nir_sep-mag.pdf',format='pdf')

pany_chu=pany[zph_ga!=-99]
Imag_chu=imag[zph_ga!=-99]

pany_und=pany[zph_ga==-99]
Imag_und=imag[zph_ga==-99]

pany_chu_zsp=pany[zsp!=-99]
Imag_chu_zsp=imag[zsp!=-99]

pany_chu_zph=pany[(zsp==-99) & (zph_ga!=-99)]
Imag_chu_zph=imag[(zsp==-99) & (zph_ga!=-99)]

bins=np.linspace(12,max(Imag_chu[(Imag_chu<99.) & (Imag_chu!=-99.)]),30)
plt.figure()
plt.hist(Imag_chu[(Imag_chu<99.) & (Imag_chu!=-99.)],bins=bins,histtype='step',color='blue',linewidth=3,label='With redshift')
plt.hist(Imag_chu[(Imag_chu<99.) & (Imag_chu!=-99.) & (pany_chu>p_any_cut)],bins=bins,histtype='step',color='blue',linewidth=3,alpha=0.2)
plt.hist(Imag_chu_zsp[(Imag_chu_zsp<99.) & (Imag_chu_zsp!=-99.)],bins=bins,histtype='step',color='cyan',linewidth=2,label='Spec-z')
plt.hist(Imag_chu_zsp[(Imag_chu_zsp<99.) & (Imag_chu_zsp!=-99.) & (pany_chu_zsp>p_any_cut)],bins=bins,histtype='step',color='cyan',linewidth=2,alpha=0.2)

plt.hist(Imag_chu_zph[(Imag_chu_zph<99.) & (Imag_chu_zph!=-99.)],bins=bins,histtype='step',color='red',linewidth=2,label='Photo-z')
plt.hist(Imag_chu_zph[(Imag_chu_zph<99.) & (Imag_chu_zph!=-99.) & (Imag_chu_zph>p_any_cut)],bins=bins,histtype='step',color='red',linewidth=2,alpha=0.2)

plt.hist(Imag_und[(Imag_und<99.) & (Imag_und!=-99.)],bins=30,histtype='step',color='orange',linewidth=3,linestyle='dashed',label='No redshift')
plt.hist(Imag_und[(Imag_und<99.) & (Imag_und!=-99.) & (pany_und>p_any_cut)],bins=30,histtype='step',color='orange',linewidth=3,linestyle='dashed',alpha=0.2)

plt.axvline(x=22.5,color='k',label='AGES limit')
plt.xlabel('I-band magnitude',fontsize=15)
plt.ylabel('Number',fontsize=15)
plt.legend()
plt.tick_params(which='major',direction='inout',length=8,labelsize=13)
plt.tight_layout()
plt.show()
#plt.savefig('/Users/alberto/Desktop/cdwfs_opt-nir_zbreakdown.pdf',format='pdf')
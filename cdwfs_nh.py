import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
from scipy.stats import kde
import matplotlib.gridspec as gridspec
from collections import OrderedDict
from matplotlib import colors
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import seaborn as sns
import shapely.geometry as geom

cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

wd='/Users/alberto/Desktop/XBOOTES/'
date = '200113'
p_any_cut=0.31

# Open the matched master catalog (-cp version contains only one CDWFS source per line) 
# 6800 sources 
#cat=fits.open('/Users/alberto/Downloads/nway-master/cdwfs_I-Ks-3.6-cp.fits')
cat=fits.open(wd+'CDWFS_I-Ks-3.6_v'+date+'.fits')
data=cat[1].data
cols=cat[1].columns
names=cols.names
pany=data['p_any']

crf = data['CDWFS_CR_F']
ff  = data['CDWFS_FLUX_F']
crs = data['CDWFS_CR_S']
fs  = data['CDWFS_FLUX_S']
crh = data['CDWFS_CR_H']
fh  = data['CDWFS_FLUX_H']
hr0 = data['CDWFS_HR']
ehr0p = data['CDWFS_E_HR_+']
ehr0n = data['CDWFS_E_HR_-']
zspec0=data['zsp']
zph_g0=data['zph_G']
chi_g0=data['chi2_G']
zph_ga0=data['zph_G+A']
chi_ga0=data['chi2_G+A']
ebv0=data['E_B-V']
chi_s0=data['chi2_S']

# Array of 6800 srcs
z = np.full_like(zph_ga0,zph_ga0)
ztype = np.full_like(zph_ga0,0)

z[(chi_g0 < chi_ga0) & (chi_g0!=-99)] = zph_g0[(chi_g0 < chi_ga0) & (chi_g0!=-99)]
z[zspec0 != -99] = zspec0[zspec0 != -99]
ztype[zspec0 != -99] = 1
ztype[z==-99] = -99
what=np.full_like(chi_g0,'A',dtype='str')
what[(chi_g0<chi_ga0) & (chi_g0<chi_s0)]='G'
what[(chi_s0<chi_ga0) & (chi_s0<chi_g0)]='S'
what[chi_g0==-99]='U'

# Photometric, spectroscopic, no redshift
print(len(ztype[ztype==0]),len(ztype[ztype==1]),len(ztype[ztype==-99]))

nh,enhp,enhn=[],[],[]
ksout,ksout_up,ksout_lo=[],[],[]
khout,khout_up,khout_lo=[],[],[]
zz=np.arange(0,5.1,0.1)
zinterp=np.linspace(0.1,5.1,1001) # To more precisely interpolate
NH = np.logspace(21,25,201)

# Compute NH interpolating the HR-Z plane for each source and find the K factor
# ratio of fobs/fint for soft and hard bands
'''
for i in range(len(z)):
	print(i+1,len(z))
	if ztype[i] == -99:
		nh.append(-99)
		enhp.append(-99)
		enhn.append(-99)
		ksout.append(-99)
		ksout_up.append(-99)
		ksout_lo.append(-99)
		khout.append(-99)
		khout_up.append(-99)
		khout_lo.append(-99)
	else:
		
		distance,distance2,distance3=[],[],[]
		
		soft_cf = fs[i]/crs[i]
		hard_cf = fh[i]/crh[i]
			
		flux_unob=np.genfromtxt('/Users/alberto/Desktop/HR_Z_MYT/grid_hr-z_gamma1.8/hr-z-pl_unobsc.dat',unpack=True,usecols=0)
		soft_unob=flux_unob[::2]/soft_cf
		hard_unob=flux_unob[1::2]/hard_cf
		
		for j in range(len(NH)):
			flux = np.genfromtxt('/Users/alberto/Desktop/HR_Z_MYT/grid_hr-z_gamma1.8/hr-z-pl_'+str(round(np.log10(NH[j]),2))+'.dat',unpack=True,usecols=0)
			
			soft=flux[::2]/soft_cf
			hard=flux[1::2]/hard_cf
			
			hr=(hard-soft)/(hard+soft)
			coords = np.array([zz,hr]).T
			
			line = geom.LineString(coords)
			point = geom.Point(z[i], hr0[i])
			distance.append(point.distance(line))
			
			point2 = geom.Point(z[i], hr0[i]+ehr0p[i])
			point3 = geom.Point(z[i], hr0[i]-ehr0n[i])
			distance2.append(point2.distance(line))
			distance3.append(point3.distance(line))
			
		
		distance=np.array(distance)
		distance2=np.array(distance2)
		distance3=np.array(distance3)
		aux = np.log10(NH[distance==np.min(distance)])
		nh.append(round(aux[0],2))
		aux1 = np.log10(NH[distance2==np.min(distance2)]) # Upper limit on NH
		aux2 = np.log10(NH[distance3==np.min(distance3)]) # Lower limit on NH
		enhp.append(round(aux1[0]-aux[0],2))
		enhn.append(round(aux[0]-aux2[0],2))

		# Once I know the NH, get the correction factors
		# IF NH = 1e21, fobs = fint so k = 1
		if round(aux[0],2) == 21.0:
			
			ksout.append(1.0)
			khout.append(1.0)
		
			ksout_up.append(1.0)
			khout_up.append(1.0)
		
			ksout_lo.append(1.0)
			khout_lo.append(1.0)
			
		else:
		
			fluxNH = np.genfromtxt('/Users/alberto/Desktop/HR_Z_MYT/grid_hr-z_gamma1.8/hr-z-pl_'+str(round(aux[0],2))+'.dat',unpack=True,usecols=0)
			softNH=fluxNH[::2]/soft_cf
			hardNH=fluxNH[1::2]/hard_cf
		
			ks = softNH/soft_unob
			kh = hardNH/hard_unob
		
			ksinterp = np.interp(zinterp,zz,ks)
			khinterp = np.interp(zinterp,zz,kh)
		
			ksout.append(ksinterp[abs(zinterp-z[i])==np.min(abs(zinterp-z[i]))][0])
			khout.append(khinterp[abs(zinterp-z[i])==np.min(abs(zinterp-z[i]))][0])
		
			fluxNH1 = np.genfromtxt('/Users/alberto/Desktop/HR_Z_MYT/grid_hr-z_gamma1.8/hr-z-pl_'+str(round(aux1[0],2))+'.dat',unpack=True,usecols=0)
			softNH1=fluxNH1[::2]/soft_cf
			hardNH1=fluxNH1[1::2]/hard_cf
			
			ks1 = softNH1/soft_unob
			kh1 = hardNH1/hard_unob
		
			ksinterp1 = np.interp(zinterp,zz,ks1)
			khinterp1 = np.interp(zinterp,zz,kh1)
		
			ksout_up.append(ksinterp1[abs(zinterp-z[i])==np.min(abs(zinterp-z[i]))][0])
			khout_up.append(khinterp1[abs(zinterp-z[i])==np.min(abs(zinterp-z[i]))][0])
		
			fluxNH2 = np.genfromtxt('/Users/alberto/Desktop/HR_Z_MYT/grid_hr-z_gamma1.8/hr-z-pl_'+str(round(aux2[0],2))+'.dat',unpack=True,usecols=0)
			softNH2=fluxNH2[::2]/soft_cf
			hardNH2=fluxNH2[1::2]/hard_cf
			
			ks2 = softNH2/soft_unob
			kh2 = hardNH2/hard_unob
		
			ksinterp2 = np.interp(zinterp,zz,ks2)
			khinterp2 = np.interp(zinterp,zz,kh2)
		
			ksout_lo.append(ksinterp2[abs(zinterp-z[i])==np.min(abs(zinterp-z[i]))][0])
			khout_lo.append(khinterp2[abs(zinterp-z[i])==np.min(abs(zinterp-z[i]))][0])
		
w=open(wd+'cdwfs_nh.dat','w')
w.write('NH \t E_NH_+ \t E_NH_- \t Ksoft \t MAX_Ksoft \t MIN_Ksoft \t Khard \t MAX_Khard \t MIN_Khard \n')
for j in range(len(nh)):
	w.write(str(nh[j])+' \t '+str(enhp[j])+' \t '+str(enhn[j])+' \t '+str(ksout[j])+' \t '+str(ksout_up[j])+' \t '+str(ksout_lo[j])+' \t '+str(khout[j])+' \t '+str(khout_up[j])+' \t '+str(khout_lo[j])+'\n')
w.close()
'''

nh, ks, kh = np.genfromtxt(wd+'cdwfs_nh.dat',unpack=True,usecols=[0,3,6],skip_header=1)
nh=np.array(nh)
binsnh=np.linspace(21,25,9)

# Cut for sources with secure counterpart
z = z[pany>p_any_cut]
ztype = ztype[pany>p_any_cut]
fh = fh[pany>p_any_cut]
fs = fs[pany>p_any_cut]
nh = nh[pany>p_any_cut]
ks = ks[pany>p_any_cut]
kh = kh[pany>p_any_cut]
hr = hr0[pany>p_any_cut]
what = what[pany>p_any_cut]

# Cut for sources with redshift and NH
myz = z[ztype!=-99]
myfh = fh[ztype!=-99]
myfs = fs[ztype!=-99]
mynh = nh[nh!=-99]
myks = ks[nh!=-99]
mykh = kh[nh!=-99]
myhr = hr[nh!=-99]
mywhat = what[ztype!=-99]

print(np.min(mynh),np.max(mynh))
print(np.min(myks),np.max(myks))
print(np.min(mykh),np.max(mykh))

GAMMA = 1.8
fh_int = myfh/mykh
fs_int = myfs/myks

fluxrestframe=myfh*(1+myz)**(GAMMA-2)
fluxrestframe_int=fh_int*(1+myz)**(GAMMA-2)
dl=cosmo.luminosity_distance(myz)
dl2=dl.value*3.086e24 # from Mpc to cm
lfull=4*3.141592*fluxrestframe*dl2**2
lfull_int=4*3.141592*fluxrestframe_int*dl2**2

lbins = np.logspace(41,46,12)
zbins = np.logspace(-1,np.log10(5),15)

obsc = mynh[mynh>=22]
unob = mynh[mynh<22]

z_obsc = myz[mynh>=22]
z_unob = myz[mynh<22]

l_obsc = lfull_int[mynh>=22]
l_obsc_SM = lfull_int[myhr>-0.2] # As done by Stefano
l_obsc_phot = lfull_int[mywhat=='G']

tot,be=np.histogram(lfull_int,bins=lbins)
tot_phot,be=np.histogram(lfull_int[(mywhat!='U') & (mywhat!='S')],bins=lbins)

obs,be=np.histogram(l_obsc,bins=lbins)
obs_SM,be=np.histogram(l_obsc_SM,bins=lbins)
obs_phot,be0=np.histogram(l_obsc_phot,bins=lbins)

bce = list((be[i+1]+be[i])/2. for i in range(len(be)-1))

obsc_fraction = obs/tot
obsc_fraction_SM = obs_SM/tot
obsc_fraction_phot = obs_phot/tot_phot

plt.figure()
plt.plot(bce, obsc_fraction, 'c*',label='NH-HR')
plt.plot(bce, obsc_fraction_SM, 'gD',label='HR cut')
plt.plot(bce, obsc_fraction_phot, 'rs',label='Photometric classification')
#plt.hist(l_obsc,bins=lbins,alpha=0.6,label='OBSCURED')
#plt.hist(l_unob,bins=lbins,alpha=0.6,label='UNOBSCURED')
plt.xscale('log')
plt.legend()
plt.axis([1e42,1e45,0,1])
plt.show()

sys.exit()

zz=np.logspace(np.log10(1e-4),np.log10(5),100)
flim=3e-15
flimrestframe=flim*((1+zz)**(GAMMA-2))
dl=cosmo.luminosity_distance(zz)
dl2=dl.value*3.086e24 # from Mpc to cm
l=flimrestframe*4*3.141592*dl2**2

(x,y)=np.genfromtxt('/Users/alberto/Desktop/XBOOTES/aird_lstar_un.dat',unpack=True)
(x2,y2)=np.genfromtxt('/Users/alberto/Desktop/XBOOTES/aird_lstar_ob.dat',unpack=True)

lstar_un=10**(y)
lstar_ob=10**(y2)
z1=10**(x)-1
z2=10**(x2)-1

plt.figure(figsize=[7,7])

#plt.plot(ztot[new==0],lfull[new==0],'.',color='tomato',zorder=-1,alpha=0.5)
#plt.plot(ztot[new==1],lfull[new==1],'.',color='dodgerblue',alpha=0.5)
#plt.plot(myz,lfull,'.',color='dodgerblue',alpha=0.5)
plt.plot(myz,lfull_int,'.',color='tomato',alpha=0.5)

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

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
from collections import OrderedDict
from matplotlib import colors
from astropy.cosmology import FlatLambdaCDM
import scipy.optimize
#import seaborn as sns
from scipy.stats import kde
import math 
from scipy.special import comb
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from matplotlib.ticker import MaxNLocator
from astropy.visualization import make_lupton_rgb
from astropy.convolution import Gaussian2DKernel
from scipy.signal import convolve as scipy_convolve
from astropy.convolution import convolve
from astropy.wcs import WCS

cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*xx*3600)**2 +((pointa[1]-pointb[1])*3600)**2)
    
def integrand(x,a):
    return x**(-a+1)

def gauss(x,mu,sigma):
	g=np.exp(-(x-mu)**2/(2*sigma**2))
	return g


def smoothstep(x, x_min, x_max, N):
    N=int(N)
    
    x = np.clip((x - x_min) / (x_max - x_min), 0, 1)
	
    result = 0
    for n in range(0, N + 1):
         result += comb(N + n, n) * comb(2 * N + 1, N - n) * (-x) ** n

    result *= x ** (N + 1)

    return result

def func(x,x0,p):
	return x0 + (x-x0)*p
	
def sigmoid(x,m,b):
	return 2/(1+exp(-b*(x-m))) - 1

wd='/Users/alberto/Desktop/XBOOTES/'

cutf,cuts,cuth=10**(-4.6),10**(-4.4),10**(-4.2)

cat=fits.open(wd+'new_mosaics_detection/cdwfs_merged_cat1_exp-psf.fits')
dat=cat[1].data
ra=dat['RA']
probf=dat['PROB_F']
efluxfp=dat['E_FLUX_F_+']

probs=dat['PROB_S']
efluxsp=dat['E_FLUX_S_+']

probh=dat['PROB_H']
efluxhp=dat['E_FLUX_H_+']

print(len(ra[(probf<cutf) & (efluxfp==0.0)]))
print(len(ra[(probs<cuts) & (efluxsp==0.0)]))
print(len(ra[(probh<cuth) & (efluxhp==0.0)]))
sys.exit()

#### WRITE OUT THE XBOOTES SOURCES WE MISS AND THE ONES WITH MULTIPLE CDWFS COUNTERPART
'''
kcat=fits.open(wd+'xbootes_kenter+05.fits')
ra_k=kcat[1].data['RAJ2000']
dec_k=kcat[1].data['DEJ2000']
name_k=kcat[1].data['CXOXB']
kcat.close()

cat=fits.open(wd+'new_mosaics_detection/cdwfs_merged_cat1_exp-psf.fits')
xb_count = cat[1].data['XB_ID']
cat.close()

xb_count0 = xb_count[xb_count != '0']
print(len(xb_count0))
xb_unique = np.unique(xb_count0)
print(len(xb_unique))
sys.exit()

w=open(wd+'xbootes_missing.dat','w')
w.write('CXOXB \t RA \t DEC\n')
for i in range(len(ra_k)):	
	if len(xb_count[xb_count == name_k[i]]) == 0:
		w.write(str(name_k[i])+' \t '+str(ra_k[i])+' \t '+str(dec_k[i])+'\n')
w.close()
sys.exit()
'''

#### HOW MANY XBOOTES SOURCES WE GET WITH A GIVEN PROBABILITY CUT?
'''
cat=fits.open(wd+'new_mosaics_detection/cdwfs_merged_cat0_exp-psf.fits')
raf=cat[1].data['RA_F']
decf=cat[1].data['DEC_F']
probf=cat[1].data['PROB_F']
r90f=cat[1].data['R90_F']
totf=cat[1].data['TOT_F']
bkgf=cat[1].data['BKG_F']
netf=cat[1].data['NET_F']
e_netf_up=cat[1].data['E_NET_F_+']
e_netf_lo=cat[1].data['E_NET_F_-']
expf=cat[1].data['EXP_F']
crf=cat[1].data['CR_F']
e_crf_up=cat[1].data['E_CR_F_+']
e_crf_lo=cat[1].data['E_CR_F_-']
fluxf=cat[1].data['FLUX_F']
e_fluxf_up=cat[1].data['E_FLUX_F_+']
e_fluxf_lo=cat[1].data['E_FLUX_F_-']

ras=cat[1].data['RA_S']
decs=cat[1].data['DEC_S']
probs=cat[1].data['PROB_S']
r90s=cat[1].data['R90_S']
tots=cat[1].data['TOT_S']
bkgs=cat[1].data['BKG_S']
nets=cat[1].data['NET_S']
e_nets_up=cat[1].data['E_NET_S_+']
e_nets_lo=cat[1].data['E_NET_S_-']
exps=cat[1].data['EXP_S']
crs=cat[1].data['CR_S']
e_crs_up=cat[1].data['E_CR_S_+']
e_crs_lo=cat[1].data['E_CR_S_-']
fluxs=cat[1].data['FLUX_S']
e_fluxs_up=cat[1].data['E_FLUX_S_+']
e_fluxs_lo=cat[1].data['E_FLUX_S_-']

rah=cat[1].data['RA_H']
dech=cat[1].data['DEC_H']
probh=cat[1].data['PROB_H']
r90h=cat[1].data['R90_H']
toth=cat[1].data['TOT_H']
bkgh=cat[1].data['BKG_H']
neth=cat[1].data['NET_H']
e_neth_up=cat[1].data['E_NET_H_+']
e_neth_lo=cat[1].data['E_NET_H_-']
exph=cat[1].data['EXP_H']
crh=cat[1].data['CR_H']
e_crh_up=cat[1].data['E_CR_H_+']
e_crh_lo=cat[1].data['E_CR_H_-']
fluxh=cat[1].data['FLUX_H']
e_fluxh_up=cat[1].data['E_FLUX_H_+']
e_fluxh_lo=cat[1].data['E_FLUX_H_-']

cutf,cuts,cuth=10**(-4.6),10**(-4.4),10**(-4.2)
#cutf,cuts,cuth=1,1,1

probf[raf==0.0]=9999
probs[ras==0.0]=9999
probh[rah==0.0]=9999

print(len(raf[probf<=cutf]))
print(len(ras[probs<=cuts]))
print(len(rah[probh<=cuth]))
print('we should have a total of',len(raf[(probf<=cutf) | (probs<=cuts) | (probh<=cuth)]),'sources.')

#sys.exit()

count=0
ra_u,dec_u,r90_u=[],[],[]
hr,e_hr_up,e_hr_lo=[],[],[]
out_name,unique,id,poserr=[],[],[],[]
kenter_match=0
kcat=fits.open(wd+'xbootes_kenter+05.fits')
ra_k=kcat[1].data['RAJ2000']
dec_k=kcat[1].data['DEJ2000']
name_k=kcat[1].data['CXOXB']
used_r90=0
for j in range(len(probf)):
	if (probf[j] <= cutf or probs[j] <= cuts or probh[j] <= cuth): # Source is above threshold in at least one band
		count=count+1
		id.append(count)
		unique.append(True)
		prob_array=[probf[j],probs[j],probh[j]]
		
		if min(prob_array) == probf[j]: # Choose the ra,dec based on most significant band
			ra_u.append(raf[j])
			dec_u.append(decf[j])
			r90_u.append(r90f[j])
			
		elif min(prob_array) == probs[j]:
			ra_u.append(ras[j])
			dec_u.append(decs[j])
			r90_u.append(r90s[j])
			
		else:
			ra_u.append(rah[j])
			dec_u.append(dech[j])
			r90_u.append(r90h[j])
			

		## Match with Kenter+05
		cdwfs=[ra_u[-1],dec_u[-1]]

		found=False
		ra_cont,dec_cont,d,name=[],[],[],[]
		for k in range(len(ra_k)): # Look up into Kenter catalog
			kenter=[ra_k[k],dec_k[k]]
			if distance(cdwfs,kenter) < 1.1*r90_u[-1]:
				if found==False: # This is the first match
					kenter_match=kenter_match+1
					found=True
					ra_cont.append(ra_k[k])
					dec_cont.append(dec_k[k])
					d.append(distance(cdwfs,kenter))
					name.append(name_k[k])
				elif found==True:
					ra_cont.append(ra_k[k])
					dec_cont.append(dec_k[k])
					d.append(distance(cdwfs,kenter))
					name.append(name_k[k])
					print('Finding a double countepart for '+str(ra_u[-1])+', '+str(dec_u[-1])+'')
		for l in range(len(ra_cont)): # retain the closest counterpart only
			if min(d)==d[l]:
				out_name.append(name[l])
		if found==False: # No kenter counterpart was found
			out_name.append(0)

	else:
		unique.append(False)
kcat.close()

print(kenter_match, 'Kenter')
sys.exit()
'''

#### CREATE REGION FILE
'''
# Open file
band='broad'
(ra,dec)=np.genfromtxt(wd+'poiss_rand_'+band+'_sim_indep_22-Nov-19_filtered_new.dat',unpack=True,skip_header=1, usecols=[1,2])
w=open(wd+'poiss_rand_'+band+'_sim_indep_22-Nov-19_filtered_new.reg','w')
for k in range(len(ra)):
	w.write('circle('+str(ra[k])+'d,'+str(dec[k])+'d,2\") #color=magenta \n')
w.close()
sys.exit()
'''

#### ANALYSIS FOR ADI FOORD

cat=fits.open(wd+'CDWFS_I-Ks-3.6_v120119.fits')
data=cat[1].data
r90=data['CDWFS_R90_F']
zspec=data['zsp']
zphot=data['zph_G+A']
totf=data['CDWFS_NET_F']

z = zspec.copy()
flag = np.full(len(z),'sp')

z[zspec==-99.]=zphot[zspec==-99.]

flag[zspec==-99.] = 'ph'

flag = flag[z!=-99.]
theta = 10.*np.sqrt((r90-1.)/10.)
theta = theta[z!=-99]
totf = totf[z!=-99.]
z = z[z!=-99.]

print(len(z[z>=2]), 'sources z>2')
print(len(z[z>=3]), 'sources z>3')
print(len(totf[totf>30.]), 'sources NET_F>30')
sys.exit()

cut = 5.

oaa5 = theta[theta<cut]
z_oaa5 = z[theta<cut]
zsp_oaa5 = z[(theta<cut) & (flag=='sp')]
zph_oaa5 = z[(theta<cut) & (flag=='ph')]
print(len(z_oaa5), len(zsp_oaa5), len(zph_oaa5))
print(len(zsp_oaa5[zsp_oaa5<1]), len(zsp_oaa5[zsp_oaa5<2]), len(zsp_oaa5[zsp_oaa5<3]))
print(len(z_oaa5[z_oaa5<1]), len(z_oaa5[z_oaa5<2]), len(z_oaa5[z_oaa5<3]))
bins = np.linspace(0,5,21)

plt.figure()
plt.hist(z_oaa5, bins = bins, histtype='step', linewidth = 3, color='k')
plt.hist(zsp_oaa5, bins=bins, histtype='step', linewidth = 1, color='b', label='Spec')
plt.hist(zph_oaa5, bins=bins, histtype='step', linewidth = 1, color='r', label='Phot')
plt.legend()
plt.xlabel('Redshift', fontsize=12)
plt.tight_layout()
plt.show()
sys.exit()

#### CREATE PDF FOR GAMMA AND CHECK IT
'''
(cf_f,cf_s,cf_h)=np.genfromtxt(wd+'cdwfs_ecf_flux-to-cr_CY18.dat',unpack=True,skip_header=1,usecols=[1,2,3])
cf_f=cf_f*1e10
cf_s=cf_s*1e10
cf_h=cf_h*1e10

# Define the Gamma PDF
xvals = np.arange(0,11,1)

mu=5
sigma=2
p=[]
for ii in range(len(xvals)):
	p.append(gauss(xvals[ii],mu,sigma))
p=np.array(p)
new=p/np.sum(p)

random_index=np.random.choice(xvals, p=new)
print((random_index+9)/10.)
print(cf_f[random_index],cf_s[random_index],cf_h[random_index])
sys.exit()

plt.figure()
plt.hist(random_index, bins = 20, histtype = 'step', density= True)
plt.plot(xvals,new,'k-')
plt.axvline(mu)
plt.axvline(mu-sigma)
plt.axvline(mu+sigma)
plt.show()

sys.exit()
'''

#### COMPARE R50 AND R90 SENSITIVITY, OR NORMAL AND GEORGAKAKIS SENSITIVITY
'''
x0,y0 = np.genfromtxt(wd+'cdwfs_hard_sens_r90.dat',unpack=True)
x1,y1 = np.genfromtxt(wd+'cdwfs_hard_sens_georgakakis_r90.dat',unpack=True)

plt.figure()
plt.plot(x0,y0,'k-',label='R50')
plt.plot(x1,y1,'r-',label='R90')
plt.xscale('log')
plt.yscale('log')
plt.show()
sys.exit()
'''

#### Plot dNdS AND FIT THE NORMALIZATION?
'''
band = 'hard'
cut = 1e-4

# Real data
cat=fits.open(wd+'new_mosaics_detection/cdwfs_merged_cat1_exp-psf.fits')
data=cat[1].data
if band == 'broad':
	ra =  data['RA']
	dec = data['DEC']
	exp = data['EXP_F']
	tot = data['TOT_F']
	bkg = data['BKG_F']
	cr = data['CR_F']
	flux1 = data['FLUX_F']
	eflux1 = data['E_FLUX_F_+']
	prob = data['PROB_F']
elif band == 'soft':
	ra =  data['RA']
	dec = data['DEC']
	exp = data['EXP_S']
	tot = data['TOT_S']
	bkg = data['BKG_S']
	cr = data['CR_S']
	flux1 = data['FLUX_S']
	eflux1 = data['E_FLUX_S_+']
	prob = data['PROB_S']
else:
	ra =  data['RA']
	dec = data['DEC']
	exp = data['EXP_H']
	tot = data['TOT_H']
	bkg = data['BKG_H']
	cr = data['CR_H']
	flux1 = data['FLUX_H']
	eflux1 = data['E_FLUX_H_+']
	prob = data['PROB_H']
	probf = data['PROB_F']
	probs = data['PROB_S']
cat.close()

# Apply probability cut
ra = ra[prob <= cut]
dec = dec[prob <= cut]
exp = exp[prob <= cut]
tot = tot[prob <= cut]
bkg = bkg[prob <= cut]
cr = cr[prob <= cut]
eflux1 = eflux1[prob <= cut]
flux1 = flux1[prob <= cut]
prob1 = prob[prob <= cut]

# Exclude upper limits
ra = ra[eflux1 != 0]
dec = dec[eflux1 != 0]
exp = exp[eflux1 != 0]
tot = tot[eflux1 != 0]
bkg = bkg[eflux1 != 0]
cr = cr[eflux1 != 0]
flux0 = flux1[eflux1 != 0]
prob0 = prob1[eflux1 != 0]
eflux0 = eflux1[eflux1 != 0]

print(len(flux0), 'sources before SNR cut')

# Apply a SNR cut

snr_cut = 2.5
net = tot-bkg
snr = net/bkg

ra = ra[snr > snr_cut]
dec = dec[snr > snr_cut]
exp = exp[snr > snr_cut]
tot = tot[snr > snr_cut]
bkg = bkg[snr > snr_cut]
cr = cr[snr > snr_cut]
flux0 = flux0[snr > snr_cut]
prob0 = prob0[snr > snr_cut]
eflux0 = eflux0[snr > snr_cut]

print(len(flux0), 'sources after SNR cut')

# Apply an exposure cut
#expo_cut = 9e4

#ra = ra[exp < expo_cut]
#dec = dec[exp < expo_cut]
#tot = tot[exp < expo_cut]
#bkg = bkg[exp < expo_cut]
#cr = cr[exp < expo_cut]
#flux0 = flux0[exp < expo_cut]
#prob0 = prob0[exp < expo_cut]
#probf = probf[exp < expo_cut]
#probs = probs[exp < expo_cut]
#eflux0 = eflux0[exp < expo_cut]
#exp = exp[exp < expo_cut]


print(np.min(flux0),np.max(flux0))

if band == 'soft':
	bins = np.logspace(-16,-12,51)
else:
	bins = np.logspace(-15,-12,41)
x = list((bins[i+1]+bins[i])/2. for i in range(len(bins)-1))
x = np.array(x)
x0,y0 = np.genfromtxt(wd+'cdwfs_'+band+'_sens_r90.dat',unpack=True)
sens = np.interp(x,x0,y0)

# Check the interpolation
#plt.figure()
#plt.plot(x0,y0,'ko')
#plt.plot(x,sens,'r+')
#plt.xscale('log')
#plt.show()
#sys.exit()


a,b = np.histogram(flux0, bins = bins)

c = a/sens

y = list(reversed(np.cumsum(list(reversed(c)))))


#x,y,ey = np.genfromtxt(wd+'cdwfs_lognlogs_'+band+'_cutexp8e4.dat',unpack=True)
(civx_lo,civy_lo)=np.genfromtxt(wd+'civano_lognlogs_'+band+'_lo.txt',unpack=True)
(civx_hi,civy_hi)=np.genfromtxt(wd+'civano_lognlogs_'+band+'_hi.txt',unpack=True)

civy_hi2 = np.interp(civx_lo, civx_hi, civy_hi)

if band == 'hard':
	x = x/0.75 # convert from 2-7 to 2-10 using Gamma = 1.8
	#x = x/6.887E-01 # convert from 2-7 to 2-10 using Gamma = 1.4
	
	#x = x*1.615 # convert from 0.5-2 to 2-10 using Gamma = 1.8
	#x = x*2.955 # convert from 0.5-2 to 2-10 using Gamma = 1.4
	
y = y*(x/1e-14)**1.5
#ey = ey*(x/1e-14)**1.5

plt.figure()
#plt.errorbar(x,y,yerr=ey,marker='o',color='k',linestyle='None',label='CDWFS')
plt.plot(x,y,marker='o',color='k',linestyle='None',label='CDWFS')
plt.fill_between(civx_lo,civy_lo,civy_hi2,color='red',alpha = 0.3,label='Civano+16')
plt.xscale('log')
plt.yscale('log')
if band == 'soft':
	plt.xlabel('0.5-2 keV Flux')
elif band == 'hard':
	plt.xlabel('2-10 keV Flux')

plt.ylabel('N(>S)*S^1.5')
plt.axis([3e-17,1e-12,1,500])
plt.legend()
plt.show()

sys.exit()


# Function to compute dN/dS given parameters
def dnds(fx,k):

	fref = 1e-14
	
	k1 = k*(fb/fref)**(b1-b2)
	
	fx = 10**fx
	
	if type(fx) == float:
		if fx <= fb:
			return k*(fx/fref)**b1
		else:
			return k1*(fx/fref)**b2
	elif (type(fx) == list)	or (type(fx) == np.ndarray):
		if type(fx) == list:
			fx = np.array(fx)
		aux = k*(fx/fref)**b1
		aux[fx > fb] = k1*(fx[fx > fb]/fref)**b2
		return np.log10(aux)

x,y,ey = np.genfromtxt(wd+'cdwfs_dnds_'+band+'.dat',unpack=True)
y = y[x<1e-12]
ey = ey[x<1e-12]
x=x[x<1e-12]

x0 = np.log10(x)
y0 = np.log10(y)
ey0 = ey/(np.log(10)*y)

b1,b2,fb = np.genfromtxt(wd+'cdwfs_dnds_bfit-pars_'+band+'.dat',unpack=True,usecols=[0,2,4],skip_header=1)

popt,pcov = curve_fit(dnds, x0, y0, sigma = ey0)
print(popt,np.sqrt(pcov))

fit = dnds(x0,popt[0])

plt.figure()
plt.errorbar(x0,y0,yerr=ey0,marker='o',color='red',linestyle='None')
plt.plot(x0,fit,'b-')
#plt.xscale('log')
#plt.yscale('log')
#plt.axis([1e-16,1e-13,1e11,1e21])
plt.show()
sys.exit()
'''

#### SENSITIVITY, COMPARISON WITH XBOOTES
'''
band = ['broad','soft','hard']

f,ax=plt.subplots(len(band),1,sharey=True,sharex=True,figsize=[4,12])
for j in range(len(band)):
	
	#if j == 0:
	#	fl,ar=np.genfromtxt(wd+'kenter05_sens_05-7keV.dat',unpack=True)
	#	x=np.log10(fl)
	#	ax[j].plot(x,ar,color='k',linestyle='--',linewidth=2)
		
	fl2,ar2=np.genfromtxt(wd+'cdwfs_'+band[j]+'_sens_georgakakis_r90.dat',unpack=True)

	x2=np.log10(fl2)
	
	ax[j].plot(x2,ar2,color='k',linestyle='-',linewidth=3)

	ax[j].set_ylabel(r'Area (deg$^2$)',fontsize=15)
	
	ax[j].axis([-16.5,-13,1e-3,10])
	ax[j].annotate(band[j].capitalize(),xy=(-16,8),fontsize=15)
	
	nbins = len(ax[j].get_yticklabels())  
	ax[j].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper')) 
		
ax[2].tick_params(axis='x',labelsize=13)
ax[2].set_xlabel(r'Log(Flux/(erg cm$^{-2}$ s$^{-1}$))',fontsize=15)
plt.tight_layout()
plt.subplots_adjust(hspace=0)
plt.show()
sys.exit()
'''

#### EXPOSURE-CORRECTED RGB IMAGE OF MOSAIC
'''
r = fits.open('/Volumes/Dupree/CDWFS/0.5-2.0_flux.img')
g = fits.open('/Volumes/Dupree/CDWFS/2.0-4.5_flux.img')
b = fits.open('/Volumes/Dupree/CDWFS/4.5-7.0_flux.img')

smig_wcs = WCS(r[0].header)

image_r = r[0].data
image_g = g[0].data
image_b = b[0].data
r.close()
g.close()
b.close()

# We smooth with a Gaussian kernel with x_stddev=1 (and y_stddev=1)
# It is a 9x9 array
kernel = Gaussian2DKernel(x_stddev=5)
conv_r = convolve(image_r, kernel)
conv_g = convolve(image_g, kernel)
conv_b = convolve(image_b, kernel)
image = make_lupton_rgb(conv_r, conv_g, conv_b, Q=5, stretch=2.5e-8)

#image = make_lupton_rgb(image_r, image_g, image_b, Q=10, stretch=5e-8)

plt.subplot(projection=smig_wcs)
plt.grid(color='white', ls='solid')
plt.imshow(image,origin = 'lower')
plt.show()

sys.exit()
'''

#### MASK THE SOURCES CLOSE TO THE BORDERS OF THE MOSAIC
'''
band = 'soft'
# take input sources in full band
(f_s,ra_s,dec_s)=np.genfromtxt(wd+'poiss_rand_'+band+'_filtered.dat',skip_header=1,unpack=True,usecols=[0,1,2])

x0,x1 = 216.35,219.5
y0,y1 = 32.475,35.725

x00,y00=218.06,32.30
x11,y11=219.31,33.48

def eqy(x):
	return (y11-y00)/(x11-x00)*(x-x00)+y00

mask = (ra_s >= x0) & (ra_s <= x1) & (dec_s >= y0) & (dec_s <= y1) & (dec_s >= eqy(ra_s))
filt_x = ra_s[~mask]
filt_y = dec_s[~mask]

print(len(ra_s),len(filt_x))

plt.figure(figsize=[4,5])
#plt.scatter(ra_s,dec_s)
plt.scatter(filt_x,filt_y)
plt.show()
sys.exit()
'''

#### CREATE REGION FILES TO COMPARE WAVDETECT AND INPUT OF SIMULATIONS (SPURIOUS AND SOURCES
#### NOT PICKED UP BY WAVDETECT)
'''
band = 'soft'
cut=1e-4

ra,dec=np.genfromtxt(wd+'poiss_rand_'+band+'_filtered_new.dat',unpack=True,skip_header=1,usecols=[1,2])
w=open(wd+'poiss_rand_'+band+'_filtered_new.reg','w')
for i in range(len(ra)):
	w.write('circle('+str(ra[i])+'d, '+str(dec[i])+'d, 1\") # color=magenta\n')
w.close()

cat0=fits.open(wd+'sim_indep/0cdwfs_'+band+'_sim_cat1.fits')
wav_ra=cat0[1].data['RA']
wav_dec=cat0[1].data['DEC']
wav_prob=cat0[1].data['PROB']
cat0.close()

# The green sources are detected by wavdetect, so the ONLY GREEN are the spurious 
wav_ra=wav_ra[wav_prob < cut]
wav_dec=wav_dec[wav_prob < cut]
w=open(wd+'sim_indep/0cdwfs_'+band+'_sim_cat1.reg','w')
for i in range(len(wav_ra)):
	w.write('circle('+str(wav_ra[i])+'d, '+str(wav_dec[i])+'d, 10\") # color=green\n')
w.close()

cat0=fits.open(wd+'sim_indep/0cdwfs_'+band+'_sim_cat1_lehmer.fits')
inp_ra=cat0[1].data['RA']
inp_dec=cat0[1].data['DEC']
inp_prob=cat0[1].data['PROB']
cat0.close()

# The red sources are the ones in input, so the ONLY RED ones are the ones missing 
inp_ra=inp_ra[inp_prob < cut]
inp_dec=inp_dec[inp_prob < cut]
w=open(wd+'sim_indep/0cdwfs_'+band+'_sim_cat1_lehmer.reg','w')
for i in range(len(inp_ra)):
	w.write('circle('+str(inp_ra[i])+'d, '+str(inp_dec[i])+'d, 5\") # color=red\n')
w.close()

sys.exit()
'''

#### SENSITIVITY, COMPARE GEORGAKAKIS WITH SIMULATIONS
band = ['broad','soft','hard']
logcut=['-4.6','-4.4','-4.2']
type=['sim_indep_22-Nov-19']
f,ax=plt.subplots(1,len(band),sharey=True,sharex=True,figsize=[15,5])
for j in range(len(band)):
	for i in range(len(type)):
	
		centers01,ratio3,eratio3=np.genfromtxt(wd+'cdwfs_'+band[j]+'_sens_'+logcut[j]+'_'+type[i]+'.dat',unpack=True)
		fl2,ar2=np.genfromtxt(wd+'cdwfs_'+band[j]+'_sens_'+logcut[j]+'_geo.dat',unpack=True)
		fl,ar=np.genfromtxt(wd+'cdwfs_'+band[j]+'_sens_'+logcut[j]+'_civ.dat',unpack=True)
		
		eratio3[np.isnan(eratio3)]=ratio3[np.isnan(eratio3)]
		ratio3[np.isnan(ratio3)]=9.3
		#uplims=np.zeros(len(ratio3))
		#uplims[eratio3>=ratio3]=1
		#print(uplims)
		
		ar2_sup=ar2*2.
		ar2_inf=ar2/2.
		
		x0=np.log10(centers01)
		x=np.log10(fl)
		x2=np.log10(fl2)
				
		ax[j].plot(x2,ar2,color='k',linestyle='-',linewidth=3,label='Analytical')
		#ax[j].plot(x,ar,color='k',linestyle='--',linewidth=3,label='No Eddington bias')
		#ax[j].fill_between(x2,ar2_inf,ar2_sup,color='gray',alpha=0.5)
		
		ax[j].errorbar(x0,ratio3,yerr=eratio3,color='C'+str(i+2),marker='o',linestyle='-',linewidth=2,label='Simulation')

		ax[j].set_xlabel(r'Log(Flux/(erg cm$^{-2}$ s$^{-1}$))',fontsize=15)
		#ax[j].set_yscale('log')
		ax[j].axis([-16.5,-13,1e-3,10])
		#ax[j].axhline(y=0.08,color='gray',linestyle='dashed')
		#ax[j].annotate('ACIS-I FoV', color='gray',xy=(-13.7,0.1))
		#ax[j].annotate(band[j].capitalize(),xy=(-13.4,1.5e-3))
		
		nbins = len(ax[j].get_xticklabels()) # added 
		ax[j].xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper')) # added 
		
ax[0].set_ylabel(r'Area (deg$^2$)',fontsize=15)
ax[0].tick_params(axis='y',labelsize=13)
#ax[0].set_yticks([1e-3,1e-2,0.1,1,10])
#ax[0].set_yticklabels([r'$10^{-3}$',1e-2,0.1,1,10])
plt.legend()
plt.tight_layout()
plt.subplots_adjust(wspace=0)
plt.show()
#plt.savefig(wd+'sens.pdf',format='pdf')
sys.exit()


#### GAMMA - ENERGY CONVERSION FACTORS ACROSS CHANDRA CYCLES
'''
#gamma=np.genfromtxt(wd+'poiss_rand_lehmer_filtered.dat',unpack=True,skip_header=1,usecols=3)
#plt.figure()
#plt.hist(gamma,bins=20)
#plt.xlabel('Photon Index',fontsize=15)
#plt.tick_params(which='both',direction='inout',length=8,labelsize=15)
#plt.tight_layout()
#plt.show()
#sys.exit()

cycles=['3','4','7','8','9','10','12','14','16','17','18']

f,ax=plt.subplots(1,sharey=True,figsize=(8,8))

gamma=np.arange(0.9,2.0,0.1)

for i in range(len(cycles)):

	(g,cf_f,cf_s,cf_h)=np.genfromtxt(wd+'cdwfs_ecf_flux-to-cr_CY'+cycles[i]+'.dat',unpack=True,skip_header=1, skip_footer=4)
	
	datax,datayf,datays,datayh=[],[],[],[]
	for j in range(len(gamma)):	
	
		datax.append(cycles[i])
		datayf.append(cf_f[j]*1e10)
		datays.append(cf_s[j]*1e10)
		datayh.append(cf_h[j]*1e10)

	scf = plt.scatter(datax,datayf,c=np.arange(0,11),cmap='Greens')
	scs = plt.scatter(datax,datays, c=np.arange(0,11),cmap='Reds')
	sch = plt.scatter(datax,datayh,c=np.arange(0,11),cmap='Blues')


cbarh = plt.colorbar(sch,pad=-0.085, ticks=np.arange(0,11))
cbarh.ax.set_yticklabels(np.around(gamma,1))
cbarh.ax.set_ylabel('Photon index', rotation=270, fontsize=12, labelpad = 15)
cbarh.ax.tick_params(axis='both', labelsize=12, direction = 'in', top=True, right=True, length=7)
cbars =plt.colorbar(scs, pad =-0.08)
cbars.set_ticks([])
cbarf = plt.colorbar(scf, pad =0.01)
cbarf.set_ticks([])
plt.yscale('log')
plt.xlabel('Chandra Cycle', fontsize=12)
plt.ylabel(r'Flux to CR (cm$^2$/erg)', fontsize=12., labelpad=-5)
plt.tick_params(axis='both', which='major', labelsize=12, direction = 'in', top=True, right=True, length=7)
plt.tick_params(axis='y', which='minor', labelsize=12, direction = 'in', top=True, right=True, length=4)
plt.show()
#plt.savefig('/Users/alberto/Desktop/ecfs2.pdf', format='pdf')

sys.exit()
'''

# Version 1, long horizontal plot
'''
f,ax=plt.subplots(1,sharey=True,figsize=(12,4))
ticks,new_tick_locations=[],[]
for i in range(len(cycles)):
	
	(gamma,cf_f,cf_s,cf_h)=np.genfromtxt(wd+'cdwfs_ecf_flux-to-cr_CY'+cycles[i]+'.dat',unpack=True,skip_header=1, skip_footer=4)
	
	gamma_hi = gamma[-1]
	mean_gamma = (gamma[0]+gamma_hi)/2.
	
	x=gamma+(i*(gamma_hi+0.1))
	ax.plot(x,cf_f,'g.',linestyle='solid')
	ax.plot(x,cf_s,'r.',linestyle='solid')
	ax.plot(x,cf_h,'b.',linestyle='solid')
	ticks.append(mean_gamma+i*(gamma_hi+0.1))
	new_tick_locations.append(x[0])
	new_tick_locations.append(x[-1])
	ax.axvline(x[0],linestyle='--',color='gray')
	ax.axvline(x[-1],linestyle='--',color='gray')
	if i == 0:
		ax.axvspan(0,min(gamma)+(i*(gamma_hi+0.1)),color='gray',hatch='x', alpha=0.3)
	elif i < len(cycles):
		ax.axvspan(max(gamma)+((i-1)*(gamma_hi+0.1)),min(gamma)+(i*(gamma_hi+0.1)),color='gray',hatch='x', alpha=0.3)
		if i == len(cycles)-1:
			#ax.axvspan(max(gamma)+((i+1)*(gamma_hi+1)),30.9,color='gray',hatch='x', alpha=0.3)
			ax.axvspan(21.9,22.9,color='gray',hatch='x', alpha=0.3)
		
plt.xticks(ticks,labels=cycles)
plt.text(28,22,s='F', color='g',fontweight='bold')
plt.text(28,21,s='S', color='r',fontweight='bold')
plt.text(28,20,s='H', color='b',fontweight='bold')
plt.xlabel('Chandra Cycle')
plt.ylabel(r'Flux to CR $\times 10^{10}$ cm$^2$/erg')
#plt.axis([0,29,4,23])

ax2 = ax.twiny()

def tick_function(X):
	s = []
	for i in range(len(X)):
		if i%2 == 0:
			s.append('0.9')
		else:
			s.append('1.9')
	return s

ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(tick_function(new_tick_locations))
ax2.set_xlabel(r'$\Gamma$')

plt.show()
sys.exit()
'''

#### RELIABILITY PLOTS - OLD, WAS FIXING 1% OF SPURIOUS FRACTION (DEPENDENT ON INPUT LOGNLOGS)
'''
color=['green','red','blue']
for source in ['sim_all_FINAL/','sim_all_new/']:
	i=0
	if source == 'sim_all_FINAL/':
		mark='o'
		for band in ['broad']:
			cut,sp_frac=np.genfromtxt(wd+source+'0cdwfs_'+band+'_sp-frac.dat',unpack=True)
			plt.plot(cut,sp_frac,marker=mark,color=color[i],linestyle='-',label=band)
			i=i+1
	elif source == 'sim/':
		mark='*'
		for band in ['broad','soft','hard']:
			spfr=[]
			for j in range(10):
				cut,sp_frac=np.genfromtxt(wd+source+str(j)+'cdwfs_'+band+'_sp-frac.dat',unpack=True)
				spfr.append(sp_frac)
				plt.plot(cut,sp_frac,marker=mark,color=color[i],linestyle='-',alpha=0.1)
			spfr=np.array(spfr)
			#print(spfr)
			#print('*'*15)
			av,med,emed=[],[],[]
			for ii in range(spfr.shape[1]):
				med.append(np.median(spfr[:,ii]))
				av.append(np.mean(spfr[:,ii]))
				emed.append(np.std(spfr[:,ii]))
			plt.errorbar(cut,med,yerr=emed,marker='s',color=color[i],linestyle='-',label=band)
			#plt.errorbar(cut,av,yerr=emed,marker='D',color=color[i],linestyle='-',label=band)
			i=i+1
	else:
		mark='D'
		for band in ['broad','soft','hard']:
			spfr=[]
			for j in range(10):
				cut,sp_frac=np.genfromtxt(wd+source+str(j)+'cdwfs_'+band+'_sp-frac.dat',unpack=True)
				spfr.append(sp_frac)
				plt.plot(cut,sp_frac,marker=mark,color=color[i],linestyle='--',alpha=0.1)
			spfr=np.array(spfr)
			av,med,emed=[],[],[]
			for ii in range(spfr.shape[1]):
				med.append(np.median(spfr[:,ii]))
				av.append(np.mean(spfr[:,ii]))
				emed.append(np.std(spfr[:,ii]))
			plt.errorbar(cut,med,yerr=emed,marker='d',color=color[i],linestyle='--',label=band)
			i=i+1
			
plt.axhline(y=1.0)
plt.axhline(y=3.0)
plt.xscale('log')
plt.xlabel('Prob cut',fontsize=15)
plt.ylabel('Sp fraction (%)',fontsize=15)
plt.tick_params(which='both',labelsize=13)
plt.axis([5e-7,0.2,0,2])
handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys())
plt.tight_layout()
plt.show()
#plt.savefig('/Users/alberto/Desktop/cdwfs_spfrac.pdf',format='pdf')

sys.exit()
'''

#### RANDOM STUFF ON LOGNLOGS
'''
band='broad'
def dnds(s):
    n=[]
    k,fb,b1,b2=562.20,8.1e-15,1.34,2.35 #full band Lehmer
    
    for i in range(len(s)):
        if s[i] <= fb:
            n.append(k*1e14*(s[i]/1.e-14)**(-b1))
        else:
            n.append(k*1e14*(fb/1.e-14)**(b2-b1)*(s[i]/1.e-14)**(-b2))
    return n

def dndlogs(s):
    n=[]
    k,fb,b1,b2=562.20,8.1e-15,1.34,2.35 #full band Lehmer
    
    for i in range(len(s)):
        if 10**s[i] <= fb:
            n.append(k*1e14*((10.**s[i])/1.e-14)**(-b1))
        else:
            n.append(k*1e14*(fb/1.e-14)**(b2-b1)*((10.**s[i])/1.e-14)**(-b2))
    return n



bins1=np.linspace(-16.5,-12,51)
centers1=list((bins1[i+1]+bins1[i])/2. for i in range(0,len(bins1)-1))
centers1=np.array(centers1) 

centers0=10**centers1

x=centers1
y=dndlogs(x)
y2=dnds(centers0)
plt.figure()
plt.plot(x,y,'k-')
plt.plot(x,y2,'r--')
plt.yscale('log')
plt.show()


bins1=np.linspace(-16.5,-12,51)
centers1=list((bins1[i+1]+bins1[i])/2. for i in range(0,len(bins1)-1))
centers1=np.array(centers1) 
ds1=np.array(ds1)
b1=-1.58
b2=-2.59
fb=2e-14

if band=='broad':
	ecf=1.411E-11 # gamma=1.4 for Chandra Cy18
elif band=='soft':
	ecf=8.707E-12 # gamma=1.4 for Chandra Cy18
elif band=='hard':
	ecf=2.022E-11 # gamma=1.4 for Chandra Cy18

part1=np.zeros(len(centers1))

#tots,bkgs,expos=4,0.057433,2e4
tots,bkgs,expos=9,0.41,1e5

T=bkgs+(10**centers1/ecf)*expos*0.7
prob0=scipy.stats.distributions.poisson.pmf(tots,T)
prob1=scipy.stats.distributions.poisson.pmf(tots,T)*(dndlogs(centers1)*10**ds1)

prob0=np.array(prob0)
prob1=np.array(prob1)

prob0=prob0/np.sum(prob0)
prob1=prob1/np.sum(prob1)

f,ratio=np.genfromtxt(wd+'ratio_dnds_georg.dat',unpack=True)
corr_prob=prob0*ratio

f,p=np.genfromtxt(wd+'georgakakis_pdf.dat',unpack=True)
f=np.log10(f)

plt.figure()
plt.plot(centers1,prob0,'k-',label=r'No $f_x^\beta$')
plt.plot(centers1,corr_prob,'k--',linewidth=2,label=r'What it should be')
plt.plot(centers1,prob1,'r-',label=r'Yes $f_x^\beta$')
plt.plot(f,p,'b*',linewidth=2,label='Georgakakis+08')
#plt.xscale('log')
plt.yscale('log')
#plt.axis([3e-17,1e-12,1e-7,1])
plt.axis([-16.5,-12,1e-7,1])
plt.legend()
plt.tight_layout()
plt.show()
sys.exit()


bins00=np.logspace(np.log10(5e-17),np.log10(1e-12),51)
centers00=list((bins00[i+1]+bins00[i])/2. for i in range(0,len(bins00)-1))
centers00=np.array(centers00) # fluxes array
band='broad'
b1=-1.58
b2=-2.59
fb=2e-14

if band=='broad':
	ecf=1.411E-11 # gamma=1.4 for Chandra Cy18
elif band=='soft':
	ecf=8.707E-12 # gamma=1.4 for Chandra Cy18
elif band=='hard':
	ecf=2.022E-11 # gamma=1.4 for Chandra Cy18

part1=np.zeros(len(centers00))

#tots,bkgs,expos=4,0.057433,2e4
tots,bkgs,expos=9,0.41,1e5

T=bkgs+(centers00/ecf)*expos*0.7
prob0=scipy.stats.distributions.poisson.pmf(tots,T)
prob1=scipy.stats.distributions.poisson.pmf(tots,T)*dnds(centers00)

#prob0,prob1=[],[]
#for j in range(len(centers00)):
#	T=bkgs+(centers00[j]/ecf)*expos*0.7
#	prob0.append(scipy.stats.distributions.poisson.pmf(tots,T))
#	prob1=scipy.stats.distributions.poisson.pmf(tots,T)*dnds(centers00)
#	if centers00[j] >= fb:
#		prob1.append(scipy.stats.distributions.poisson.pmf(tots,T)*(fb/1e-14)**(b1-b2)*(centers00[j]/1e-14)**b2)
#	else:
#		prob1.append(scipy.stats.distributions.poisson.pmf(tots,T)*(centers00[j]/1e-14)**b1)	

prob0=np.array(prob0)
prob1=np.array(prob1)

prob0=prob0/np.sum(prob0)
prob1=prob1/np.sum(prob1)

f,ratio=np.genfromtxt(wd+'ratio_dnds_georg.dat',unpack=True)
corr_prob=prob0*ratio

f,p=np.genfromtxt(wd+'georgakakis_pdf.dat',unpack=True)
p=10.**p

plt.figure()
plt.plot(centers00,prob0,'k-',label=r'No $f_x^\beta$')
plt.plot(centers00,corr_prob,'k--',linewidth=2,label=r'What it should be')
plt.plot(centers00,prob1,'r--',label=r'Yes $f_x^\beta$')
plt.plot(f,p,'b*',linewidth=2,label='Georgakakis+08')
plt.xscale('log')
plt.yscale('log')
plt.axis([3e-17,1e-12,1e-7,1])
plt.legend()
plt.tight_layout()
plt.show()
sys.exit()

# Normalize them (this step should be correct)
prob0=prob0/np.sum(prob0)
prob1=prob1/np.sum(prob1)

# Store in the sum
part0=part0+prob0
part1=part1+prob1

# Part 0/1 contains the sum of the PDFs
part0b=part0/(sens*nsim)
part1b=part1/(sens*nsim)

part0c=list(part0b[i]/(bins00[i+1]-bins00[i]) for i in range(len(bins00)-1))
part1c=list(part1b[i]/(bins00[i+1]-bins00[i]) for i in range(len(bins00)-1))

plt.figure()
plt.plot(centers00,logp1,'r-',label='Mine')
plt.plot(centers00,logp2,'r--',label='Mine, with fx')
plt.plot(f,p,'k-',linewidth=2,label='Georgakakis+08')
plt.axvline(x=1.35e-15)
plt.xscale('log')
#plt.axis([3e-17,2e-12,-5.5,0])
plt.legend()
plt.show()

sys.exit()
'''

#### OLD SENSITIVITY PLOTS, XBOOTES VS CDWFS (?)
'''
(f0,a0)=np.genfromtxt(wd+'cdwfs_broad_sens.dat',unpack=True)
(f1,a1)=np.genfromtxt(wd+'cdwfs_broad_sens_georgakakis.dat',unpack=True)
(f1b,a1b)=np.genfromtxt(wd+'xbootes_broad_sens_georgakakis.dat',unpack=True)
(f1c,a1c)=np.genfromtxt(wd+'xbootes_2-10keV_sens_geo.dat',unpack=True)
f1c=f1c*1.024 #convert from 2-10 to 0.5-7 keV, Gamma 1.4
(f2,a2)=np.genfromtxt(wd+'skycov.dat',unpack=True)
plt.figure()
#plt.plot(f0,a0,'k-',label='old')
#plt.plot(f1,a1,'b-',label='cdwfs_new')
plt.plot(f1b,a1b,'b--',label='xbootes_new')
plt.plot(f1c,a1c,'r--',label='xbootes_geo')
#plt.plot(f2,a2,'g-',label='sim')
plt.xscale('log')
plt.yscale('log')
plt.axis([3.e-17,3.e-13,1.e-3,12])
plt.grid()
plt.legend()
plt.show()
sys.exit()
'''

#### DIVIDE EACH 4X4 PIXEL OF EXPOMAP BY 16 AND DIVIDE PSFMAP BY EXPOMAP TO HAVE WEIGHTED AVERAGE
'''
field='xbootes'
band='soft'

with fits.open(wd+'new_mosaics_detection/'+field+'_'+band+'_expomap_4reb.fits', mode='update') as hdul:
	# Change something in hdul.
	hdul[0].data=hdul[0].data/16.0
	hdul.flush()  # changes are written back

expmap=fits.open(wd+'new_mosaics_detection/'+field+'_'+band+'_expomap_4reb.fits')
exp=expmap[0].data

s.call('cp '+wd+'psfmaps/'+field+'_'+band+'_r90sq-x-exp_4reb.fits '+wd+'psfmaps/'+field+'_'+band+'_r90sq-x-exp_4reb-cp.fits',shell=True)
with fits.open(wd+'psfmaps/'+field+'_'+band+'_r90sq-x-exp_4reb-cp.fits', mode='update') as hdul:
	# Change something in hdul.
	hdul[0].data=hdul[0].data/(16.0*exp)
	hdul[0].data[np.isnan(hdul[0].data)]=0.0
	hdul.flush()  # changes are written back
sys.exit()
'''

#### PLOTS WITH THE MERGED X-RAY CATALOG
'''
cat=fits.open(wd+'new_mosaics_detection/cdwfs_merged_cat1_exp-psf.fits')
ra=cat[1].data['RA']
dec=cat[1].data['DEC']
r90f=cat[1].data['R90_F']
r90s=cat[1].data['R90_S']
r90h=cat[1].data['R90_H']
totf=cat[1].data['TOT_F']
tots=cat[1].data['TOT_S']
toth=cat[1].data['TOT_H']
bkgf=cat[1].data['BKG_F']
bkgs=cat[1].data['BKG_S']
bkgh=cat[1].data['BKG_H']
ctsf=cat[1].data['NET_F']
ectsf_u=cat[1].data['E_NET_F_+']
ctss=cat[1].data['NET_S']
ectss_u=cat[1].data['E_NET_S_+']
ctsh=cat[1].data['NET_H']
ectsh_u=cat[1].data['E_NET_H_+']
poserr=cat[1].data['POS_ERR']
hr=cat[1].data['HR']
ehr_n=cat[1].data['E_HR_-']
ehr_p=cat[1].data['E_HR_+']
fluxf0=cat[1].data['FLUX_F']
efluxfp0=cat[1].data['E_FLUX_F_+']
fluxs0=cat[1].data['FLUX_S']
efluxsp0=cat[1].data['E_FLUX_S_+']
fluxh0=cat[1].data['FLUX_H']
efluxhp0=cat[1].data['E_FLUX_H_+']
probf0=cat[1].data['PROB_F']
probs0=cat[1].data['PROB_S']
probh0=cat[1].data['PROB_H']

print('The CDWFS catalog contains',len(ra),'X-ray point sources.')
#poserr=poserr[fluxf>0]
#ectsf_u=ectsf_u[ctsf>0]
#ctsf=ctsf[ctsf>0]
#ectss_u=ectss_u[ctss>0]
#ctss=ctss[ctss>0]
#ectsh_u=ectsh_u[ctsh>0]
#ctsh=ctsh[ctsh>0]
#probf=probf[fluxf>0]

#hrb=hr[fluxf>0]
#fluxf=fluxf[fluxf>0]
#fluxs=fluxs[fluxs>0]
#fluxh=fluxh[fluxh>0]

#fluxf2=fluxf[hrb!=-99]
#hr2=hrb[hrb!=-99]


# X-ray Flux distributions
fluxf=fluxf0[efluxfp0!=0]
fluxs=fluxs0[efluxsp0!=0]
fluxh=fluxh0[efluxhp0!=0]
probf = probf0[efluxfp0!=0]
probs = probs0[efluxsp0!=0]
probh = probh0[efluxhp0!=0]

fluxf=fluxf[probf < 1e-4]
fluxs=fluxs[probs < 1e-4]
fluxh=fluxh[probh < 1e-4]

print('F flux range:',np.min(fluxf),np.max(fluxf))
print('S flux range:',np.min(fluxs),np.max(fluxs))
print('H flux range:',np.min(fluxh),np.max(fluxh))

print(len(fluxf[fluxf > 5e-13]),'sources in the F band with F>5e-13')
print(len(fluxs[fluxs > 5e-13]),'sources in the S band with F>5e-13')
print(len(fluxh[fluxh > 5e-13]),'sources in the H band with F>5e-13')

bins=np.logspace(np.log10(1e-16),np.log10(5e-13),40)

f,ax=plt.subplots(3,sharex=True,sharey=True)
ax[0].hist(fluxf,bins=bins,histtype='step',linewidth=3,color='k')
ax[0].axvline(x = np.min(fluxf),linewidth = 3, linestyle ='dashed', color='k')
#ax[0].axvline(x = np.max(fluxf),linewidth = 3, linestyle ='dashed', color='k')
ax[0].tick_params(which='major',direction='inout',length=7,labelsize=13)
ax[0].tick_params(which='minor',direction='inout',length=4,labelsize=13)
ax[0].annotate('0.5-7 keV',xy=(7e-14,500),fontsize=13)

ax[1].hist(fluxs,bins=bins,histtype='step',linewidth=3,color='r')
ax[1].axvline(x = np.min(fluxs),linewidth = 3, linestyle ='dashed', color='r')
#ax[1].axvline(x = np.max(fluxs),linewidth = 3, linestyle ='dashed', color='r')
ax[1].tick_params(which='major',direction='inout',top=True,length=7,labelsize=13)
ax[1].tick_params(which='minor',direction='inout',top=True,length=4,labelsize=13)
ax[1].annotate('0.5-2 keV',xy=(7e-14,400),fontsize=13)

ax[2].hist(fluxh,bins=bins,histtype='step',linewidth=3,color='b')
ax[2].axvline(x = np.min(fluxh),linewidth = 3, linestyle ='dashed', color='b')
#ax[2].axvline(x = np.max(fluxh),linewidth = 3, linestyle ='dashed', color='b')
ax[2].set_xlabel(r'Flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=12)
ax[2].set_xscale('log')
ax[2].tick_params(which='major',direction='inout',top=True,length=8,labelsize=13)
ax[2].tick_params(which='minor',direction='inout',top=True,length=4,labelsize=13)
ax[2].annotate('2-7 keV',xy=(7e-14,400),fontsize=13)
plt.tight_layout()
plt.subplots_adjust(hspace=0)
plt.show()

print('min(pos_err):',min(poserr),'max(pos_err):',max(poserr))
print('median(pos_err):',np.median(poserr))

bins=np.logspace(np.log10(min(poserr)),np.log10(max(poserr)),50)
xlabels=[0.1,1]
ylabels=[1,10,100]

a,b=np.histogram(poserr,bins=bins)
bc=list((b[i+1]+b[i])/2. for i in range(len(b)-1))
bc=np.array(bc)
cum=np.cumsum(a)/np.sum(a)

#print(np.min(abs(cum-0.9)))
print(cum[abs(bc-1)==np.min(abs(bc-1))],'% of sources have pos err less than 1 arcsec')

print('68% of sources have pos err less than',bc[abs(cum-0.68)==np.min(abs(cum-0.68))])
print('90% of sources have pos err less than',bc[abs(cum-0.9)==np.min(abs(cum-0.9))])
print('99% of sources have pos err less than',bc[abs(cum-0.99)==np.min(abs(cum-0.99))])

f,ax=plt.subplots(figsize=[6,3])
ax.hist(poserr,bins=bins,histtype='step',hatch='/',linewidth=3, color='k')
ax.set_xscale('log')
ax.set_xlabel(r'Positional error ["]',fontsize=12)
ax.set_ylabel(r'Number of sources',fontsize=12)
ax.tick_params(axis ='both', which='major', right=True, top=True, direction ='in', length = 6, labelsize = 12)
ax.tick_params(axis ='both', which='minor', right=True, top=True, direction ='in', length = 4, labelsize = 12)
ax.set_xticks(xlabels)
ax.set_xticklabels(xlabels)

ax2=ax.twinx()
ax2.plot(bc,cum,'b-')
ax2.tick_params(axis ='both', which='major', right=True, top=True, direction ='in', length = 6, labelsize = 12)
ax2.tick_params(axis ='both', which='minor', right=True, top=True, direction ='in', length = 4, labelsize = 12)
ax2.set_ylabel(r'Cumulative fraction',fontsize=12)

plt.tight_layout()
#plt.show()
plt.savefig(wd+'poserr.pdf',format='pdf')

sys.exit()

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
plt.axvline(x=1e-4,linestyle='dashed')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Prob of being spurious',fontsize=13)
plt.ylabel(r'0.5-7 keV flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=13)
plt.axis([1e-30,1,5e-17,5e-13])
plt.tick_params(which='major',labelsize=13)
plt.tight_layout()
plt.show()
#plt.savefig(wd+'cdwfs_full_flux.pdf',format='pdf')

#this cuts out the UpperLims
ctsf=ctsf[ectsf_u!=0.0]
ctss=ctss[ectss_u!=0.0]
ctsh=ctsh[ectsh_u!=0.0]

bins=np.logspace(np.log10(min(ctsf)),np.log10(max(ctsf)),30)
plt.figure()
plt.hist(ctsf,bins=bins,label='F',histtype='step')
plt.hist(ctss,bins=bins,label='S',histtype='step')
plt.hist(ctsh,bins=bins,label='H',histtype='step')
plt.xlabel('Net cts',fontsize=13)
plt.tick_params(which='major',labelsize=13)
plt.xscale('log')
plt.legend()
plt.tight_layout()
plt.show()
#plt.savefig(wd+'cdwfs_net_distr.pdf',format='pdf')
sys.exit()
'''

#### SECTION ON THE APEC AND DIFFUSE SOFT BACKGROUND
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
'''

#### ANOTHER SPURIOUS FRACTION PLOT (OLD)
'''
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

#### MAKE THE CHUNG+14 REGION FILE
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
'''

#### BREAKDOWN OF THE COMPOSITION OF THE X-RAY CATALOG IN F,S,H COMBINATIONS
'''
cat=fits.open(wd+'new_mosaics_detection/cdwfs_merged_cat1_exp-psf.fits')
ra=cat[1].data['RA']
dec=cat[1].data['DEC']
r90=cat[1].data['R90_S']
probf=cat[1].data['PROB_F']
probs=cat[1].data['PROB_S']
probh=cat[1].data['PROB_H']
cutf,cuts,cuth=10**(-4.5),10**(-4.3),10**(-4.5)
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
				print('something wrong')
print('fsh:',fsh)
print('fs:',fs)
print('fh:',fh)
print('sh:',sh)
print('f:',f)
print('s:',s)
print('h:',h)
print('Total of combinations:',fsh+fs+fh+sh+f+s+h)
#print(len(totf),len(totf[totf>=30.]))
sys.exit()
'''

#### TEST THE INTERPOLATION OF THE GAMMA - CONVERSION FACTOR AND PLOT THE GAUSSIAN DISTRIB
'''
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

#### OTHER STUFF ON THE APEC AND SOFT DIFFUSE BACKGROUND (OLD)
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

#### CHECK WITH DS9 THAT EXTRACTION REGION FILE IS CENTERED IN ACIS FOV - SOFT BKG RELATED
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

#### OLD STUFF ON RELIABILITY AND SPURIOUS FRACTIONS
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

#### EXTREMLY OLD STUFF ON COMPLETENESS AND RELIABILITY AT 97% AND 99% RELIABILITY (!)
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
plt.hist(d,bins=50,alpha=0.7,density=1,label='Full')
plt.hist(d2,bins=50,alpha=0.7,density=1,label='Soft')
plt.hist(d3,bins=50,alpha=0.7,density=1,label='Hard')
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
'''

#### DANGEROUS - REMOVE STUFF FROM FOLDERS (COMMENTED THE COMMAND JUST TO BE SURE)
'''
obsid=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=1,dtype='str')
for i in range(len(obsid)):
	if len(obsid[i]) == 4:
		stem='0'+obsid[i]
	elif len(obsid[i]) == 3:
		stem='00'+obsid[i]
	elif len(obsid[i]) == 5:
		stem=obsid[i]
	#s.call('rm -f '+wd+'data/'+obsid[i]+'/repro_new_asol/out/*soft*',shell=True)
	
print('Done')
sys.exit()
'''

#### OLD STUFF ON CHECKING THE DOUBLE MATCHES WITH XBOOTES (?)
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

#### DISTRIBUTION OF XBOOTES DETECTIONS MADE BY ME 
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

#### COMPARISON OLD SIMULATIONS AND DATA COUNTS AND PROBABILITIES
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
'''

#### ANALYSIS OF OLD CHUINKS OF MOSAICS IN FULL RESOLUTION, AND MATCH WITH KENTER CATALOG
#### HERE I ALSO SAVED THE FITS FILE VERSION OF THE KENTER CATALOG
'''
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

bins=np.logspace(np.log10(1e-16),np.log10(1e-12),20)

plt.figure()
plt.hist(flux,bins=bins,histtype='step',linewidth=3,color='red')
plt.xscale('log')
plt.xlabel(r'0.5-7 keV flux (erg cm$^{-2}$ s$^{-1}$)',fontsize=13)
plt.ylabel('N',fontsize=13)
plt.tight_layout()
plt.savefig(wd+'cdwfs_broad_flux_distr.pdf',format='pdf',dpi=1000)


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
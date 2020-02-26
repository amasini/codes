import numpy as np
import matplotlib.pyplot as plt
import sys
import subprocess as s
import os
from astropy.io import fits
import scipy.stats.distributions
import time
from scipy.optimize import minimize
from sklearn.utils import resample
from scipy.special import gammainc

# Function to compute distances in arcsec
def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*xx*3600)**2 +((pointa[1]-pointb[1])*3600)**2)
    
# Function to compute dN/dS given parameters
def dnds(fx,params):
	
	b1,b2,fb = params[0],params[1],params[2]
	k,fref = 1.6956e16,1e-14
		
	k1 = k*(fb/fref)**(b1-b2)
	
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
		return aux

# Function to compute the sum of log Likelihood for dN/dS fitting - normalization doesn't 
# matter
def func(params):
	fx = centers00
	k,b1,b2,fb=5e16,params[0],params[1],params[2]
	fref=1e-14
	k1 = k*(fb/fref)**(b1-b2)
	
	if type(fx) == list:
		fx = np.array(fx)
	aux = k*(fx/fref)**b1
	aux[fx > fb] = k1*(fx[fx > fb]/fref)**b2
	
	lnpconi=[]
	for i in range(len(exp)):
		ecf=flux0[i]/cr[i]
		T=bkg[i]+(fx/ecf)*exp[i]*0.9
		
		pb=scipy.stats.distributions.poisson.pmf(tot[i],T)*aux*ds00
		pconi = np.sum(pb)/np.sum(aux*sens*ds00)
		
		lnpconi.append(-np.log(pconi))
	
	lnpconi=np.array(lnpconi)
	return np.sum(lnpconi)

# Build a structure
def build_struct(params):
	s=[]
	for j in range(len(params[0])):
		aux=[]
		for i in range(len(params)):
			aux.append(params[i][j])
		s.append(aux)
	return s

# Function to compute the sum of log Likelihood for dN/dS fitting - normalization doesn't 
# matter - FOR BOOTSTRAPPING
def boot_func(params):
	fx = centers00
	k,b1,b2,fb=5e16,params[0],params[1],params[2]
	fref=1e-14
	k1 = k*(fb/fref)**(b1-b2)
	
	if type(fx) == list:
		fx = np.array(fx)
	aux = k*(fx/fref)**b1
	aux[fx > fb] = k1*(fx[fx > fb]/fref)**b2
	
	lnpconi=[]
	for i in range(len(exp)):
	
		T=boot[i][0]+(centers00/boot[i][1])*boot[i][2]*0.9
		
		pb=scipy.stats.distributions.poisson.pmf(boot[i][3],T)*aux*ds00
		pconi = np.sum(pb)/np.sum(aux*sens*ds00)
		
		lnpconi.append(-np.log(pconi))
	
	lnpconi=np.array(lnpconi)
	return np.sum(lnpconi)

wd='/Users/alberto/Desktop/XBOOTES/'

band='hard'
bootstrap = True
nboot = 100
write_output = True

if band == 'broad':
	cut = 10**(-4.63)
	bb = 'F'
elif band == 'soft':
	cut = 10**(-4.57)
	bb = 'S'
else:
	cut = 10**(-4.40)
	bb = 'H'
	
logcut = np.log10(cut)

# Real data
cat=fits.open(wd+'CDWFS_I-Ks-3.6_v200113.fits')
data=cat[1].data
ra =  data['CDWFS_RA']
dec = data['CDWFS_DEC']
r90 = data['CDWFS_R90_'+bb]
exp = data['CDWFS_EXP_'+bb]
tot = data['CDWFS_TOT_'+bb]
bkg = data['CDWFS_BKG_'+bb]
cr = data['CDWFS_CR_'+bb]
flux = data['CDWFS_FLUX_'+bb]
eflux = data['CDWFS_E_FLUX_'+bb+'_+']
prob = data['CDWFS_PROB_'+bb]
cat.close()

# Apply probability cut (automatically excludes upperlimits)
ra = ra[prob <= cut]
dec = dec[prob <= cut]
r90 = r90[prob <= cut]
exp = exp[prob <= cut]
tot = tot[prob <= cut]
bkg = bkg[prob <= cut]
cr = cr[prob <= cut]
eflux1 = eflux[prob <= cut]
flux0 = flux[prob <= cut]
prob1 = prob[prob <= cut]

what = 'shal'

if what == 'shal':
	center = [217.86,34.42]
	#center = [218.28,33.03]
	#center = [217.6,33.41]
	#center = [217.83,33.7]
elif what == 'deep':
	center = [216.41,35.61]
# Select sources within radius from center
#center = [216.41,35.61] # ~LaLa survey [216.41,35.61]; other center could be [217.86,34.42]
radius = 670 # arcsec 

# Filter sources computing distance to center of interesting region
d=[]
for k in range(len(ra)):
	d.append(distance(center,[ra[k],dec[k]]))

d = np.array(d)

mask = d <= radius

ra = ra[mask]
dec = dec[mask]
exp = exp[mask]
tot = tot[mask]
bkg = bkg[mask]
cr = cr[mask]
eflux1 = eflux1[mask]
flux0 = flux0[mask]
prob1 = prob1[mask]
		
# Filter maps based on selected regions and make sensitivities
'''
map = ['expomap','bkgmap','r90sq','ecfmap']
# Cut Bkgmap,Expomap, file with sources to consider - DEEP field
for i in range(len(map)):
	s.call('dmcopy \''+wd+'new_mosaics_detection/cdwfs_'+band+'_'+map[i]+'_4reb.fits[sky=circle('+str(center[0])+'d,'+str(center[1])+'d,'+str(radius)+'\")]\' '+wd+'test_hard_lognlogs/cdwfs_'+band+'_'+map[i]+'_'+what+'.fits clobber=yes',shell=True)

# Compute sensitivity for these regions following Georgakakis

########### PARAMETERS ############
rpsf='r90'
logcut=-4.40
pthresh=10**(logcut)

rebin_factor=4.
scale=(0.492/3600.)*rebin_factor #pixel size in deg
arcsec2pix=scale*3600.
###################################

print('Doing',what,'in the',band,'band, using '+rpsf+' and '+str(round(logcut,2)))


### Take expo map (exposure has average exposure in 4x4 pixels in s)
expmap=fits.open(wd+'test_hard_lognlogs/cdwfs_'+band+'_expomap_'+what+'.fits')
exp=expmap[0].data
exp[np.isnan(exp)]=0.0 #put nans to 0
expmap.close()

### Take average psfmap squared (in arcsec)
psfmap=fits.open(wd+'test_hard_lognlogs/cdwfs_'+band+'_r90sq_'+what+'.fits')
psf=psfmap[0].data
psfmap.close()

# No need to convert to arcsec if psfmap in arsec is used
sourcearea=np.pi*psf # Use r90
fpsf=0.9

### Take bkg map (bkgmap has SUMMED bkg in 4x4 pixels, like data)
bkgmap=fits.open(wd+'test_hard_lognlogs/cdwfs_'+band+'_bkgmap_'+what+'.fits')
bkg=bkgmap[0].data
backheader=bkgmap[0].header
bkg[np.isnan(bkg)]=0.0 #put nans to 0
bkgmap.close()

### Take energy conversion factors map (weighted average of countrate to flux factor)
ecfmap=fits.open(wd+'test_hard_lognlogs/cdwfs_'+band+'_ecfmap_'+what+'.fits')
ecf=ecfmap[0].data
ecfmap.close()

### Bkg is in cts/pixel; bkg2 is in cts [(cts/arcsec2)*arcsec2]
pixarea=arcsec2pix*arcsec2pix
bkg2=bkg*(sourcearea/pixarea)

pcts=np.zeros_like(bkg2,dtype=float)


# GAMMA INCOMPLETE METHOD
# Given background and threshold probability, find how many counts are needed to have P < threhsold
tin=time.time()
tot0=np.arange(1,101)
f=np.logspace(np.log10(1e-17),np.log10(2e-10),101)
totprob=np.zeros_like(f)
for i in range(bkg2.shape[0]):
	for j in range(bkg2.shape[1]):
		if exp[i][j] != 0.0:
			trial=gammainc(tot0,bkg2[i][j])

			pcts[i][j]=np.min(tot0[trial < pthresh])

			T=bkg2[i][j]+f*exp[i][j]/ecf[i][j]*fpsf

			prob=gammainc(pcts[i][j],T)
			totprob=totprob+prob
		else:
			pcts[i][j]=0
				
print((time.time()-tin),'seconds for the georkakais map.')

totprob=totprob*2.988e-7 # convert to area (1 4x4 pix is 2.988e-7 deg2)

### Write out result
w=open(wd+'test_hard_lognlogs/cdwfs_'+band+'_sens_geo_'+what+'.dat','w')
for j in range(len(f)):
	w.write(str(f[j])+' \t '+str(totprob[j])+'\n')
w.close()


# FOLLOWING CIVANO, EASIER METHOD

# Given background and threshold probability, find how many counts are needed to have P < threhsold
tin=time.time()
tot0=np.arange(1,101)
f=np.logspace(np.log10(1e-17),np.log10(2e-10),101)
fcen = list((f[i+1]+f[i])/2. for i in range(len(f)-1))
fcen = np.array(fcen)
#totprob=np.zeros_like(f)
for i in range(bkg2.shape[0]):
	for j in range(bkg2.shape[1]):
		if exp[i][j] != 0.0:
			trial=gammainc(tot0,bkg2[i][j])
			
			pct = np.min(tot0[trial < pthresh])
			
			pcts[i][j] = ((pct - bkg2[i][j])*ecf[i][j])/(exp[i][j]*fpsf)
			
		else:
			pcts[i][j]=0
			
print((time.time()-tin)/60.,'minutes for the normal map.')

a,b = np.histogram(pcts, bins = f)
totprob = np.cumsum(a)

totprob=totprob*2.988e-7 # convert to area (1 4x4 pix is 2.988e-7 deg2)

### Write out result
w=open(wd+'test_hard_lognlogs/cdwfs_'+band+'_sens_civ_'+what+'.dat','w')
for j in range(len(fcen)):
	w.write(str(fcen[j])+' \t '+str(totprob[j])+'\n')
w.close()
'''

# Plot the sensitivities
'''
f0,a0 = np.genfromtxt(wd+'test_hard_lognlogs/cdwfs_hard_sens_geo_deep.dat',unpack=True)
f1,a1 = np.genfromtxt(wd+'test_hard_lognlogs/cdwfs_hard_sens_geo_shal.dat',unpack=True)

plt.figure()
plt.plot(f0,a0,'k-',label='Deep')
plt.plot(f1,a1,'b-',label='Shal')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()
'''

# Apply a flux cut - if I cut in flux like this, I lose information on the faint end of the lognlogs and on the flux break.
flux_cut = 5e-13
mask = flux0 < flux_cut

ra = ra[mask]
dec = dec[mask]
exp = exp[mask]
tot = tot[mask]
bkg = bkg[mask]
cr = cr[mask]
eflux1 = eflux1[mask]
flux0 = flux0[mask]
prob1 = prob1[mask]

print(np.min(flux0),np.max(flux0))
print('Using',len(exp),'sources for the following computation.')

bins00=np.logspace(np.log10(5e-17),np.log10(5e-13),51) # Check the limits here based on detected sources?
centers00=list((bins00[i+1]+bins00[i])/2. for i in range(0,len(bins00)-1))
centers00=np.array(centers00)
ds00 = list((bins00[i+1]-bins00[i]) for i in range(0,len(bins00)-1))
ds00 = np.array(ds00)

# EASY WAY

# Take sensitivity curve made with Civano method
(rawf,rawa)=np.genfromtxt(wd+'test_hard_lognlogs/cdwfs_'+band+'_sens_civ_'+what+'.dat',unpack=True)
sens=np.interp(centers00,rawf,rawa)

part1,bc = np.histogram(flux0, bins = bins00)
part1b=part1/sens

cumpart0=list(reversed(np.cumsum(list(reversed(part1b)))))

if write_output == True:
		w=open(wd+'test_hard_lognlogs/cdwfs_lognlogs_'+band+'_'+what+'_easy.dat','w')
		for jj in range(len(centers00)):
			w.write(str(centers00[jj])+' \t '+str(cumpart0[jj])+' \n')
		w.close()

geos,geon = np.genfromtxt(wd+'geo_lognlogs_'+band+'.txt',unpack=True)
#geos= 6.887E-01*geos # convert from 2-10 to 2-7, Gamma = 1.4
if band == 'hard':
	geos= 0.75*geos # convert from 2-10 to 2-7, Gamma = 1.8
elif band == 'broad':
	geos = 8.455E-01*geos # convert from 0.5-10 to 0.5-7, Gamma = 1.8

ff,yy = np.genfromtxt(wd+'test_hard_lognlogs/cdwfs_lognlogs_'+band+'_deep_easy.dat',unpack=True)

plt.figure()
plt.plot(geos,geon,'k-',label='Georgakakis')
plt.plot(centers00,cumpart0,'ro',ms=5,label='CDWFS-Shal')
plt.plot(centers00,yy,'bo',ms=5,label='CDWFS-Deep')
plt.xlabel(r''+band+' band flux (erg cm$^{-2}$ s$^{-1}$)')
plt.ylabel(r'$N(>S)$ (deg$^{-2}$)')
plt.xscale('log')
plt.yscale('log')
#plt.axis([5e-18,1e-12,0.1,200])
plt.legend()
plt.show()


# Check the interpolation
#plt.figure()
#plt.plot(rawf,rawa,'ko')
#plt.plot(centers00,sens,'r+')
#plt.xscale('log')
#plt.show()
#sys.exit()

######################################
# RECOVER THE dN/dS WITH THE PARAMETERS
'''
print('Now fitting the dN/dS...')

# Take sensitivity curve made with Georgakakis method
(rawf,rawa)=np.genfromtxt(wd+'test_hard_lognlogs/cdwfs_'+band+'_sens_geo_'+what+'.dat',unpack=True)
sens=np.interp(centers00,rawf,rawa)

if bootstrap == True:
	# Compute dN/dS following Georgakakis+08 and use nboot bootstrap for uncertainties
	tin=time.time()
	ecf = flux0/np.array(cr)

	pars=[bkg,ecf,exp,tot]
	
	data = build_struct(pars)
	bootstrap_dnds,p0,p1,p2,bootstrap_lognlogs=[],[],[],[],[]
	print('Starting bootstrap...')
	for bb in range(nboot):
		print(bb+1,'/',nboot)

		boot = resample(data, replace=True, n_samples=len(data), random_state=None)
		
		# Minimize the -Likelihood function starting from a guess of the parameters
		guess=[-1,-2,8e-15]
		res=minimize(boot_func, guess, method='nelder-mead')
		#print(res)
		
		parameters=res.x
		
		p0.append(parameters[0])
		p1.append(parameters[1])
		p2.append(parameters[2])

		part1=np.zeros_like(centers00)
		for i in range(len(exp)):	

			# Compute the PDFs
			T=boot[i][0]+(centers00/boot[i][1])*boot[i][2]*0.9
			prob1=scipy.stats.distributions.poisson.pmf(boot[i][3],T)*(dnds(centers00,parameters)*ds00)

			# Normalize them
			prob1=prob1/np.sum(prob1)

			# Store in the sum
			part1=part1+prob1

		# Part1 contains the sum of the PDFs
		part1b=part1/sens

		cumpart1b=list(reversed(np.cumsum(list(reversed(part1b)))))
		bootstrap_lognlogs.append(cumpart1b)

		# This is the effective dN/dS (sources/deg2/flux_bin)
		part1c=part1b/ds00
		print(part1c)
		if part1c[0] == np.nan:
			sys.exit()
		bootstrap_dnds.append(part1c)

	mu_dnds,sigma_dnds=[],[]
	mu_lognlogs,sigma_lognlogs=[],[]
	for i in range(len(bootstrap_dnds[0])):
		a,b=[],[]
		for j in range(len(bootstrap_dnds)):
			a.append(bootstrap_dnds[j][i])
		mu_dnds.append(np.median(a))
		sigma_dnds.append(np.std(a))

		for k in range(len(bootstrap_lognlogs)):
			b.append(bootstrap_lognlogs[k][i])
		mu_lognlogs.append(np.median(b))
		sigma_lognlogs.append(np.std(b))

	b1=np.median(p0)
	eb1=np.std(p0)
	b2=np.median(p1)
	eb2=np.std(p1)
	fb=np.median(p2)
	efb=np.std(p2)

	print('/'*10)
	print(b1,eb1)
	print(b2,eb2)
	print(fb,efb)
	print('/'*10)

	bfit=dnds(centers00,[b1,b2,fb])
	bfit_lognlogs = list(reversed(np.cumsum(list(reversed(bfit*ds00)))))

	print(round(float(time.time()-tin)/60.,1),' minutes for the bootstrap.')
	
	cumpart0=mu_lognlogs
	e_cumpart0 = sigma_lognlogs
	
	if write_output == True:
		w=open(wd+'test_hard_lognlogs/cdwfs_dnds_'+band+'_'+what+'.dat','w')
		for jj in range(len(mu_dnds)):
			w.write(str(centers00[jj])+' \t '+str(mu_dnds[jj])+' \t '+str(sigma_dnds[jj])+'\n')
		w.close()
	
		w=open(wd+'test_hard_lognlogs/cdwfs_dnds_bfit-pars_'+band+'_'+what+'.dat','w')
		w.write('Beta1 \t eBeta1 \t Beta2 \t eBeta2 \t Fb \t eFb \n')
		w.write(str(b1)+' \t '+str(eb1)+' \t '+str(b2)+' \t '+str(eb2)+' \t '+str(fb)+' \t '+str(efb)+'\n')
		w.close()
	
		w=open(wd+'test_hard_lognlogs/cdwfs_lognlogs_'+band+'_'+what+'.dat','w')
		for jj in range(len(mu_lognlogs)):
			w.write(str(centers00[jj])+' \t '+str(mu_lognlogs[jj])+' \t '+str(sigma_lognlogs[jj])+'\n')
		w.close()

'''

x,y,ey = np.genfromtxt(wd+'test_hard_lognlogs/cdwfs_lognlogs_hard_shal.dat',unpack=True)
x,y1,ey1 = np.genfromtxt(wd+'test_hard_lognlogs/cdwfs_lognlogs_hard_deep.dat',unpack=True)

geos,geon = np.genfromtxt(wd+'geo_lognlogs_'+band+'.txt',unpack=True)
if band == 'hard':
	geos= 0.75*geos # convert from 2-10 to 2-7, Gamma = 1.8
	#geos= 6.887E-01*geos # convert from 2-10 to 2-7, Gamma = 1.4
elif band == 'broad':
	geos = 8.455E-01*geos # convert from 0.5-10 to 0.5-7, Gamma = 1.8

plt.figure()
if bootstrap == True:
	#plt.errorbar(centers00,cumpart0,yerr=e_cumpart0,color='b',marker='.',ms=10,label='CDWFS-Shal')
	plt.errorbar(centers00,y,yerr=ey,color='r',marker='.',ms=10,label='CDWFS-Shal')
	plt.errorbar(centers00,y1,yerr=ey1,color='b',marker='.',ms=10,label='CDWFS-Deep')
else:
	plt.plot(centers00,cumpart0,color='r',marker='.',ms=10,label='CDWFS')
plt.plot(geos,geon,'k-',ms=10,label='Georgakakis')
plt.xlabel(r''+band+' band flux (erg cm$^{-2}$ s$^{-1}$)')

plt.ylabel(r'$N(>S)$ (deg$^{-2}$)')
plt.axis([5e-17,1e-12,0.1,3e4])
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()


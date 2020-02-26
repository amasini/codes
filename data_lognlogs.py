# Logn-logS of the data
import numpy as np
import sys
from astropy.io import fits
import matplotlib.pyplot as plt
import time
import scipy.stats.distributions
from scipy.optimize import minimize
#import numdifftools as nd
from sklearn.utils import resample

# Function of the differential number counts
def dnds(fx,params):
	
	b1,b2,fb = params[0],params[1],params[2]
	k,fref = 1e16,1e-14
		
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

# Function to maximize likelihood L (= minimize -L)
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
	
		T=boot[i][0]+(centers00/boot[i][1])*boot[i][2]*0.9
		
		pb=scipy.stats.distributions.poisson.pmf(boot[i][3],T)*aux*ds00
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

wd = '/Users/alberto/Desktop/XBOOTES/'	
band = 'broad'

#define cuts (may be revised based on spurious_sources.py)
if band=='broad':
	cut=10**(-4.63)
	bb = 'F'
elif band=='soft':
	cut=10**(-4.57)
	bb = 'S'
elif band=='hard':
	cut=10**(-4.40)
	bb = 'H'

logcut = np.log10(cut)

cat1=fits.open(wd+'CDWFS_I-Ks-3.6_v200113.fits')

ra_d=cat1[1].data['CDWFS_RA']
dec_d=cat1[1].data['CDWFS_DEC']
tot=cat1[1].data['CDWFS_TOT_'+bb]
bkg=cat1[1].data['CDWFS_BKG_'+bb]
exp=cat1[1].data['CDWFS_EXP_'+bb]
prob=cat1[1].data['CDWFS_PROB_'+bb]
r90=cat1[1].data['CDWFS_R90_'+bb]
cr=cat1[1].data['CDWFS_CR_'+bb]
flux_d=cat1[1].data['CDWFS_FLUX_'+bb]

cat1.close()

# Cut detected sample at a given probability threshold 
ra_d=ra_d[prob<=cut]
dec_d=dec_d[prob<=cut]
tot=tot[prob<=cut]
bkg=bkg[prob<=cut]
exp=exp[prob<=cut]
cr=cr[prob<=cut]
r90=r90[prob<=cut]
flux_d=flux_d[prob<=cut]
prob=prob[prob<=cut]

ecf=flux_d/cr

# Define the flux range and bins
bins00=np.logspace(np.log10(5e-16),np.log10(1e-12),51)
centers00=list((bins00[i+1]+bins00[i])/2. for i in range(0,len(bins00)-1))
ds00=list((bins00[i+1]-bins00[i]) for i in range(0,len(bins00)-1))
centers00=np.array(centers00)
ds00=np.array(ds00)

# Take the sensitivity made with Georgakakis method (sens.py code)
fl2,ar2=np.genfromtxt(wd+'cdwfs_'+band+'_sens_'+str(round(logcut,1))+'_geo.dat',unpack=True)
sens=np.interp(centers00,fl2,ar2)

#######################################################
# OUTPUT LOGNLOGS

# Minimize the -Likelihood function starting from a guess of the parameters
#guess=[-1,-2,5e-15]
#res=minimize(func, guess, method='nelder-mead')
#print(res)
#parameters=res.x

#H=nd.Hessian(func)(res.x)
#print(H)

# Compute dN/dS following Georgakakis+08 and use 100 bootstrap for uncertainties
pars=[bkg,ecf,exp,tot]
data = build_struct(pars)
bootstrap,p0,p1,p2=[],[],[],[]
nboot=10
for bb in range(nboot):
	
	boot = resample(data, replace=True, n_samples=len(data), random_state=None)
	
	# Minimize the -Likelihood function starting from a guess of the parameters
	guess=[-1,-2,5e-15]
	res=minimize(func, guess, method='nelder-mead')
	print(res)
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

	#print('I cycled on',len(exp),'sources, and the total number of sources in the sum is',round(np.sum(part1),0))

	# Part1 contains the sum of the PDFs
	part1b=part1/sens

	# This is the effective dN/dS (sources/deg2/flux_bin)
	#part1c=list(part1b[i]/(bins00[i+1]-bins00[i]) for i in range(len(bins00)-1))
	part1c=part1b/ds00

	bootstrap.append(part1c)

mu,sigma=[],[]
for i in range(len(bootstrap[0])):
	a=[]
	for j in range(len(bootstrap)):
		a.append(bootstrap[j][i])
	mu.append(np.median(a))
	sigma.append(np.std(a))

mu=np.array(mu)
sigma=np.array(sigma)

b1=np.median(p0)
eb1=np.std(p0)
b2=np.median(p1)
eb2=np.std(p1)
fb=np.median(p2)
efb=np.std(p2)

print('/'*10)
print(b1,eb1/np.sqrt(nboot))
print(b2,eb2/np.sqrt(nboot))
print(fb,efb/np.sqrt(nboot))
print('/'*10)

bfit=dnds(centers00,[b1,b2,fb])
bfit/ds00
#print(bfit)

cumpart1=list(reversed(np.cumsum(list(reversed(part1b)))))
#######################################################

# Plot the results
f,(ax1,ax2)=plt.subplots(2,1,sharex=True,figsize=[7,9])

ax1.plot(centers00,cumpart1,'c-',linewidth=2,label=r'Georgakakis method, with fx$^{\beta}$')
ax1.set_ylabel(r'N(>S) (deg$^{-2}$)',fontsize=13)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.axis([5e-17,2e-12,1,2e4])
ax1.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=12)

ax2.errorbar(centers00,mu,yerr=sigma,fmt='.',color='k',label=r'Georgakakis method, with fx$^{\beta}$')
ax2.plot(centers00,5*bfit,'r--')
ax2.set_ylabel(r'dN/dS (deg$^{-2}$ [erg cm$^{-2}$ s$^{-1}$]$^{1.5}$)',fontsize=13)
ax2.set_xlabel(r'S (erg cm$^{-2}$ s$^{-1}$)',fontsize=13)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.axis([5e-17,2e-12,1e11,1e22])
#ax2.axis([5e-17,2e-12,5e-22,2e-18])
ax2.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=12)

ax1.legend()
ax2.legend()
plt.subplots_adjust(hspace=0)
plt.show()
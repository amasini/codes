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

band='broad'
cut = 5e-5
bootstrap = True
nboot = 20

cat=fits.open(wd+'new_mosaics_detection/cdwfs_merged_cat1_exp-psf.fits')
data=cat[1].data
if band == 'broad':
	exp = data['EXP_F']
	tot = data['TOT_F']
	bkg = data['BKG_F']
	cr = data['CR_F']
	flux0 = data['FLUX_F']
	eflux0 = data['E_FLUX_F_+']
	prob = data['PROB_F']
elif band == 'soft':
	exp = data['EXP_S']
	tot = data['TOT_S']
	bkg = data['BKG_S']
	cr = data['CR_S']
	flux0 = data['FLUX_S']
	eflux0 = data['E_FLUX_S_+']
	prob = data['PROB_S']
else:
	exp = data['EXP_H']
	tot = data['TOT_H']
	bkg = data['BKG_H']
	cr = data['CR_H']
	flux0 = data['FLUX_H']
	eflux0 = data['E_FLUX_H_+']
	prob = data['PROB_H']
cat.close()

exp = exp[prob <= cut]
tot = tot[prob <= cut]
bkg = bkg[prob <= cut]
cr = cr[prob <= cut]
eflux0 = eflux0[prob <= cut]
flux0 = flux0[prob <= cut]
prob = prob[prob <= cut]

exp = exp[eflux0 != 0]
tot = tot[eflux0 != 0]
bkg = bkg[eflux0 != 0]
cr = cr[eflux0 != 0]
flux0 = flux0[eflux0 != 0]
prob = prob[eflux0 != 0]
eflux0 = eflux0[eflux0 != 0]

print('Using',len(exp),'sources for the following computation.')

#ra = ra[flux0 < 1e-13]
#dec = dec[flux0 < 1e-13]
#tot = tot[flux0 < 1e-13]
#exp = exp[flux0 < 1e-13]
#bkg = bkg[flux0 < 1e-13]
#eflux0 = eflux0[flux0 < 1e-13]
#prob = prob[flux0 < 1e-13]
#flux0 = flux0[flux0 < 1e-13]

'''
# Take the CDWFS catalog
cat=fits.open(wd+'CDWFS_I-Ks-3.6_v2.fits')
data=cat[1].data

ra = data['CHA_RA']
dec = data['CHA_DEC']
if band == 'broad':
	exp = data['CHA_EXP_F']
	tot = data['CHA_TOT_F']
	bkg = data['CHA_BKG_F']
	cr = data['CHA_CR_F']
	flux0 = data['CHA_FLUX_F']
	eflux0 = data['CHA_E_FLUX_F_+']
	prob = data['CHA_PROB_F']
elif band == 'soft':
	exp = data['CHA_EXP_S']
	tot = data['CHA_TOT_S']
	bkg = data['CHA_BKG_S']
	cr = data['CHA_CR_S']
	flux0 = data['CHA_FLUX_S']
	eflux0 = data['CHA_E_FLUX_S_+']
	prob = data['CHA_PROB_S']
else:
	exp = data['CHA_EXP_H']
	tot = data['CHA_TOT_H']
	bkg = data['CHA_BKG_H']
	cr = data['CHA_CR_H']
	flux0 = data['CHA_FLUX_H']
	eflux0 = data['CHA_E_FLUX_H_+']
	prob = data['CHA_PROB_H']
cat.close()


ecf = flux0/cr


ra = ra[prob <= cut]
dec = dec[prob <= cut]
tot = tot[prob <= cut]
exp = exp[prob <= cut]
bkg = bkg[prob <= cut]
eflux0 = eflux0[prob <= cut]
flux0 = flux0[prob <= cut]
prob = prob[prob <= cut]

ra = ra[flux0 < 1e-13]
dec = dec[flux0 < 1e-13]
tot = tot[flux0 < 1e-13]
exp = exp[flux0 < 1e-13]
bkg = bkg[flux0 < 1e-13]
eflux0 = eflux0[flux0 < 1e-13]
prob = prob[flux0 < 1e-13]
flux0 = flux0[flux0 < 1e-13]
'''

bins = np.logspace(-16,-11,50)
plt.figure()
plt.hist(flux0, bins = bins)
plt.xscale('log')
plt.show()

#print('The source with a flux above 1E-12 in the soft band is at',ra[flux0 > 1e-12],dec[flux0 > 1e-12])

# Civano+16 figure 14, top left panel
if band =='soft':
	(civf,civn)=np.genfromtxt(wd+'civano_05-2keV_lognlogs.dat',unpack=True)
else:
	(civf,civn)=np.genfromtxt(wd+'civano_2-10keV_lognlogs.dat',unpack=True)
	(civf2,civn2)=np.genfromtxt(wd+'civano_lognlogs_hard.txt',unpack=True)

bins00=np.logspace(np.log10(5e-17),np.log10(5e-12),50) # Check the limits here based on detected sources
#bins00=np.logspace(np.log10(np.min(flux0)),np.log10(np.max(flux0)),50) 
centers00=list((bins00[i+1]+bins00[i])/2. for i in range(0,len(bins00)-1))
centers00=np.array(centers00)
ds00 = list((bins00[i+1]-bins00[i]) for i in range(0,len(bins00)-1))
ds00 = np.array(ds00)

# take sensitivity curve made with Georgakakis method
(rawf,rawa)=np.genfromtxt(wd+'cdwfs_'+band+'_sens_georgakakis_r90_5e-5.dat',unpack=True)
sens=np.interp(centers00,rawf,rawa)

#plt.figure()
#plt.plot(rawf,rawa,'ko')
#plt.plot(centers00,sens,'r+')
#plt.xscale('log')
#plt.show()
#sys.exit()

######################################
# RECOVER THE dN/dS WITH THE PARAMETERS
print('Now fitting the dN/dS...')

if bootstrap == True:
	# Compute dN/dS following Georgakakis+08 and use 100 bootstrap for uncertainties
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
		guess=[-1,-2,5e-15]
		res=minimize(boot_func, guess, method='nelder-mead')

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
	print(b1,eb1/np.sqrt(nboot))
	print(b2,eb2/np.sqrt(nboot))
	print(fb,efb/np.sqrt(nboot))
	print('/'*10)

	bfit=dnds(centers00,[b1,b2,fb])
	bfit_lognlogs = list(reversed(np.cumsum(list(reversed(bfit*ds00)))))

	print(round(float(time.time()-tin)/60.,1),' minutes for the bootstrap.')
	'''
	f,(ax1,ax2)=plt.subplots(2,1,sharex=True,figsize=[7,9])
	
	ax1.plot(centers00,bfit_lognlogs,'r-',linewidth=2, label='Best Fit')
	ax1.errorbar(centers00,mu_lognlogs,yerr=sigma_lognlogs,color='c',marker='.',linewidth=2,label='Recovered')
	
	ax1.set_ylabel(r'N(>S) (deg$^{-2}$)',fontsize=13)
	ax1.set_xscale('log')
	ax1.set_yscale('log')
	ax1.axis([5e-17,2e-12,0.1,2e4])
	ax1.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=12)


	ax2.plot(centers00,bfit*centers00**2.5,'r-',linewidth=2, label='Best Fit')
	ax2.errorbar(centers00,mu_dnds*centers00**2.5,yerr=sigma_dnds*centers00**2.5,color='c',marker='.',linewidth=2,label='Recovered')
	
	ax2.set_ylabel(r'dN/dS (deg$^{-2}$ [erg cm$^{-2}$ s$^{-1}$]$^{1.5}$)',fontsize=13)
	ax2.set_xlabel(r'S (erg cm$^{-2}$ s$^{-1}$)',fontsize=13)
	ax2.set_xscale('log')
	ax2.set_yscale('log')
	#ax2.axis([5e-17,2e-12,1e11,1e21])
	ax2.axis([5e-17,2e-12,1e-21,1e-18])
	ax2.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=12)

	ax1.legend()
	ax2.legend()
	plt.subplots_adjust(hspace=0)
	plt.show()
	'''
	#cumpart0=mu_lognlogs*(centers00/1e-14)**1.5
	#e_cumpart0 = sigma_lognlogs*(centers00/1e-14)**1.5
	
	cumpart0=mu_lognlogs
	e_cumpart0 = sigma_lognlogs
	
	w=open(wd+'cdwfs_lognlogs_'+band+'.dat','w')
	for jj in range(len(cumpart0)):
		w.write(str(centers00[jj])+' \t '+str(cumpart0[jj])+' \t '+str(e_cumpart0[jj])+'\n')
	w.close()
	
	geos,geon = np.genfromtxt(wd+'geo_lognlogs_'+band+'.txt',unpack=True)
	#geos= 6.887E-01*geos # convert from 2-10 to 2-7, Gamma = 1.4
	if band == 'hard':
		geos= 0.75*geos # convert from 2-10 to 2-7, Gamma = 1.8
	elif band == 'broad':
		geos = 8.455E-01*geos # convert from 0.5-10 to 0.5-7, Gamma = 1.8
	
	plt.figure()
	plt.errorbar(centers00,cumpart0,yerr=e_cumpart0,color='r',marker='.',ms=10,label='CDWFS')
	plt.plot(geos,geon,'k-',ms=10,label='Georgakakis')
	#plt.plot(civf,civn,'gs',ms=10,label='Civano')
	if band == 'hard':
		plt.plot(civf2,civn2,'cs',ms=10)
	plt.xlabel(r''+band+' band flux (erg cm$^{-2}$ s$^{-1}$)')
	
	#plt.ylabel(r'$N(>S)\times S^{1.5}$ (deg$^{-2}$)')
	#plt.axis([5e-17,1e-12,0.1,400])
	plt.ylabel(r'$N(>S)$ (deg$^{-2}$)')
	plt.axis([5e-17,1e-12,0.1,3e4])
	
	plt.xscale('log')
	plt.yscale('log')
	plt.legend()
	plt.show()
else:
	
	tin=time.time()
	guess = [-1.3,-2.5,5e-15]
	res = minimize(func, guess, method='nelder-mead', options={'adaptive':True})
	print(round(float(time.time()-tin),1),' seconds for the -likelihood minimization.')
	print('guess:',guess)
	print('output:',res.x)

	pars = res.x
	
	#pars = [-1.3,-2.5,8e-15]
	
	part0=np.zeros_like(centers00)
	part1=np.zeros_like(centers00)

	for i in range(len(exp)):

		observed_flux=flux0[i]

		ecf=observed_flux/cr[i]

		# Compute the PDFs
		T=bkg[i]+(centers00/ecf)*exp[i]*0.9
		prob1=scipy.stats.distributions.poisson.pmf(tot[i],T)*(dnds(centers00,pars)*ds00)

		# Normalize them (this step should be correct)
		prob1=prob1/np.sum(prob1)

		# Store in the sum
		part1=part1+prob1

	# Part 0/1 contains the sum of the PDFs
	part1b=part1/sens
	
	# dN/dS
	part1c=list(part1b[i]/(bins00[i+1]-bins00[i]) for i in range(len(bins00)-1))

	# logN-logS
	#cumpart1=list(reversed(np.cumsum(list(reversed(part1b)))))

	######################################
	
	#flux1 = flux0[prob < 2e-6]
	#part1,bc = np.histogram(flux0, bins = bins00)
	#part1b=part1/sens

	cumpart0=list(reversed(np.cumsum(list(reversed(part1b)))))*(centers00/1e-14)**1.5

	plt.figure()
	plt.plot(centers00,cumpart0,'ko',ms=10,label='CDWFS')
	plt.plot(civf,civn,'gs',ms=10,label='Civano')
	if band == 'hard':
		plt.plot(civf2,civn2,'cs',ms=10)
	plt.xlabel(r''+band+' band flux (erg cm$^{-2}$ s$^{-1}$)')
	plt.ylabel(r'$N(>S)\times S^{1.5}$ (deg$^{-2}$)')
	plt.xscale('log')
	plt.yscale('log')
	#plt.axis([5e-18,1e-12,0.1,200])
	plt.legend()
	plt.show()
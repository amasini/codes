# NEW VERSION, CLEANER AND HOPEFULLY FASTER
import numpy as np
import sys
from astropy.io import fits
import matplotlib.pyplot as plt
import time
from sklearn.utils import resample
from scipy.stats import gaussian_kde
import scipy.stats.distributions
from scipy.optimize import minimize
import matplotlib as mpl
from scipy.stats import kde

'''##################### FUNCTIONS #####################'''
# Function to compute distances in arcsec
def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*xx*3600)**2 +((pointa[1]-pointb[1])*3600)**2)

# Function to choose the best counterpart based on probability (assume I know the position
# of prob in the array)
def choose_best(array):
	if len(array) == 1:
		return array
	else:
		# Assuming prob is element [2]
		probs = list(zip(*array))[2]
		for i in range(len(array)):
			if array[i][2] == np.min(probs):
				return array[i]
					
# Function to compute dN/dS given parameters
def dnds(fx,params):
	
	b1,b2,fb = params[0],params[1],params[2]
	fref = 1e-14
	if band == 'broad':
		k = 5.622e16
	elif band == 'soft':
		k = 1.6956e16
	else:
		k = 5.7313e16
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
	for i in range(len(expos)):
		ecf=flux_out[i]/crs[i]
		T=bkgs[i]+(fx/ecf)*expos[i]*0.9
		
		pb=scipy.stats.distributions.poisson.pmf(tots[i],T)*aux*ds00
		pconi = np.sum(pb)/np.sum(aux*sens*ds00)
		
		lnpconi.append(-np.log(pconi))
	
	lnpconi=np.array(lnpconi)
	return np.sum(lnpconi)

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
	for i in range(len(expos)):
	
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


'''##################### INPUTS #####################'''
# Important parameters
wd = '/Users/alberto/Desktop/XBOOTES/'
band = 'hard'
simfolder = 'sim_indep_22-Nov-19/'
use_outflux = True
fit_dnds = True
bootstrap = True
write_out = False
nboot = 10
nsim = 10

if band == 'broad':
	logcut = -4.63
	inppars=[-1.34,-2.35,8.1e-15] # broad Lehmer's dN/dS params
elif band == 'soft':
	logcut = -4.57
	inppars=[-1.49,-2.48,6.0e-15] # soft
else:
	logcut = -4.40
	inppars=[-1.32,-2.55,6.4e-15] # hard
cut = 10**logcut


'''##################### MAIN #####################'''
print('Starting to match '+str(nsim)+' sims in the '+band+' band...')

# Define the binning for everything 
bins00 = np.logspace(np.log10(5e-17),np.log10(1e-12),51)
centers00 = list((bins00[i+1]+bins00[i])/2. for i in range(0,len(bins00)-1))
ds00 = list((bins00[i+1]-bins00[i]) for i in range(0,len(bins00)-1))
centers00 = np.array(centers00)
ds00 = np.array(ds00)
area = 9.3 

t_in=time.time()

for k in range(nsim):

	input,flux_inp,flux_out=[],[],[]
	
	tots,bkgs,expos,crs=[],[],[],[]
	probabilities = []

	print(k+1, end=" ")
	
	##########################
	##### Start matching #####
	# Take catalog of detected sources
	cat1 = fits.open(wd+simfolder+str(k)+'cdwfs_'+band+'_sim_cat1_exp-psf.fits')
	
	simdata = cat1[1].data
	
	prob = simdata['PROB']
	simdata = np.array(simdata)
	
	# Cut detected sample at a given probability threshold 
	simdata = simdata[prob <= cut]
	print(len(simdata),'above threshold ('+str(round(logcut,2))+')')
	
	# Input sources depend on band and simfolder
	(flux_cdwfs,ra_cdwfs,dec_cdwfs,exp_cdwfs)=np.genfromtxt(wd+'inp_list_sim/poiss_rand_'+band+'_'+simfolder[:-1]+'_filtered_exp.dat',unpack=True,skip_header=1,usecols=[0,1,2,5])	

	ra_cdwfs = ra_cdwfs[exp_cdwfs >= 1800.]
	dec_cdwfs = dec_cdwfs[exp_cdwfs >= 1800.]
	flux_cdwfs = flux_cdwfs[exp_cdwfs >= 1800.]
	
	# This is the array of the input (precise) fluxes
	input.append(flux_cdwfs)
	
	# Sort them to start from the bright ones
	ra_cdwfs=ra_cdwfs[::-1]
	dec_cdwfs=dec_cdwfs[::-1]
	flux_cdwfs=flux_cdwfs[::-1]
	
	unmatched,blendings=0,0

	r90eff,deff=[],[]
	for i in range(len(ra_cdwfs)):
		input_source=[ra_cdwfs[i],dec_cdwfs[i]]
		
		found=0
		counterparts=[]
		count=0
		
		delta = 0.0056 #(0.028 ~100")
		
		filt_sim=simdata[(simdata['RA'] >= ra_cdwfs[i]-delta) & (simdata['RA'] <= ra_cdwfs[i]+delta) & (simdata['DEC'] >= dec_cdwfs[i]-delta) & (simdata['DEC'] <= dec_cdwfs[i]+delta)]
		
		counterparts,distances=[],[]
		for j in range(len(filt_sim['RA'])):
			cdwfs_source=[filt_sim['RA'][j],filt_sim['DEC'][j]]
			
			match_rad=1.1*filt_sim['AV_R90'][j]
			d=distance(input_source,cdwfs_source)
			
			if d <= match_rad: #found a match
				if found==0: #it's the first
					found=1
					
					counterparts.append(filt_sim[j])
					distances.append(d)
					
				else: #it's not the first match	
					count=count+2
					blendings=blendings+1
					
					counterparts.append(filt_sim[j])
					distances.append(d)
		
		if found == 1:
			if len(counterparts) == 1:
				
				flux_inp.append(flux_cdwfs[i])

				deff.append(distances[0])
				
				flux_out.append(counterparts[0][13])
				r90eff.append(counterparts[0][3])
				tots.append(counterparts[0][4])
				bkgs.append(counterparts[0][5])
				expos.append(counterparts[0][9])
				crs.append(counterparts[0][10])
				probabilities.append(counterparts[0][2])
				simdata=np.delete(simdata,np.where(simdata['RA']==counterparts[0][0]),0)
								
			else:
			
				flux_inp.append(flux_cdwfs[i])
				
				counterparts=np.array(counterparts)
				counterparts=choose_best(counterparts)
				
				flux_out.append(counterparts[13])
				r90eff.append(counterparts[3])
				tots.append(counterparts[4])
				bkgs.append(counterparts[5])
				expos.append(counterparts[9])
				crs.append(counterparts[10])
				probabilities.append(counterparts[2])
				simdata=np.delete(simdata,np.where(simdata['RA']==counterparts[0]),0)
				
		else:

			unmatched=unmatched+1
	
	t_out=time.time()
	print(round(float(t_out-t_in),1),' seconds for the match.')
	##### End of matching #####
	###########################

	#flux_inp2 = np.array(flux_inp)
	#flux_out2 = np.array(flux_out)
	#probabilities2 = np.array(probabilities)
	#tots2 = np.array(tots)
	#bkgs2 = np.array(bkgs)
	#snr2 = tots2/bkgs2
	
	# INPUT LOGNLOGS
	quante,bincenters_s = np.histogram(input,bins=bins00)
	quante_perarea = quante/(area)
	quante_perarea_perflux = list(quante_perarea[i]/(bins00[i+1]-bins00[i]) for i in range(len(bins00)-1))
	quante_perarea_perflux = np.array(quante_perarea_perflux)
	ncum_in = list(reversed(np.cumsum(list(reversed(quante_perarea)))))

	# SENSITIVITY CURVE
	# Analytical
	fl3,ar3=np.genfromtxt(wd+'cdwfs_'+band+'_sens_'+str(round(logcut,1))+'_geo.dat',unpack=True)
	sens=np.interp(centers00,fl3,ar3)


	x0,y0 = np.genfromtxt(wd+'cdwfs_'+band+'_sens_'+str(round(logcut,1))+'_civ.dat',unpack=True)
	sens_easy = np.interp(centers00,x0,y0)

	a,b = np.histogram(flux_out, bins = bins00)

	c = a/(sens_easy)
	ncum_in = list(reversed(np.cumsum(list(reversed(dnds(centers00,inppars)*ds00)))))
	ncum_in_poiss = list(reversed(np.cumsum(list(reversed(quante_perarea)))))
	
	y = list(reversed(np.cumsum(list(reversed(c)))))
	y = y*(centers00/1e-14)**1.5
	y_in = ncum_in*(centers00/1e-14)**1.5
	y_in_poiss = ncum_in_poiss*(centers00/1e-14)**1.5

	w=open(wd+'lognlogs_'+band+'_'+simfolder[:-1]+'_input.dat','w')
	for jj in range(len(y_in)):
		w.write(str(centers00[jj])+' \t '+str(y_in[jj])+' \t '+str(y_in_poiss[jj])+'\n')
	w.close()

	w=open(wd+str(k)+'lognlogs_'+band+'_'+simfolder[:-1]+'_easyway.dat','w')
	for jj in range(len(y)):
		w.write(str(centers00[jj])+' \t '+str(y[jj])+'\n')
	w.close()

	# FIT (OR NOT) THE dN/dS double power law
	if fit_dnds == True:
	
		#////////////bootstrap/////////////#
		if bootstrap == True:
			# Compute dN/dS following Georgakakis+08 and use 100 bootstrap for uncertainties
			tin=time.time()
			ecfs = flux_out/np.array(crs)
			pars=[bkgs,ecfs,expos,tots]
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
				for i in range(len(expos)):	
	
					# Compute the PDFs
					T=boot[i][0]+(centers00/boot[i][1])*boot[i][2]*0.9
					prob1=scipy.stats.distributions.poisson.pmf(boot[i][3],T)*(dnds(centers00,parameters)*ds00)
	
					# Normalize them
					prob1=prob1/np.sum(prob1)
	
					# Store in the sum
					part1=part1+prob1

				# Part1 contains the sum of the PDFs
				part1b=part1/(sens)

				cumpart1b=list(reversed(np.cumsum(list(reversed(part1b)))))
				bootstrap_lognlogs.append(cumpart1b)
	
				# This is the effective dN/dS (sources/deg2/flux_bin)
				part1c=part1b/ds00
	
				bootstrap_dnds.append(part1c)

			mu_dnds,sigma_dnds=[],[]
			mu_lognlogs,sigma_lognlogs=[],[]
			for ii in range(len(bootstrap_dnds[0])):
				a,b=[],[]
				for ji in range(len(bootstrap_dnds)):
					a.append(bootstrap_dnds[ji][ii])
				mu_dnds.append(np.median(a))
				sigma_dnds.append(np.std(a))

				for ki in range(len(bootstrap_lognlogs)):
					b.append(bootstrap_lognlogs[ki][ii])
				mu_lognlogs.append(np.median(b))
				sigma_lognlogs.append(np.std(b))

			#mu=np.array(mu)
			#sigma=np.array(sigma)

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
		#////////////bootstrap/////////////#
	
		else: # No bootstrap for uncertainties
	
			print('Now fitting the dN/dS...')
			tin=time.time()
			guess = [-1,-2,5e-15]
			res = minimize(func, guess, method='nelder-mead')

			print(round(float(time.time()-tin),1),' seconds for the -likelihood minimization.')
			print('input:',inppars)
			print('guess:',guess)
			print('output:',res.x)

			pars = res.x


			if use_outflux == True:
	
				print('use_outflux is set to true; computing PDFs...')
	
				part0=np.zeros_like(centers00)
				part1=np.zeros_like(centers00)

				for i in range(len(expos)):
		
					intrinsic_flux=flux_inp[i]
					observed_flux=flux_out[i]

					ecf=observed_flux/crs[i]

					# Compute the PDFs
					T=bkgs[i]+(centers00/ecf)*expos[i]*0.9
					prob1=scipy.stats.distributions.poisson.pmf(tots[i],T)*(dnds(centers00,pars)*ds00)
	
					# Normalize them (this step should be correct)
					prob1=prob1/np.sum(prob1)

					#print(observed_flux)
					#plt.plot(centers00,prob1,'-',color='gray')

					# Store in the sum
					part1=part1+prob1
			else:
				quante_out2,bincenters_s = np.histogram(flux_inp, bins = bins00)
				part1 = quante_out2

			# Part 0/1 contains the sum of the PDFs
			part1b=part1/(sens*nsim)

			# dN/dS
			part1c=list(part1b[i]/(bins00[i+1]-bins00[i]) for i in range(len(bins00)-1))

			# logN-logS
			cumpart1=list(reversed(np.cumsum(list(reversed(part1b)))))
	else:
		print('Fit is set to false; using input parameters...')
		pars = inppars 
	
		if use_outflux == True:
	
			print('use_outflux is set to true; computing PDFs...')

			part0=np.zeros_like(centers00)
			part1=np.zeros_like(centers00)

			for i in range(len(expos)):
	
				intrinsic_flux=flux_inp[i]
				observed_flux=flux_out[i]

				ecf=observed_flux/crs[i]

				# Compute the PDFs
				T=bkgs[i]+(centers00/ecf)*expos[i]*0.9
				prob1=scipy.stats.distributions.poisson.pmf(tots[i],T)*(dnds(centers00,pars)*ds00)

				# Normalize them (this step should be correct)
				prob1=prob1/np.sum(prob1)

				# Store in the sum
				part1=part1+prob1
		else:
			quante_out2,bincenters_s = np.histogram(flux_inp, bins = bins00)
			part1 = quante_out2

		# Part 0/1 contains the sum of the PDFs
		part1b=part1/(sens*nsim)

		# dN/dS
		part1c=list(part1b[i]/(bins00[i+1]-bins00[i]) for i in range(len(bins00)-1))

		# logN-logS
		cumpart1=list(reversed(np.cumsum(list(reversed(part1b)))))


	#######################################################
	y2 = mu_lognlogs*(centers00/1e-14)**1.5
	y2_err = sigma_lognlogs*(centers00/1e-14)**1.5
	w=open(wd+str(k)+'lognlogs_'+band+'_'+simfolder[:-1]+'_bootstr.dat','w')
	for j in range(len(y2)):
		w.write(str(centers00[j])+' \t '+str(y2[j])+' \t '+str(y2_err[j])+'\n')
	w.close()


sys.exit()	

# Plot the results

plt.figure()
#plt.errorbar(x,y,yerr=ey,marker='o',color='k',linestyle='None',label='CDWFS')
#plt.plot(x,y,marker='o',color='k',linestyle='None',label='Simulation')
plt.errorbar(centers00,y2,yerr=y2_err,marker='o',color='k',linestyle='None',label='Simulation')
plt.plot(x,y_in_poiss,marker='o',color='r',linestyle='None',label='Input-Poiss')
plt.plot(x,y_in,'b-',label='Input')
plt.xscale('log')
plt.yscale('log')
if band == 'soft':
	plt.xlabel('0.5-2 keV Flux')
elif band == 'hard':
	plt.xlabel('2-10 keV Flux')

plt.ylabel('N(>S)*S^1.5')
#plt.axis([3e-17,1e-12,1,500])
plt.legend()
plt.show()

sys.exit()
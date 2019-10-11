#completeness and flux_in VS flux_out following Georgakakis+08

import numpy as np
import sys
from astropy.io import fits
import matplotlib.pyplot as plt
import time
#from ciao_contrib.region.check_fov import FOVFiles
from scipy.stats import gaussian_kde
import scipy.stats.distributions
from scipy.optimize import minimize
import matplotlib as mpl
#import numdifftools as nd

def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*xx*3600)**2 +((pointa[1]-pointb[1])*3600)**2)

def build_struct(a,b,c,d,e,f,g,h,hh):
	s=[]
	for j in range(len(a)):
		s.append([a[j],b[j],c[j],d[j],e[j],f[j],g[j],h[j],hh[j]])
	s=np.array(s)
	return s

# define some gammas and list the conversion factors from PIMMS, plus the ratio to the full band flux
gamma=np.arange(0.9,2.4,0.1)

fluxrat_s=[0.2050,0.2267,0.2500,0.2749,0.3013,0.3292,0.3584,0.3888,0.4201,0.4521,0.4846,0.5173,0.5499,0.5822,0.6138]

# now interpolate these binned functions to get more finely sampled curves
xvals = np.arange(0.9, 2.31, 0.01)

for i in range(len(xvals)):
	xvals[i]=round(xvals[i],2)
fluxinterp_s = np.interp(xvals, gamma, fluxrat_s)


wd='/Users/alberto/Desktop/XBOOTES/'

band='soft'
simfolder='sim_indep/'
nsim=10

flux_inp,flux_out=[],[]
input=[]
tots,bkgs,expos,crs=[],[],[],[]
print('Starting to match '+str(nsim)+' sims in the '+band+' band...')
t_in=time.time()
for k in range(nsim):
	
	print(k+1, end=" ")
	#take catalog of detected sources (wavdetect, full mosaic 4x4, 5e-5, only bkgmap)
	cat1=fits.open(wd+simfolder+str(k)+'cdwfs_'+band+'_sim_cat1.fits')

	ra_d=cat1[1].data['RA']
	dec_d=cat1[1].data['DEC']
	tot=cat1[1].data['TOT']
	bkg=cat1[1].data['BKG']
	exp=cat1[1].data['EXP'] # Exposure in secs
	prob=cat1[1].data['PROB']
	r90=cat1[1].data['AV_R90']
	cr=cat1[1].data['CR']
	flux_d=cat1[1].data['FLUX']
	#print(len(ra_d),'total sample')

	cut = 1e-4
	'''
	#define cut at 1% of spurious fraction
	if band=='broad':
		cut=8e-5
	elif band=='soft':
		cut=6e-4
	elif band=='hard':
		cut=4e-5
	'''
	#cut detected sample at a given probability threshold 
	ra_d=ra_d[prob<=cut]
	dec_d=dec_d[prob<=cut]
	tot=tot[prob<=cut]
	bkg=bkg[prob<=cut]
	exp=exp[prob<=cut]
	cr=cr[prob<=cut]
	r90=r90[prob<=cut]
	flux_d=flux_d[prob<=cut]
	prob=prob[prob<=cut]
	print(len(ra_d),'above threshold (',cut,')')
	
	pool=build_struct(ra_d,dec_d,r90,flux_d,prob,tot,bkg,exp,cr)
	#print(len(pool),'above rel threshold')

	### NEED TO FILTER THESE SOURCES WITH THE TOTAL FOV OF THE CDWFS, SOME OF THEM ARE OUTSIDE 
	### AND CANNOT BE MATCHED BY DEFINITION
	#take list of sources in input to simulation
	#(flux_cdwfs,ra_cdwfs,dec_cdwfs)=np.genfromtxt(wd+'poiss_rand_'+band+'_lehmer.dat',unpack=True,skip_header=1)
	#w=open(wd+'poiss_rand_'+band+'_lehmer_filtered.dat','w')
	#w.write('Flux \t RA \t DEC \n')
	#my_obs = FOVFiles('@'+wd+'fov.lis')
	#for i in range(len(ra_cdwfs)):
	#	myobs = my_obs.inside(ra_cdwfs[i], dec_cdwfs[i])
	#	if len(myobs) > 0:
	#		w.write(str(flux_cdwfs[i])+' \t '+str(ra_cdwfs[i])+' \t '+str(dec_cdwfs[i])+' \n')
	#w.close()
	#print(len(ra_cdwfs))

	#take filtered list of sources in input to simulation
	if simfolder == 'sim_indep/':
		(flux_cdwfs,ra_cdwfs,dec_cdwfs)=np.genfromtxt(wd+'poiss_rand_'+band+'_filtered_new.dat',unpack=True,skip_header=1,usecols=[0,1,2])
		'''
		if band != 'soft':
			ra_cdwfs=ra_cdwfs[flux_cdwfs>3e-16]
			dec_cdwfs=dec_cdwfs[flux_cdwfs>3e-16]
			flux_cdwfs=flux_cdwfs[flux_cdwfs>3e-16]
		else:
			ra_cdwfs=ra_cdwfs[flux_cdwfs>9e-17]
			dec_cdwfs=dec_cdwfs[flux_cdwfs>9e-17]
			flux_cdwfs=flux_cdwfs[flux_cdwfs>9e-17]
		'''
	elif simfolder == 'sim_all_FINAL/':
		(flux_cdwfs,ra_cdwfs,dec_cdwfs,gamma_cdwfs)=np.genfromtxt(wd+'poiss_rand_lehmer_filtered.dat',unpack=True,skip_header=1,usecols=[0,1,2,3])
		for m in range(len(gamma_cdwfs)):
			if band=='soft':
				flux_ratio=fluxinterp_s[xvals==gamma_cdwfs[m]]
				flux_cdwfs[m]=flux_ratio*flux_cdwfs[m]
			elif band=='hard':
				flux_ratio=1.-fluxinterp_s[xvals==gamma_cdwfs[m]]
				flux_cdwfs[m]=flux_ratio*flux_cdwfs[m]
	elif simfolder == 'sim_newgamma/':
		(flux_cdwfs,ra_cdwfs,dec_cdwfs,gamma_cdwfs)=np.genfromtxt(wd+'poiss_rand_lehmer_newgamma_filtered.dat',unpack=True,skip_header=1,usecols=[0,1,2,3])
		#if band=='soft':
		#	flux_cdwfs=0.33*flux_cdwfs
		#	# Gamma =1.8
		#	#flux_cdwfs=flux_cdwfs*4.521E-01
		#elif band=='hard':
		#	flux_cdwfs=0.67*flux_cdwfs
		for m in range(len(gamma_cdwfs)):
			if band=='soft':
				flux_ratio=fluxinterp_s[xvals==gamma_cdwfs[m]]
				flux_cdwfs[m]=flux_ratio*flux_cdwfs[m]
			elif band=='hard':
				flux_ratio=1.-fluxinterp_s[xvals==gamma_cdwfs[m]]
				flux_cdwfs[m]=flux_ratio*flux_cdwfs[m]
	elif simfolder == 'sim_all_new/':
		(flux_cdwfs,ra_cdwfs,dec_cdwfs,gamma_cdwfs)=np.genfromtxt(wd+'poiss_rand_lehmer_filtered.dat',unpack=True,skip_header=1,usecols=[0,1,2,3])
		if band=='soft':
			flux_cdwfs=0.33*flux_cdwfs
			# Gamma =1.8
			#flux_cdwfs=flux_cdwfs*4.521E-01
		elif band=='hard':
			flux_cdwfs=0.67*flux_cdwfs
			#flux_cdwfs=flux_cdwfs*5.479E-01
	
	elif simfolder == 'sim_indep_new/':
		(flux_cdwfs,ra_cdwfs,dec_cdwfs)=np.genfromtxt(wd+'poiss_rand_'+band+'_filtered_new.dat',unpack=True,skip_header=1,usecols=[0,1,2])
	
	'''
	(flux_cdwfs,ra_cdwfs,dec_cdwfs,gamma_cdwfs)=np.genfromtxt(wd+'poiss_rand_lehmer_newgamma_filtered.dat',unpack=True,skip_header=1,usecols=[0,1,2,3])
	for m in range(len(gamma_cdwfs)):
		if band=='soft':
			flux_ratio=fluxinterp_s[xvals==gamma_cdwfs[m]]
			flux_cdwfs[m]=flux_ratio*flux_cdwfs[m]
		elif band=='hard':
			flux_ratio=1.-fluxinterp_s[xvals==gamma_cdwfs[m]]
			flux_cdwfs[m]=flux_ratio*flux_cdwfs[m]
	'''
	input.append(flux_cdwfs)
	
	# Assuming one fixed Gamma=1.4
	#if band=='soft':
	#	flux_cdwfs=0.33*flux_cdwfs
	#	# Gamma =1.8
	#	#flux_cdwfs=flux_cdwfs*4.521E-01
	#elif band=='hard':
	#	flux_cdwfs=0.67*flux_cdwfs
	#	#flux_cdwfs=flux_cdwfs*5.479E-01

	# Sort them to start from the bright ones
	ra_cdwfs=ra_cdwfs[::-1]
	dec_cdwfs=dec_cdwfs[::-1]
	flux_cdwfs=flux_cdwfs[::-1]
	
	unmatched,blendings=0,0
	newpool=pool

	r90eff,deff=[],[]
	for i in range(len(ra_cdwfs)):
		input_source=[ra_cdwfs[i],dec_cdwfs[i]]
	
		found=0
		counterparts=[]
		count=0

		delta = 0.0056 #(0.028 ~100")
		ra_d_filt=newpool[:,0][(newpool[:,0]>=ra_cdwfs[i]-delta) & (newpool[:,0]<=ra_cdwfs[i]+delta) & (newpool[:,1]>=dec_cdwfs[i]-delta) & (newpool[:,1]<=dec_cdwfs[i]+delta)]
		dec_d_filt=newpool[:,1][(newpool[:,0]>=ra_cdwfs[i]-delta) & (newpool[:,0]<=ra_cdwfs[i]+delta) & (newpool[:,1]>=dec_cdwfs[i]-delta) & (newpool[:,1]<=dec_cdwfs[i]+delta)]
		flux_d_filt=newpool[:,3][(newpool[:,0]>=ra_cdwfs[i]-delta) & (newpool[:,0]<=ra_cdwfs[i]+delta) & (newpool[:,1]>=dec_cdwfs[i]-delta) & (newpool[:,1]<=dec_cdwfs[i]+delta)]
		r90_d_filt=newpool[:,2][(newpool[:,0]>=ra_cdwfs[i]-delta) & (newpool[:,0]<=ra_cdwfs[i]+delta) & (newpool[:,1]>=dec_cdwfs[i]-delta) & (newpool[:,1]<=dec_cdwfs[i]+delta)]
		prob_d_filt=newpool[:,4][(newpool[:,0]>=ra_cdwfs[i]-delta) & (newpool[:,0]<=ra_cdwfs[i]+delta) & (newpool[:,1]>=dec_cdwfs[i]-delta) & (newpool[:,1]<=dec_cdwfs[i]+delta)]
		tot_d_filt=newpool[:,5][(newpool[:,0]>=ra_cdwfs[i]-delta) & (newpool[:,0]<=ra_cdwfs[i]+delta) & (newpool[:,1]>=dec_cdwfs[i]-delta) & (newpool[:,1]<=dec_cdwfs[i]+delta)]
		bkg_d_filt=newpool[:,6][(newpool[:,0]>=ra_cdwfs[i]-delta) & (newpool[:,0]<=ra_cdwfs[i]+delta) & (newpool[:,1]>=dec_cdwfs[i]-delta) & (newpool[:,1]<=dec_cdwfs[i]+delta)]
		exp_d_filt=newpool[:,7][(newpool[:,0]>=ra_cdwfs[i]-delta) & (newpool[:,0]<=ra_cdwfs[i]+delta) & (newpool[:,1]>=dec_cdwfs[i]-delta) & (newpool[:,1]<=dec_cdwfs[i]+delta)]
		cr_d_filt=newpool[:,8][(newpool[:,0]>=ra_cdwfs[i]-delta) & (newpool[:,0]<=ra_cdwfs[i]+delta) & (newpool[:,1]>=dec_cdwfs[i]-delta) & (newpool[:,1]<=dec_cdwfs[i]+delta)]
		counterparts,probabilities,distances,mat_rad,counts,countsbkg,exposures,countrates=[],[],[],[],[],[],[],[]
		for j in range(len(ra_d_filt)):
			cdwfs_source=[ra_d_filt[j],dec_d_filt[j]]
			match_rad=1.1*r90_d_filt[j]
			d=distance(input_source,cdwfs_source)
			if d <= match_rad: #found a match
				if found==0: #it's the first
					found=1
					counterparts.append(flux_d_filt[j])
					probabilities.append(prob_d_filt[j])
					distances.append(d)
					mat_rad.append(match_rad)
					counts.append(tot_d_filt[j])
					countsbkg.append(bkg_d_filt[j])
					exposures.append(exp_d_filt[j])
					countrates.append(cr_d_filt[j])
				else: #it's not the first match	
					count=count+2
					blendings=blendings+1
					counterparts.append(flux_d_filt[j])
					probabilities.append(prob_d_filt[j])
					distances.append(d)
					mat_rad.append(match_rad)
					counts.append(tot_d_filt[j])
					countsbkg.append(bkg_d_filt[j])
					exposures.append(exp_d_filt[j])
					countrates.append(cr_d_filt[j])
		if found == 1:
			if len(counterparts) == 1:
				flux_inp.append(flux_cdwfs[i])
				flux_out.append(counterparts[0])
				r90eff.append(mat_rad[0])
				deff.append(distances[0])
				tots.append(counts[0])
				bkgs.append(countsbkg[0])
				expos.append(exposures[0])
				crs.append(countrates[0])
				#print(len(newpool[:,0]))
				#print(len(newpool[newpool[:,3]==counterparts[0]]),newpool[newpool[:,3]==counterparts[0]])
				newpool=np.delete(newpool,np.where(newpool[:,3]==counterparts[0]),0)
				#print(len(newpool[:,0]))
			else:
				counterparts=np.array(counterparts)
				probabilities=np.array(probabilities)
				mat_rad=np.array(mat_rad)
				distances=np.array(distances)
				counts=np.array(counts)
				countsbkg=np.array(countsbkg)
				exposures=np.array(exposures)
				countrates=np.array(countrates)
				
				flux_inp.append(flux_cdwfs[i])
				flux_out.append(counterparts[probabilities==np.min(probabilities)][0])
				r90eff.append(mat_rad[probabilities==np.min(probabilities)][0])
				deff.append(distances[probabilities==np.min(probabilities)][0])
				tots.append(counts[probabilities==np.min(probabilities)][0])
				bkgs.append(countsbkg[probabilities==np.min(probabilities)][0])
				expos.append(exposures[probabilities==np.min(probabilities)][0])
				crs.append(countrates[probabilities==np.min(probabilities)][0])
#print(len(newpool[newpool[:,3]==counterparts[probabilities==np.min(probabilities)]]),newpool[newpool[:,3]==counterparts[probabilities==np.min(probabilities)]])
				#print(len(newpool[:,0]))
				newpool=np.delete(newpool,np.where(newpool[:,3]==counterparts[probabilities==np.min(probabilities)]),0)

	
		else:
			unmatched=unmatched+1
	
	#print(len(ra_cdwfs), 'in input')
	#print(len(flux_inp), 'matched')
	#print(unmatched, 'unmatched')
	#print(blendings,' blendings')

	'''
	n,be=np.histogram(deff,bins=25)
	be2=list((be[i+1]+be[i])/2. for i in range(len(be)-1))
	nsum=np.cumsum(n)/float(np.sum(n))

	f,(ax1,ax2)=plt.subplots(2,1,sharex=True)
	ax1.hist(deff,bins=25)
	ax1.set_ylabel('N')
	ax2.plot(be2,nsum,'b-')
	ax2.set_ylabel('Cumulative fraction')
	ax2.set_xlabel('Distance [arcsec]')
	plt.show()
	#plt.savefig(wd+'cdwfs_dist-cum.pdf',format='pdf')
	'''

	'''
	r90eff=np.array(r90eff)
	deff=np.array(deff)

	n,bedges=np.histogram(r90eff,bins=15)
	bcenters=list((bedges[i+1]+bedges[i])/2. for i in range(len(bedges)-1))
	av=[]
	for j in range(len(bcenters)):
		flag=[]
		for i in range(len(deff)):
			if (r90eff[i] > bedges[j] and r90eff[i] < bedges[j+1]): 
				flag.append(1)
			else:
				flag.append(0)
		flag=np.array(flag)
		deff2=deff[flag==1]
		if len(deff2) > 1:
			av.append(np.median(deff2))
		else:
			av.append(deff2)

	plt.figure()
	plt.plot(r90eff,deff,'r.')
	plt.plot(bcenters,av,'ks')
	plt.xlabel('R90 [arcsec]')
	plt.ylabel('Distance [arcsec]')
	plt.tight_layout()
	#plt.savefig(wd+'cdwfs_dist-r90.pdf',format='pdf')
	plt.show()
	'''
		
t_out=time.time()
print(round(float(t_out-t_in)/60.,1),' minutes for the match.')
##### End of matching #####
###########################


bins00=np.logspace(np.log10(5e-17),np.log10(1e-12),51)
centers00=list((bins00[i+1]+bins00[i])/2. for i in range(0,len(bins00)-1))
ds00=list((bins00[i+1]-bins00[i]) for i in range(0,len(bins00)-1))
centers00=np.array(centers00)
ds00=np.array(ds00)
area=9.3 # total area of the survey 

# Flux-flux plot
'''
flux_out=np.array(flux_out)
x=np.log10(flux_inp)
y=np.log10(flux_out.astype(np.float64))
k = gaussian_kde(np.vstack([x, y]))
xi, yi = np.mgrid[x.min():x.max():x.size**0.5*1j,y.min():y.max():y.size**0.5*1j]
zi = k(np.vstack([xi.flatten(), yi.flatten()]))

#set zi to 0-1 scale
zi = (zi-zi.min())/(zi.max() - zi.min())
zi =zi.reshape(xi.shape)
#plt.imshow(zi,origin='lower')
#plt.show()

plt.figure()
#set up plot
origin = 'lower'
levels = [0.1,0.25,0.5,0.68,0.95,0.975,1]

CS = plt.contour(xi, yi, zi,levels = levels,
			  colors=('k',),
			  linewidths=(1,),
			  origin=origin)

plt.clabel(CS, fmt='%.3f', colors='b', fontsize=8)
plt.plot(np.linspace(-16,-12,20),np.linspace(-16,-12,20),'r--')
plt.scatter(x,y,color='gray',marker='.',zorder=-1,alpha=0.8)
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Input flux')
plt.ylabel('Output flux')
plt.show()
'''

# INPUT LOGNLOGS
quante,bincenters_s=np.histogram(input,bins=bins00)
quante_perarea=quante/(area*nsim)
quante_perarea_perflux=list(quante_perarea[i]/(bins00[i+1]-bins00[i]) for i in range(len(bins00)-1))
quante_perarea_perflux=np.array(quante_perarea_perflux)
ncum_in=list(reversed(np.cumsum(list(reversed(quante_perarea)))))

# OUTPUT LOGNLOGS, USING FIXED AREA AND TRUE, INPUT FLUXES
quante_out,bincenters_s=np.histogram(flux_inp,bins=bins00) # Use the "True" fluxes or the observed ones?
quante_out_perarea=quante_out/(area*nsim)
quante_out_perarea_perflux=list(quante_out_perarea[i]/(bins00[i+1]-bins00[i]) for i in range(len(bins00)-1))
quante_out_perarea_perflux=np.array(quante_out_perarea_perflux)
ncum_out=list(reversed(np.cumsum(list(reversed(quante_out_perarea)))))

# Take the sensitivity made with Georgakakis method (sens.py code)
if simfolder == 'sim_all_FINAL/':
	fl2,ar2=np.genfromtxt(wd+'cdwfs_'+band+'_sens_georgakakis_1.8.dat',unpack=True)
elif (simfolder == 'sim_indep/') or (simfolder == 'sim_newgamma/'):
	fl2,ar2=np.genfromtxt(wd+'cdwfs_'+band+'_sens_georgakakis.dat',unpack=True)
elif simfolder == 'sim_all_new/':
	fl2,ar2=np.genfromtxt(wd+'cdwfs_'+band+'_sens_georgakakis_cy18.dat',unpack=True)
	ar2_norm=ar2/np.max(ar2)

#fl3,ar3=np.genfromtxt(wd+'cdwfs_'+band+'_sens_georgakakis_r90.dat',unpack=True)
if band != 'hard':
	fl3,ar3=np.genfromtxt(wd+'cdwfs_'+band+'_interpolation_sens.dat',unpack=True)
	#fl3,ar3=np.genfromtxt(wd+'cdwfs_'+band+'_sens_georgakakis_r90.dat',unpack=True)
else:
	fl3,ar3=np.genfromtxt(wd+'cdwfs_'+band+'_sens_georgakakis_r90.dat',unpack=True)
	
# Make the sensitivity plots
ratio=np.array(ncum_out)/np.array(ncum_in)*area
num=quante_out
den=quante
ratio3=(num/den)*area
enum=np.sqrt(num)*area
eden=np.sqrt(den)*area
eratio=np.sqrt((enum/den)**2+(ratio3*eden/den)**2)

ratio3[np.isnan(ratio3)]=9.3

f,ax1=plt.subplots(1,1)
ax1.plot(centers00,ratio,'g-',linewidth=2,label='Cum(matched)/Cum(input)')
ax1.plot(centers00,ratio3,'b*',markersize=15,label='matched/input')
#ax1.errorbar(centers00,ratio3,yerr=eratio,color='blue',marker='*',ms=15,linewidth=1,label='matched/input')
#ax1.plot(fl2,ar2,'r--',linewidth=2,label='Following Georgakakis, R50')
ax1.plot(fl3,ar3,'k--',linewidth=2,label='Following Georgakakis')
ax1.set_xlabel(r'Flux (erg cm$^{-2}$ s$^{-1}$)',fontsize=15)
ax1.set_ylabel(r'Area (deg$^2$)',fontsize=15)
ax1.tick_params(axis='both',labelsize=13)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.axis([5e-17,5e-13,1e-3,10])

plt.legend()
plt.tight_layout()
plt.show()

'''
if simfolder == 'sim_indep/':
	w=open(wd+'cdwfs_'+band+'_sens_sim_indep.dat','w')
elif simfolder == 'sim_all_FINAL/':
	w=open(wd+'cdwfs_'+band+'_sens_sim_all_FINAL.dat','w')
elif simfolder == 'sim_newgamma/':
	w=open(wd+'cdwfs_'+band+'_sens_sim_newgamma.dat','w')
elif simfolder == 'sim_all_new/':
	w=open(wd+'cdwfs_'+band+'_sens_sim_all_new.dat','w')
'''
#w=open(wd+'cdwfs_'+band+'_sens_sim_indep.dat','w')
#for i in range(len(ratio3)):
#	w.write(str(centers00[i])+' \t '+str(ratio3[i])+'\n')
#w.close()

#######################################################
# OUTPUT LOGNLOGS, CORRECTING FOR SENSITIVITY

# Use the Georgakakis sens, but need to interpolate - standard method
sens=np.interp(centers00,fl3,ar3)

quante_out_perarea_corr=quante_out/(sens*nsim)
quante_out_perarea_perflux_corr=list(quante_out_perarea_corr[i]/(bins00[i+1]-bins00[i]) for i in range(len(bins00)-1))
ncum_out_corr=list(reversed(np.cumsum(list(reversed(quante_out_perarea_corr)))))

# Use the Georgakakis sens AND methodology to correct for flux uncertainty (and Eddington
# bias)
# Assume I know beta (b1=-1.34, b2=-2.35) for the full band
#b1=-1.34
#b2=-2.35
#fb=8.1e-15

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

'''
########
def func(params):
	fx = centers00
	k,b1,b2,fb=params[0],params[1],params[2],params[3]
	fref=1e-14	
	k1 = k*(fb/fref)**(b1-b2)
	
	if type(fx) == list:
		fx = np.array(fx)
	aux = k*(fx/fref)**b1
	aux[fx > fb] = k1*(fx[fx > fb]/fref)**b2
	return aux


#p0=np.logspace(np.log10(4e11),np.log10(4e20),40)
p0=562e14
p1=np.linspace(-2,-1,11)
#p1=-1.34
#p2=np.linspace(-3,-2,40)
p2=-2.35
#p3=8.1e-15
p3=np.logspace(np.log10(1e-15),np.log10(1e-13),11)
cmap = mpl.cm.autumn
plt.figure()
for h in range(len(p3)):
	res=[]
	for k in range(len(p1)):
	
		guess=[p0,p1[k],p2,p3[h]]
		lnpconi=[]
		for i in range(len(expos)):
			ecf=flux_out[i]/crs[i]
			T=bkgs[i]+(centers00/ecf)*expos[i]*0.9
		
			pb=scipy.stats.distributions.poisson.pmf(tots[i],T)*func(guess)*ds00
			pconi = np.sum(pb)/np.sum(func(guess)*sens*ds00)
		
			lnpconi.append(-np.log(pconi))
	
		lnpconi=np.array(lnpconi)
		res.append(np.sum(lnpconi))


	plt.plot(p1,res,linestyle='-',color=cmap(h / float(len(p3))),label=str(p3[h]))

plt.axvline(x=-1.34)
#plt.xscale('log')
plt.legend()
plt.show()

sys.exit()
########
'''

tin=time.time()
guess=[-1.3,-2.5,6e-15]
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

res=minimize(func, guess, method='nelder-mead')

if band == 'broad':
	inppars=[-1.34,-2.35,8.1e-15] # broad Lehmer's dN/dS params
elif band == 'soft':
	inppars=[-1.49,-2.48,6.0e-15] # soft
else:
	inppars=[-1.32,-2.55,6.4e-15] # hard

print(round(float(time.time()-tin)/60.,1),' minutes for the -likelihood minimization.')
print('input:',inppars)
print('guess:',guess)
print('output:',res.x)

#H=nd.Hessian(func)(res.x)
#print(H)

pars=res.x
part0=np.zeros_like(centers00)
part1=np.zeros_like(centers00)

for i in range(len(expos)):
	if (expos[i] > 0) and (expos[i] < 4e7):
		
		intrinsic_flux=flux_inp[i]
		observed_flux=flux_out[i]
		
		ecf=flux_out[i]/crs[i]
		# Compute the PDFs
		T=bkgs[i]+(centers00/ecf)*expos[i]*0.9
		prob1=scipy.stats.distributions.poisson.pmf(tots[i],T)*(dnds(centers00,pars)*ds00)
		
		# Normalize them (this step should be correct)
		prob1=prob1/np.sum(prob1)
	
		# Store in the sum
		part1=part1+prob1

print('I cycled on',len(expos),'sources, and the total number of sources in the sum is',round(np.sum(part1),0))

# Part 0/1 contains the sum of the PDFs
part1b=part1/(sens*nsim)

#part0c=list(part0b[i]/(bins00[i+1]-bins00[i]) for i in range(len(bins00)-1))
part1c=list(part1b[i]/(bins00[i+1]-bins00[i]) for i in range(len(bins00)-1))

# If I first cumulate the PDFs and THEN divide for the CUMULATED sensitivity
cumsens=ratio # This is the "old" flatter curve
cumsens[-1]=1.
cumsens=np.array(cumsens)

#cumpart0=list(reversed(np.cumsum(list(reversed(part0)))))
#cumpart0=np.array(cumpart0)
#cumpart0=cumpart0/(cumsens*nsim)

#cumpart1=list(reversed(np.cumsum(list(reversed(part1)))))
cumpart1=list(reversed(np.cumsum(list(reversed(part1b)))))

#######################################################

'''
print('dN/dS with Georgakakis method, no correction fx^beta:')
print(part0c)
print('dN/dS in input:')
print(quante_perarea_perflux)
print('Their ratio is:')
RATIO=quante_perarea_perflux/part0c
for i in range(len(RATIO)):
	print(centers00[i],RATIO[i])

plt.figure()
plt.plot(centers00,RATIO,'ro')
plt.axhline(y=1)
plt.axvline(x=1e-14)
plt.xscale('log')
plt.yscale('log')
plt.show()
'''

f,(ax1,ax2)=plt.subplots(2,1,sharex=True,figsize=[7,9])
#ax1.plot(centers00,ncum_out,'r-',linewidth=2,label='Using input fluxes, fixed area')
#ax1.plot(centers00,ncum_out_corr,'b-',linewidth=2,label='Using input fluxes, correcting for sens')
#ax1.plot(centers00,cumpart0,'g-',linewidth=2,label=r'Georgakakis method, NO fx$^{\beta}$')
ax1.plot(centers00,cumpart1,'cs',linewidth=2,label=r'Georgakakis method, with fx$^{\beta}$')
ax1.plot(centers00,ncum_in,'k--',linewidth=3, label='In input file')
ax1.set_ylabel(r'N(>S) (deg$^{-2}$)',fontsize=13)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.axis([5e-17,2e-12,1,2e4])
ax1.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=12)

#ax2.plot(centers00,quante_out_perarea_perflux,'r*',linewidth=2,label='Using input fluxes, fixed area')
#ax2.plot(centers00,quante_out_perarea_perflux_corr,'bs',linewidth=2,label='Using input fluxes, correcting for sens')
#ax2.plot(centers00,part0c,'gs',linewidth=2,label=r'Georgakakis method, NO fx$^{\beta}$')
ax2.plot(centers00,part1c,'cs',linewidth=2,label=r'Georgakakis method, with fx$^{\beta}$')
ax2.plot(centers00,quante_perarea_perflux,'ko',linewidth=2, label='In input file')
#ax2.plot(np.logspace(-16,-14,20),(np.logspace(-16,-14,20)/1e-14)**(-1.34)*560e14,'g--')
#ax2.plot(np.logspace(-14,-12,20),(np.logspace(-14,-12,20)/1e-14)**(-2.35)*560e14,'c--')
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
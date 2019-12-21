# Trying to get a spurious fraction which is independent of the number of input sources
# (i.e., logn-logs) but depends only on how background and exposure behave

import numpy as np
import sys
from astropy.io import fits
import matplotlib.pyplot as plt
import time
#from ciao_contrib.region.check_fov import FOVFiles
from scipy.stats import gaussian_kde
import scipy.stats.distributions

def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*xx*3600)**2 +((pointa[1]-pointb[1])*3600)**2)

def build_struct(a,b,c,d,e):
	s=[]
	for j in range(len(a)):
		s.append([a[j],b[j],c[j],d[j],e[j]])
	s=np.array(s)
	return s

def choose_best(array):
	if len(array) == 1:
		return array
	else:
		# Assuming prob is element [2]
		probs = list(zip(*array))[2]
		for i in range(len(array)):
			if array[i][2] == np.min(probs):
				return array[i]
					
wd='/Users/alberto/Desktop/XBOOTES/'

band='hard'
min=-6
max=0
nbins=(max-min)*20+1
bins=np.linspace(min,max,nbins)
'''
cc=0
for simfolder in ['sim_indep_22-Nov-19/', 'sim_indep_06-Dec-19/']:
	
	inp_src_filename='poiss_rand_'+band+'_'+simfolder[:-1]+'_filtered_new.dat'
	(flux_sim,ra_sim,dec_sim)=np.genfromtxt(wd+inp_src_filename,unpack=True,skip_header=1,usecols=[0,1,2])
	if simfolder == 'sim_indep_22-Nov-19/':
		nsim = 10
		sp_prob=[]
		cc=cc+nsim
	elif simfolder == 'sim_indep_06-Dec-19/':
		nsim=10
		sp_prob2=[]
		cc=cc+nsim
	elif simfolder == 'sim_indep_12-Dec-19/':
		nsim=10
		sp_prob3=[]
		cc=cc+nsim
		
	spurious, real=[],[]
	e_spurious, e_real=[],[]
	input=[]
	for logcut in bins:

		print(logcut)
		cut=10**logcut
		
		t_in=time.time()
		##########################
		##### Start matching #####
		sp,re = [],[]
		c=1
		for k in range(nsim):
			print(c,'/',cc)
	
			# Take catalog of detected sources (wavdetect, full mosaic 4x4, 5e-5, only bkgmap)
			cat1 = fits.open(wd+simfolder+str(k)+'cdwfs_'+band+'_sim_cat1_exp-psf.fits')
	
			simdata = cat1[1].data
	
			prob = simdata['PROB']
			simdata = np.array(simdata)
			
			# Cut detected sample at a given probability threshold 
			simdata = simdata[prob <= cut]
			tot = len(simdata)
			print(tot,'above threshold ('+str(round(logcut,1))+')')
	
			# Input sources depend on band and simfolder
			(flux_cdwfs,ra_cdwfs,dec_cdwfs)=np.genfromtxt(wd+inp_src_filename,unpack=True,skip_header=1,usecols=[0,1,2])	
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
		
				counterparts,distances,mat_rad,counts,countsbkg,exposures,countrates=[],[],[],[],[],[],[]
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
				
						simdata=np.delete(simdata,np.where(simdata['RA']==counterparts[0][0]),0)
								
					else:
						counterparts=np.array(counterparts)
						counterparts=choose_best(counterparts)
						
						simdata=np.delete(simdata,np.where(simdata['RA']==counterparts[0]),0)
				
				else:
					unmatched=unmatched+1
	
			sp.append(len(simdata))
			re.append(tot-len(simdata))
			print(len(simdata),'spurious sources left.')
			c=c+1
			
		real.append(np.median(re))
		e_real.append(np.std(re))
		spurious.append(np.median(sp))
		e_spurious.append(np.std(sp))
		t_out=time.time()
		print(round(float(t_out-t_in),1),' seconds for the match.')
		##### End of matching #####
		###########################	
		
		
	t_out=time.time()
	print(spurious)

	w = open(wd+'spurious_sources_'+band+'_'+simfolder[:-1]+'_finer.dat','w')
	for j in range(len(bins)):
		#w.write(str(round(bins[j],1))+' \t '+str(spurious[j])+' \t '+str(e_spurious[j])+' \t '+str(real[j])+' \t '+str(e_real[j])+'\n')
		w.write(str(bins[j])+' \t '+str(spurious[j])+' \t '+str(e_spurious[j])+' \t '+str(real[j])+' \t '+str(e_real[j])+'\n')
	w.close()

sys.exit()
'''

band = ['broad','soft','hard']
xsize = 6
f,ax=plt.subplots(len(band), 1, figsize=[xsize,xsize*len(band)], sharex=True, sharey=True)
for i in range(len(band)):

	bins,spurious,e_spurious,real,e_real=np.genfromtxt(wd+'spurious_sources_'+band[i]+'_sim_indep_22-Nov-19_finer.dat', unpack=True)
	
	bins2,spurious2,e_spurious2,real2,e_real2=np.genfromtxt(wd+'spurious_sources_'+band[i]+'_sim_indep_06-Dec-19_finer.dat', unpack=True)
	
	bins3,spurious3,e_spurious3,real3,e_real3=np.genfromtxt(wd+'spurious_sources_'+band[i]+'_sim_indep_12-Dec-19_finer.dat', unpack=True)
	
	print(np.max(spurious))
	norm_sp = spurious/np.max(spurious)
	e_norm_sp = e_spurious/np.max(spurious)
	norm_re = real/np.max(real)
	e_norm_re = e_real/np.max(real)

	norm_sp2 = spurious2/np.max(spurious2)
	e_norm_sp2 = e_spurious2/np.max(spurious2)
	norm_re2 = real2/np.max(real2)
	e_norm_re2 = e_real2/np.max(real2)
	
	norm_sp3 = spurious3/np.max(spurious3)
	e_norm_sp3 = e_spurious3/np.max(spurious3)
	norm_re3 = real3/np.max(real3)
	e_norm_re3 = e_real3/np.max(real3)
	
	diff = norm_re-norm_sp
	e_diff = np.sqrt(e_norm_re**2+e_norm_sp**2)

	diff2 = norm_re2-norm_sp2
	e_diff2 = np.sqrt(e_norm_re2**2+e_norm_sp2**2)
	
	diff3 = norm_re3-norm_sp3
	e_diff3 = np.sqrt(e_norm_re3**2+e_norm_sp3**2)
	
	peak_prob = bins[diff==np.max(diff)]
	
	peak_prob2 = bins2[diff2==np.max(diff2)]
	
	peak_prob3 = bins3[diff3==np.max(diff3)]
	
	print(peak_prob)
	print(spurious[bins==peak_prob])
	
	print(peak_prob2)
	print(spurious2[bins2==peak_prob2])
	
	print(peak_prob3)
	print(spurious3[bins3==peak_prob3])
	
	xmin = (peak_prob+8.)/8.
	if len(band) > 1:
	
		ax[i].errorbar(bins,norm_sp,yerr=e_norm_sp,color='k')
		ax[i].errorbar(bins,norm_re,yerr=e_norm_re,color='k')
		ax[i].errorbar(bins,diff,yerr=e_diff,color='b')
		
		ax[i].errorbar(bins2,norm_sp2,yerr=e_norm_sp2,color='k',linestyle='dashed')
		ax[i].errorbar(bins2,norm_re2,yerr=e_norm_re2,color='k',linestyle='dashed')
		ax[i].errorbar(bins2,diff2,yerr=e_diff2,color='b',linestyle='dashed')
		
		ax[i].errorbar(bins3,norm_sp3,yerr=e_norm_sp3,color='k',linestyle='dotted')
		ax[i].errorbar(bins3,norm_re3,yerr=e_norm_re3,color='k',linestyle='dotted')
		ax[i].errorbar(bins3,diff3,yerr=e_diff3,color='b',linestyle='dotted')
		
		ax[i].axvline(x=peak_prob,color='red',linestyle='-')
		ax[i].axvline(x=peak_prob2,color='red',linestyle='--')
		ax[i].axvline(x=peak_prob3,color='red',linestyle='dotted')
		ax[i].axhline(y=norm_sp[bins==peak_prob],xmin=xmin, color='red',linestyle='--')
		ax[i].tick_params(top=True, direction='in')
		ax[i].set_ylabel('Fraction')
		ax[i].annotate(band[i].capitalize(), xy=(-7.5,0.9))
		ax[i].annotate('"Real"', xy=(-3.7,0.8))
		ax[i].annotate('Spurious', xy=(-2.2,0.6))
		ax[i].annotate('Difference', xy=(-2.0,0.25), color='blue')
		ax[i].set_xlim([-8,0])
	
	
		locs, labels = plt.yticks() 
	
		labels = (locs*np.max(spurious)).astype(int)
		ax2=ax[i].twinx()
		lim1 = ax[i].get_ylim()
		ax2.set_ylim(lim1)
		ax2.set_ylabel('Number of spurious')
		ax2.set_yticklabels(labels)
		ax2.tick_params(direction='in')
	else:
		
		ax.errorbar(bins,norm_sp,yerr=e_norm_sp,color='k')
		ax.errorbar(bins,norm_re,yerr=e_norm_re,color='k')
		ax.errorbar(bins,diff,yerr=e_diff,color='b')
		
		ax.errorbar(bins2,norm_sp2,yerr=e_norm_sp2,color='k',linestyle='dashed')
		ax.errorbar(bins2,norm_re2,yerr=e_norm_re2,color='k',linestyle='dashed')
		ax.errorbar(bins2,diff2,yerr=e_diff2,color='b',linestyle='dashed')
		
		ax.errorbar(bins3,norm_sp3,yerr=e_norm_sp3,color='k',linestyle='dotted')
		ax.errorbar(bins3,norm_re3,yerr=e_norm_re3,color='k',linestyle='dotted')
		ax.errorbar(bins3,diff3,yerr=e_diff3,color='b',linestyle='dotted')
		
		ax.axvline(x=peak_prob,color='red',linestyle='-')
		ax.axvline(x=peak_prob2,color='red',linestyle='--')
		#ax.axhline(y=norm_sp[bins==peak_prob],xmin=xmin, color='red',linestyle='--')
		ax.tick_params(top=True, direction='in')
		ax.set_ylabel('Fraction')
		#ax.annotate(band[i].capitalize(), xy=(-7.5,0.9))
		#ax.annotate('"Real"', xy=(-3.7,0.8))
		#ax.annotate('Spurious', xy=(-2.2,0.6))
		#ax.annotate('Difference', xy=(-2.0,0.25), color='blue')
		ax.set_xlim([-8,0])
		ax.set_ylim([-0.1,1.1])
		'''
		locs, labels = plt.yticks() 
	
		labels = (locs*np.max(spurious)).astype(int)
		ax2=ax.twinx()
		lim1 = ax.get_ylim()
		ax2.set_ylim(lim1)
		ax2.set_ylabel('Number of spurious')
		ax2.set_yticklabels(labels)
		ax2.tick_params(direction='in')
		'''
if len(band) > 1:
	ax[i].set_xlabel(r'$\log{P}$')
	plt.subplots_adjust(hspace=0)
else:
	ax.set_xlabel(r'$\log{P}$')
	plt.tight_layout()
plt.show()
#plt.savefig('/Users/alberto/Desktop/3.png')

sys.exit()


######### OLD VERSION WITH ALL THE KINDS OF SIMULATIONS TRIED
'''
band='soft'
cut=1

t_in=time.time()
spurious=[]
c=1
cc=0
for simfolder in ['sim_indep_28-Oct-19/','sim_all_new/','sim_newgamma/','sim_indep/']:
	if simfolder == 'sim_indep_28-Oct-19/':
		nsim = 1
		(flux_sim,ra_sim,dec_sim)=np.genfromtxt(wd+'poiss_rand_'+band+'_filtered_new.dat',unpack=True,skip_header=1)
		sp_prob=[]
		cc=cc+nsim
	elif simfolder == 'sim_all_new/':
		nsim = 10
		(flux_sim,ra_sim,dec_sim,gamma_sim)=np.genfromtxt(wd+'poiss_rand_lehmer_filtered.dat',unpack=True,skip_header=1)
		sp_prob0=[]
		cc=cc+nsim
	elif simfolder == 'sim_newgamma/':
		nsim = 10
		(flux_sim,ra_sim,dec_sim,gamma_sim)=np.genfromtxt(wd+'poiss_rand_lehmer_newgamma_filtered.dat',unpack=True,skip_header=1)
		sp_prob1=[]
		cc=cc+nsim
	elif simfolder == 'sim_indep/':
		nsim=10
		(flux_sim,ra_sim,dec_sim)=np.genfromtxt(wd+'poiss_rand_'+band+'_filtered_new.dat',unpack=True,skip_header=1)
		sp_prob2=[]
		cc=cc+nsim
		
	match=[]
	for k in range(nsim):
		print(c,'/',cc)
		#take catalog of detected sources (wavdetect, full mosaic 4x4, 5e-5, only bkgmap) cleaned by multiple sources (cat1)
		if (simfolder != 'sim_indep/') and (simfolder != 'sim_indep_28-Oct-19/'):
			cat1=fits.open(wd+simfolder+str(k)+'cdwfs_'+band+'_sim_cat1.fits')
		else:
			if band != 'hard':
				cat1=fits.open(wd+simfolder+str(k)+'cdwfs_'+band+'_sim_cat1_exp-psf.fits')
			else:
				cat1=fits.open(wd+simfolder+str(k)+'cdwfs_'+band+'_sim_cat1_exp-psf.fits')
		detected=[]
		
	
		ra_det=cat1[1].data['RA']
		dec_det=cat1[1].data['DEC']
		ctsf_det=cat1[1].data['TOT']
		prob_det=cat1[1].data['PROB']
		r90_det=cat1[1].data['AV_R90']
		flux_det=cat1[1].data['FLUX']
	
		ra_det=ra_det[prob_det<=cut]
		dec_det=dec_det[prob_det<=cut]
		ctsf_det=ctsf_det[prob_det<=cut]
		r90_det=r90_det[prob_det<=cut]
		flux_det=flux_det[prob_det<=cut]
		prob_det=prob_det[prob_det<=cut]
	
		detected.append(len(ra_det))
		pool=build_struct(ra_det,dec_det,r90_det,flux_det,prob_det)
	
		#print('Starting to match...')

		
		flux_inp,flux_out=[],[]
		unmatched,blendings=0,0
		newpool=pool
		for i in range(len(ra_sim)):
			input_source=[ra_sim[i],dec_sim[i]]
	
			found=0
			counterparts=[]
	

			delta = 0.007 #(0.028 ~100")
			ra_d_filt=newpool[:,0][(newpool[:,0]>=ra_sim[i]-delta) & (newpool[:,0]<=ra_sim[i]+delta) & (newpool[:,1]>=dec_sim[i]-delta) & (newpool[:,1]<=dec_sim[i]+delta)]
			dec_d_filt=newpool[:,1][(newpool[:,0]>=ra_sim[i]-delta) & (newpool[:,0]<=ra_sim[i]+delta) & (newpool[:,1]>=dec_sim[i]-delta) & (newpool[:,1]<=dec_sim[i]+delta)]
			flux_d_filt=newpool[:,3][(newpool[:,0]>=ra_sim[i]-delta) & (newpool[:,0]<=ra_sim[i]+delta) & (newpool[:,1]>=dec_sim[i]-delta) & (newpool[:,1]<=dec_sim[i]+delta)]
			r90_d_filt=newpool[:,2][(newpool[:,0]>=ra_sim[i]-delta) & (newpool[:,0]<=ra_sim[i]+delta) & (newpool[:,1]>=dec_sim[i]-delta) & (newpool[:,1]<=dec_sim[i]+delta)]
			prob_d_filt=newpool[:,4][(newpool[:,0]>=ra_sim[i]-delta) & (newpool[:,0]<=ra_sim[i]+delta) & (newpool[:,1]>=dec_sim[i]-delta) & (newpool[:,1]<=dec_sim[i]+delta)]
			counterparts,probabilities,distances,mat_rad=[],[],[],[]
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
					else: #it's not the first match	
						blendings=blendings+1
						counterparts.append(flux_d_filt[j])
						probabilities.append(prob_d_filt[j])
						distances.append(d)
						mat_rad.append(match_rad)
			if found == 1:
				if len(counterparts) == 1:
					flux_inp.append(flux_sim[i])
					flux_out.append(counterparts[0])
					
					match.append(probabilities[0])
					
					newpool=np.delete(newpool,np.where(newpool[:,3]==counterparts[0]),0)

				else:
					counterparts=np.array(counterparts)
					probabilities=np.array(probabilities)
					mat_rad=np.array(mat_rad)
					distances=np.array(distances)
				
					match.append(np.min(probabilities))
					
					flux_inp.append(flux_sim[i])
					flux_out.append(counterparts[probabilities==np.min(probabilities)])
				
					newpool=np.delete(newpool,np.where(newpool[:,3]==counterparts[probabilities==np.min(probabilities)]),0)

	
			else:
				unmatched=unmatched+1
	
		if simfolder == 'sim_indep_28-Oct-19/':
			matched=match
			for kk in range(len(newpool[:,4])):
				sp_prob.append(newpool[:,4][kk])	
		elif simfolder == 'sim_all_new/':
			matched0=match
			for kk in range(len(newpool[:,4])):
				sp_prob0.append(newpool[:,4][kk])
		elif simfolder == 'sim_newgamma/':
			matched1=match
			for kk in range(len(newpool[:,4])):
				sp_prob1.append(newpool[:,4][kk])
		elif simfolder == 'sim_indep/':
			matched2=match
			for kk in range(len(newpool[:,4])):
				sp_prob2.append(newpool[:,4][kk])
		
		
		#print('*'*20)
		#print('P cut:',cut)
		#print(len(ra_sim), 'in input of simulation.')
		#print(len(ra_det),'detected with P < P_cut.')
		#print(len(flux_inp), 'matched.')
		#print(len(newpool), 'spurious.')
		#print(unmatched, 'unmatched from input.')
		#print(blendings,' blendings.')
		#print(float(t_out-t_in),' seconds for the match.')
		#print('*'*20)
		
		spurious.append(len(newpool))
		c=c+1
		
t_out=time.time()
print(float(t_out-t_in),' seconds for the match.')
print(spurious)
print(np.median(spurious),np.mean(spurious))
plt.figure()
plt.hist(spurious,bins=10)
plt.show()

min=-15
max=0
nbins=(max-min)*10+1
bins=np.linspace(min,max,nbins)

matched=np.array(matched)	
logmat=np.log10(matched)

matched0=np.array(matched0)	
logmat0=np.log10(matched0)

matched1=np.array(matched1)	
logmat1=np.log10(matched1)

matched2=np.array(matched2)	
logmat2=np.log10(matched2)


m,be=np.histogram(logmat,bins=bins)
bc=list((be[i+1]+be[i])/2.for i in range(len(be)-1))
bc=np.array(bc)
mcum=np.cumsum(m)/np.sum(m)

m0,be=np.histogram(logmat0,bins=bins)
mcum0=np.cumsum(m0)/np.sum(m0)

m1,be=np.histogram(logmat1,bins=bins)
mcum1=np.cumsum(m1)/np.sum(m1)

m2,be=np.histogram(logmat2,bins=bins)
mcum2=np.cumsum(m2)/np.sum(m2)

#print(np.sum(m[bc<np.log10(pthresh)])/np.sum(m))
#print(np.sum(m0[bc<np.log10(pthresh)])/np.sum(m0))
#print(np.sum(m1[bc<np.log10(pthresh)])/np.sum(m1))
#print(np.sum(m2[bc<np.log10(pthresh)])/np.sum(m2))

sp_prob=np.array(sp_prob)	
logsp=np.log10(sp_prob)

sp_prob0=np.array(sp_prob0)	
logsp0=np.log10(sp_prob0)

sp_prob1=np.array(sp_prob1)	
logsp1=np.log10(sp_prob1)

sp_prob2=np.array(sp_prob2)	
logsp2=np.log10(sp_prob2)


n,be=np.histogram(logsp,bins=bins)
ncum=np.cumsum(n)/np.sum(n)

n0,be=np.histogram(logsp0,bins=bins)
ncum0=np.cumsum(n0)/np.sum(n0)

n1,be=np.histogram(logsp1,bins=bins)
ncum1=np.cumsum(n1)/np.sum(n1)

n2,be=np.histogram(logsp2,bins=bins)
ncum2=np.cumsum(n2)/np.sum(n2)

#print(np.sum(n[bc<np.log10(pthresh)])/np.sum(n))
#print(np.sum(n0[bc<np.log10(pthresh)])/np.sum(n0))
#print(np.sum(n1[bc<np.log10(pthresh)])/np.sum(n1))
#print(np.sum(n2[bc<np.log10(pthresh)])/np.sum(n2))


diff=mcum-ncum
diff0=mcum0-ncum0
diff1=mcum1-ncum1
diff2=mcum2-ncum2

sum=mcum+ncum
sum0=mcum0+ncum0
sum1=mcum1+ncum1
sum2=mcum2+ncum2


# Write out output
w=open(wd+'spurious_sources_'+band+'.dat','w')
w.write('# LogP 	 Cumfrac_mat_sim_all_indep_28-Oct-19(yes ECF, <Gamma>=1.4, different logn-logs for each band)) 	 Cumfrac_mat_sim_all_new(yes ECF, <Gamma>=1.8) 	 Cumfrac_mat_sim_newgamma(yes ECF, <Gamma>=1.5) 	 Cumfrac_mat_sim_indep(yes ECF, Gamma=1.4, different logn-logs for each band) 	 Cumfrac_sp_sim_all_FINAL 	 Cumfrac_sp_sim_all_new 	 Cumfrac_sp_sim_newgamma 	 Cumfrac_sp_sim_indep 	 Diff_sim_all_FINAL 	 Diff_sim_all_new 	 Diff_sim_newgamma 	 Diff_sim_indep 	 Sum_sim_all_FINAL 	 Sum_sim_all_new 	 Sum_sim_newgamma 	 Sum_sim_indep\n')
for u in range(len(bc)):
	w.write(str(bc[u])+' \t '+str(mcum[u])+' \t '+str(mcum0[u])+' \t '+str(mcum1[u])+' \t '+str(mcum2[u])+' \t '+str(ncum[u])+' \t '+str(ncum0[u])+' \t '+str(ncum1[u])+' \t '+str(ncum2[u])+' \t '+str(diff[u])+' \t '+str(diff0[u])+' \t '+str(diff1[u])+' \t '+str(diff2[u])+' \t '+str(sum[u])+' \t '+str(sum0[u])+' \t '+str(sum1[u])+' \t '+str(sum2[u])+'\n')
w.close()

sys.exit()
'''

f,ax=plt.subplots(2,3,sharex=True, figsize=[18,5], gridspec_kw={'height_ratios': [3, 1]})

band=['broad','soft','hard']
#band=['soft']
for i in range(len(band)):
	if band[i] == 'broad':
		#pthresh = 8e-5
		pthresh = 1e-4
	elif band[i] == 'soft':
		#pthresh = 6e-4
		pthresh = 1e-4
	else:
		#pthresh = 4e-5
		pthresh = 1e-4

	bc,mcum,mcum0,mcum1,mcum2,ncum,ncum0,ncum1,ncum2,diff,diff0,diff1,diff2,sum,sum0,sum1,sum2=np.genfromtxt(wd+'spurious_sources_'+band[i]+'.dat',unpack=True,skip_header=1)
	
	print(bc[diff==np.max(diff)])
	print(bc[diff0==np.max(diff0)])
	print(bc[diff1==np.max(diff1)])
	print(bc[diff2==np.max(diff2)])
	
	print('*'*15)
	
	ax[0][i].plot(bc,mcum,'r--')
	ax[0][i].plot(bc,mcum0,'g--')
	ax[0][i].plot(bc,mcum1,'b--')
	ax[0][i].plot(bc,mcum2,'k--',label='Real sources')
	

	ax[0][i].plot(bc,ncum,'r-')
	ax[0][i].plot(bc,ncum0,'g-')
	ax[0][i].plot(bc,ncum1,'b-')
	ax[0][i].plot(bc,ncum2,'k-',label='Spurious sources')
	
	
	ax[0][i].legend()
	ax[0][0].set_ylabel('Cumulative fraction')
	ax[0][i].tick_params(axis='both',direction='in',top=True,right=True)
	ax[0][i].axvline(x=np.log10(pthresh),color='k')
	ax[0][i].axvline(x=bc[diff==np.max(diff)],color='gray',linestyle='dashed')
	ax[0][i].axvline(x=bc[diff0==np.max(diff0)],color='gray',linestyle='dashed')
	ax[0][i].axvline(x=bc[diff1==np.max(diff1)],color='gray',linestyle='dashed')
	ax[0][i].axvline(x=bc[diff2==np.max(diff2)],color='gray',linestyle='dashed')
	ax[0][i].text(x=-14.5,y=0.7,s=band[i].capitalize())
	if i != 0:
		ax[0][i].set_yticklabels([])
	#ax[0][i].axvline(x=-4.3,color='red')
	
	ax[1][i].plot(bc,diff,'r-')
	ax[1][i].plot(bc,diff0,'g-')
	ax[1][i].plot(bc,diff1,'b-')
	ax[1][i].plot(bc,diff2,'k-')

	
	ax[1][0].set_ylabel('Difference')
	ax[1][i].tick_params(axis='both',direction='in',top=True,right=True)
	ax[1][i].axvline(x=np.log10(pthresh),color='k')
	ax[1][i].axvline(x=bc[diff==np.max(diff)],color='gray',linestyle='dashed')
	ax[1][i].axvline(x=bc[diff0==np.max(diff0)],color='gray',linestyle='dashed')
	ax[1][i].axvline(x=bc[diff1==np.max(diff1)],color='gray',linestyle='dashed')
	ax[1][i].axvline(x=bc[diff2==np.max(diff2)],color='gray',linestyle='dashed')
	if i != 0:
		ax[1][i].set_yticklabels([])
	#ax[1][i].axvline(x=-4.3,color='red')
	ax[1][i].set_xlabel('Log(P) of being spurious')
	
	#ax[2][i].plot(bc,sum,'r-')
	#ax[2][i].plot(bc,sum0,'g-')
	#ax[2][i].plot(bc,sum1,'b-')
	#ax[2][i].plot(bc,sum2,'k-')

	
	#ax[2][i].axvline(x=np.log10(pthresh),color='k')
	#ax[2][i].axvline(x=bc[diff==np.max(diff)],color='gray',linestyle='dashed')
	#ax[2][i].axvline(x=bc[diff0==np.max(diff0)],color='gray',linestyle='dashed')
	#ax[2][i].axvline(x=bc[diff1==np.max(diff1)],color='gray',linestyle='dashed')
	#ax[2][i].axvline(x=bc[diff2==np.max(diff2)],color='gray',linestyle='dashed')
	#ax[2][i].tick_params(axis='both',direction='in',top=True,right=True)
	#ax[2][i].set_xlabel('Log(P) of being spurious')
	#ax[2][i].set_ylabel('Sum')
	#ax[2][i].axvline(x=-4.3,color='red')
	
plt.subplots_adjust(hspace=0,wspace=0)
#plt.tight_layout()
plt.show()

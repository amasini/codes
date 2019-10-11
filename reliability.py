# reliability
import numpy as np
import sys
from astropy.io import fits
import matplotlib.pyplot as plt
import time
from ciao_contrib.region.check_fov import FOVFiles

def build_struct(a,b,c,d,e):
	s=[]
	for j in range(len(a)):
		s.append([a[j],b[j],c[j],d[j],e[j]])
	s=np.array(s)
	return s
	
def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*xx*3600)**2 +((pointa[1]-pointb[1])*3600)**2)

wd='/Users/alberto/Desktop/XBOOTES/'

simfolder = 'sim_indep/'

band=str(sys.argv[1])

###########################

# define some gammas and list the conversion factors from PIMMS, plus the ratio to the full band flux
gamma=np.arange(0.9,2.4,0.1)

fluxrat_s=[0.2050,0.2267,0.2500,0.2749,0.3013,0.3292,0.3584,0.3888,0.4201,0.4521,0.4846,0.5173,0.5499,0.5822,0.6138]

# now interpolate these binned functions to get more finely sampled curves
xvals = np.arange(0.9, 2.31, 0.01)

for i in range(len(xvals)):
	xvals[i]=round(xvals[i],2)
fluxinterp_s = np.interp(xvals, gamma, fluxrat_s)

'''
# define some gammas and list the conversion factors from PIMMS, plus the ratio to the full band flux
gamma=np.arange(1.3,2.4,0.1)
#cf_f=[6.243,6.355,6.456,6.544,6.617,6.674,6.712,6.731,6.731,6.709,6.668]
#cf_s=[9.592,9.391,9.186,8.976,8.763,8.547,8.328,8.107,7.885,7.662,7.438]
#cf_h=[4.798,4.864,4.930,4.996,5.061,5.126,5.191,5.255,5.318,5.381,5.442]
#cf_f=np.array(cf_f)*1e10
#cf_s=np.array(cf_s)*1e10
#cf_h=np.array(cf_h)*1e10

fluxrat_s=[0.3016,0.3295,0.3587,0.3890,0.4203,0.4524,0.4849,0.5176,0.5502,0.5825,0.6141]

# now interpolate these binned functions to get more finely sampled curves
xvals = np.linspace(1.3, 2.3, 101)
#yinterp_f = np.interp(xvals, gamma, cf_f)
#yinterp_s = np.interp(xvals, gamma, cf_s)
#yinterp_h = np.interp(xvals, gamma, cf_h)
for i in range(len(xvals)):
	xvals[i]=round(xvals[i],2)
fluxinterp_s = np.interp(xvals, gamma, fluxrat_s)
'''
#############################


#take list of sources in input to simulation
(flux_cdwfs,ra_cdwfs,dec_cdwfs)=np.genfromtxt(wd+'poiss_rand_'+band+'.dat',unpack=True,skip_header=1)
### NEED TO FILTER THESE SOURCES WITH THE TOTAL FOV OF THE CDWFS, SOME OF THEM ARE OUTSIDE 
### AND CANNOT BE MATCHED BY DEFINITION
w=open(wd+'poiss_rand_'+band+'_filtered.dat','w')
w.write('Flux \t RA \t DEC \t Gamma \n')
my_obs = FOVFiles('@'+wd+'fov.lis')
for i in range(len(ra_cdwfs)):
	myobs = my_obs.inside(ra_cdwfs[i], dec_cdwfs[i])
	if len(myobs) > 0:
		w.write(str(flux_cdwfs[i])+' \t '+str(ra_cdwfs[i])+' \t '+str(dec_cdwfs[i])+' \n')
w.close()
#print(len(ra_cdwfs))

#take filtered list of sources in input to simulation
(flux_cdwfs,ra_cdwfs,dec_cdwfs)=np.genfromtxt(wd+'poiss_rand_'+band+'_filtered.dat',unpack=True,skip_header=1)
'''
(flux_cdwfs,ra_cdwfs,dec_cdwfs,gamma_cdwfs)=np.genfromtxt(wd+'poiss_rand_'+band+'_filtered.dat',unpack=True,skip_header=1)
for m in range(len(gamma_cdwfs)):
	if band=='soft':
		flux_ratio=fluxinterp_s[xvals==gamma_cdwfs[m]]
		flux_cdwfs[m]=flux_ratio*flux_cdwfs[m]
	elif band=='hard':
		flux_ratio=1-fluxinterp_s[xvals==gamma_cdwfs[m]]
		flux_cdwfs[m]=flux_ratio*flux_cdwfs[m]
'''

# Sort them to start from the bright ones
ra_cdwfs=ra_cdwfs[::-1]
dec_cdwfs=dec_cdwfs[::-1]
flux_cdwfs=flux_cdwfs[::-1]

#take catalog of detected sources (wavdetect, full mosaic 4x4, 5e-5, only bkgmap) cleaned by multiple sources (cat1)
cat1=fits.open(wd+simfolder+str(sys.argv[2])+'cdwfs_'+band+'_sim_cat1.fits')

#Loop on the cuts
#cut detected sample at a given probability threshold 
cut=np.logspace(np.log10(1e-6),np.log10(1e-1),15)

detected,spurious=[],[]
for l in range(len(cut)):
	ra_k=cat1[1].data['RA']
	dec_k=cat1[1].data['DEC']
	cts_full=cat1[1].data['TOT']
	prob=cat1[1].data['PROB']
	r90=cat1[1].data['AV_R90']
	flux_k=cat1[1].data['FLUX']
	#detml=cat1[1].data['DETML']
	#prob=np.e**(-detml)
	
	ra_k=ra_k[prob<=cut[l]]
	dec_k=dec_k[prob<=cut[l]]
	cts_full=cts_full[prob<=cut[l]]
	r90=r90[prob<=cut[l]]
	#detml=detml[prob<=cut[l]]
	flux_k=flux_k[prob<=cut[l]]
	prob=prob[prob<=cut[l]]
	
	detected.append(len(ra_k))
	pool=build_struct(ra_k,dec_k,r90,flux_k,prob)
	
	#take list of sources in input to simulation
	#(flux_cdwfs,ra_cdwfs,dec_cdwfs,gamma_cdwfs)=np.genfromtxt(wd+'poiss_rand_lehmerx20.dat',unpack=True,skip_header=1)
	### NEED TO FILTER THESE SOURCES WITH THE TOTAL FOV OF THE CDWFS, SOME OF THEM ARE OUTSIDE 
	### AND CANNOT BE MATCHED BY DEFINITION
	#w=open(wd+'poiss_rand_lehmerx20_filtered.dat','w')
	#w.write('Full flux \t RA \t DEC \t Gamma \n')
	#my_obs = FOVFiles('@'+wd+'fov.lis')
	#for i in range(len(ra_cdwfs)):
	#	myobs = my_obs.inside(ra_cdwfs[i], dec_cdwfs[i])
	#	if len(myobs) > 0:
	#		w.write(str(flux_cdwfs[i])+' \t '+str(ra_cdwfs[i])+' \t '+str(dec_cdwfs[i])+' \t '+str(gamma_cdwfs[i])+'\n')
	#w.close()
	#print(len(ra_cdwfs))
	print('Starting to match...')
	t_in=time.time()
	flux_inp,flux_out=[],[]
	unmatched,blendings=0,0
	newpool=pool
	for i in range(len(ra_cdwfs)):
		input_source=[ra_cdwfs[i],dec_cdwfs[i]]
	
		found=0
		counterparts=[]
		count=0
	
		#print(len(newpool[:,0]))

		delta = 0.014 #(0.028 ~100")
		ra_d_filt=newpool[:,0][(newpool[:,0]>=ra_cdwfs[i]-delta) & (newpool[:,0]<=ra_cdwfs[i]+delta) & (newpool[:,1]>=dec_cdwfs[i]-delta) & (newpool[:,1]<=dec_cdwfs[i]+delta)]
		dec_d_filt=newpool[:,1][(newpool[:,0]>=ra_cdwfs[i]-delta) & (newpool[:,0]<=ra_cdwfs[i]+delta) & (newpool[:,1]>=dec_cdwfs[i]-delta) & (newpool[:,1]<=dec_cdwfs[i]+delta)]
		flux_d_filt=newpool[:,3][(newpool[:,0]>=ra_cdwfs[i]-delta) & (newpool[:,0]<=ra_cdwfs[i]+delta) & (newpool[:,1]>=dec_cdwfs[i]-delta) & (newpool[:,1]<=dec_cdwfs[i]+delta)]
		r90_d_filt=newpool[:,2][(newpool[:,0]>=ra_cdwfs[i]-delta) & (newpool[:,0]<=ra_cdwfs[i]+delta) & (newpool[:,1]>=dec_cdwfs[i]-delta) & (newpool[:,1]<=dec_cdwfs[i]+delta)]
		prob_d_filt=newpool[:,4][(newpool[:,0]>=ra_cdwfs[i]-delta) & (newpool[:,0]<=ra_cdwfs[i]+delta) & (newpool[:,1]>=dec_cdwfs[i]-delta) & (newpool[:,1]<=dec_cdwfs[i]+delta)]
		counterparts,probabilities,distances,mat_rad=[],[],[],[]
		for j in range(len(ra_d_filt)):
			cdwfs_source=[ra_d_filt[j],dec_d_filt[j]]
			match_rad=1.1*r90_d_filt[j]
			d=distance(input_source,cdwfs_source)
			if d <= match_rad: #found a match
				if found==0: #it's the first
					found=1
					#flux_inp.append(flux_cdwfs_filt[j])
					#flux_out.append(flux_d[i])
					counterparts.append(flux_d_filt[j])
					probabilities.append(prob_d_filt[j])
					distances.append(d)
					mat_rad.append(match_rad)
				else: #it's not the first match	
					count=count+2
					blendings=blendings+1
					counterparts.append(flux_d_filt[j])
					probabilities.append(prob_d_filt[j])
					distances.append(d)
					mat_rad.append(match_rad)
					#print('Found the '+str(count)+'nd/rd/th counterpart to '+str(kenter_source))
		if found == 1:
			if len(counterparts) == 1:
				flux_inp.append(flux_cdwfs[i])
				flux_out.append(counterparts[0])
				#r90eff.append(mat_rad[0])
				#deff.append(distances[0])
				#print(len(newpool[:,0]))
				#print(len(newpool[newpool[:,3]==counterparts[0]]),newpool[newpool[:,3]==counterparts[0]])
				newpool=np.delete(newpool,np.where(newpool[:,3]==counterparts[0]),0)
				#print(len(newpool[:,0]))
			else:
				counterparts=np.array(counterparts)
				probabilities=np.array(probabilities)
				mat_rad=np.array(mat_rad)
				distances=np.array(distances)
				flux_inp.append(flux_cdwfs[i])
				flux_out.append(counterparts[probabilities==np.min(probabilities)])
				#r90eff.append(mat_rad[probabilities==np.min(probabilities)])
				#deff.append(distances[probabilities==np.min(probabilities)])
				#print(len(newpool[newpool[:,3]==counterparts[probabilities==np.min(probabilities)]]),newpool[newpool[:,3]==counterparts[probabilities==np.min(probabilities)]])
				#print(len(newpool[:,0]))
				newpool=np.delete(newpool,np.where(newpool[:,3]==counterparts[probabilities==np.min(probabilities)]),0)
				#print(len(newpool[:,0]))
				#sys.exit()
	
		else:
			unmatched=unmatched+1

	t_out=time.time()
	print('*'*20)
	print('P cut:',cut[l])
	print(len(ra_cdwfs), 'in input')
	print(len(flux_inp), 'matched')
	print(unmatched, 'unmatched from input')
	print(len(newpool), 'spurious')
	print(blendings,' blendings')
	print(float(t_out-t_in),' seconds for the match.')
	print('*'*20)
	spurious.append(len(newpool))
	'''
	print('Starting to match...')
	t_in=time.time()
	flux_inp,flux_out=[],[]
	unmatched=0
	blendings=0
	cts=[]
	for i in range(len(ra_k)):
		kenter_source=[ra_k[i],dec_k[i]]
		match_rad=1.1*r90[i]
		found=0
		counterparts=[]
		at_least_one=False
		count=0
		for j in range(len(ra_cdwfs)):
			cdwfs_source=[ra_cdwfs[j],dec_cdwfs[j]]
			d=distance(kenter_source,cdwfs_source)
			if d <= match_rad: #found a match
				if found==0: #it's the first
					if flux_k[i]/flux_cdwfs[j] < 10.: #and the flux ratio is less than a factor of 10
						found=1
						flux_inp.append(flux_cdwfs[j])
						flux_out.append(flux_k[i])
						at_least_one=True
				else: #it's not the first match	
					if flux_k[i]/flux_cdwfs[j] < 10.: #and the flux ratio is less than a factor of 10
						count=count+2
						blendings=blendings+1
						#print('Found the '+str(count)+'nd/rd/th counterpart to '+str(kenter_source))
		if found==0:
			unmatched=unmatched+1
			cts.append(cts_full[i])
			#print(ra_k[i],dec_k[i],prob[i])
			#if cts_full[i] > 30.0:
				#print(ra_k[i],dec_k[i])

	#w.close()
	t_out=time.time()
	spurious.append(unmatched)
	print('*'*15)
	print(unmatched, 'unmatched')
	print(len(ra_k), 'detected')
	print('*'*15)
	#print(blendings,' blendings')
	#print(float(t_out-t_in)/60.,' minutes for the match.')
	'''
sp_frac=[]
for k in range(len(spurious)):
	sp_frac.append(float(spurious[k])/float(detected[k])*100.)

# Write out sp fraction as a function of probability cut
w=open(wd+simfolder+str(sys.argv[2])+'cdwfs_'+band+'_sp-frac.dat','w')
for i in range(len(cut)):
	w.write(str(cut[i])+' \t '+str(sp_frac[i])+'\n')
w.close()

#plt.figure()
#plt.plot(cut,sp_frac,'go',linestyle='-')
#plt.axhline(y=1.0)
#plt.axhline(y=3.0)
#plt.xscale('log')
#plt.xlabel('Prob cut')
#plt.ylabel('Sp fraction (%)')
#plt.tight_layout()
#plt.show()

#sp_frac in the hard band
#sp_frac=np.array([41./3999.,58./4145.,68./4288.,83./4437.,98./4571.,113./4677.,131./4780.,146./4878.,169./4974.,191./5039.])*100.

####################
sys.exit()         #
####################

### Once the cuts are established, make checks and plots. Use 99% reliability
if band=='soft':
	cut=1.4e-4
else:
	cut=6e-5
	
ra_k=cat1[1].data['RA']
dec_k=cat1[1].data['DEC']
cts_full=cat1[1].data['TOT']
prob=cat1[1].data['PROB']
r90=cat1[1].data['AV_R90']
flux_k=cat1[1].data['FLUX']
	
ra_k=ra_k[prob<=cut]
dec_k=dec_k[prob<=cut]
cts_full=cts_full[prob<=cut]
r90=r90[prob<=cut]
flux_k=flux_k[prob<=cut]
prob=prob[prob<=cut]

detected=len(ra_k)

print('Starting to match...')
t_in=time.time()
flux_inp,flux_out=[],[]
r90eff,deff=[],[]
unmatched=0
blendings=0
cts=[]
for i in range(len(ra_k)):
	kenter_source=[ra_k[i],dec_k[i]]
	match_rad=1.1*r90[i]
	found=0
	counterparts=[]
	at_least_one=False
	count=0
	for j in range(len(ra_cdwfs)):
		cdwfs_source=[ra_cdwfs[j],dec_cdwfs[j]]
		d=distance(kenter_source,cdwfs_source)
		if d <= match_rad: #found a match
			if found==0: #it's the first
				if flux_k[i]/flux_cdwfs[j] < 10.: #and the flux ratio is less than a factor of 10
					found=1
					r90eff.append(r90[i])
					deff.append(d)
					flux_inp.append(flux_cdwfs[j])
					flux_out.append(flux_k[i])
					at_least_one=True
			else: #it's not the first match	
				if flux_k[i]/flux_cdwfs[j] < 10.: #and the flux ratio is less than a factor of 10
					count=count+2
					blendings=blendings+1
					#print('Found the '+str(count)+'nd/rd/th counterpart to '+str(kenter_source))
	if found==0:
		unmatched=unmatched+1
		cts.append(cts_full[i])
		#print(ra_k[i],dec_k[i],prob[i])
		#if cts_full[i] > 30.0:
			#print(ra_k[i],dec_k[i])

t_out=time.time()

print('*'*15)
print(unmatched, 'unmatched')
print(len(ra_k), 'detected')
print('*'*15)
print(blendings,' blendings')
print(float(t_out-t_in)/60.,' minutes for the match.')

r90eff=np.array(r90eff)
deff=np.array(deff)

print(len(deff),len(deff[deff>2.0]),'Distances and N > 2 arsec')

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
plt.savefig(wd+'cdwfs_dist-r90.pdf',format='pdf')
#plt.show()

n1,o1=np.histogram(deff,bins=15)
bcenters1=list((o1[i+1]+o1[i])/2. for i in range(len(o1)-1))
cum=np.cumsum(n1)/float(np.sum(n1))
f,(ax1,ax2)=plt.subplots(2,1,sharex=True)
ax1.hist(deff,bins=15)
ax1.set_ylabel('N')
ax2.plot(bcenters1,cum,'b-')
ax2.set_xlabel('Distance [arcsec]')
ax2.set_ylabel('Cumulative fraction')
#plt.subplots_adjust(hspace=0.1)
plt.tight_layout()
plt.savefig(wd+'cdwfs_dist-cum.pdf',format='pdf')
#plt.show()

x=np.logspace(np.log10(1e-17),np.log10(1e-11),30)
plt.figure()
plt.plot(x,x,'k--')
plt.plot(flux_inp,flux_out,'g.')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$F_{in}$ (cgs)')
plt.ylabel(r'$F_{out}$ (cgs)')
plt.axis([3e-16,3e-12,3e-16,3e-12])
plt.tight_layout()
plt.show()
#plt.savefig(wd+'cdwfs_'+band+'_fin-fout.pdf',format='pdf')
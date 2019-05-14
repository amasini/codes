#completeness and flux_in VS flux_out

#match simulated sources detected to input ones, computing the reliability and completeness
# of the sample
import numpy as np
import sys
from astropy.io import fits
import matplotlib.pyplot as plt
import time
from ciao_contrib.region.check_fov import FOVFiles

def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*xx*3600)**2 +((pointa[1]-pointb[1])*3600)**2)

def build_struct(a,b,c,d,e):
	s=[]
	for j in range(len(a)):
		s.append([a[j],b[j],c[j],d[j],e[j]])
	s=np.array(s)
	return s
	
wd='/Users/alberto/Desktop/XBOOTES/'

band='broad'

#take catalog of detected sources (wavdetect, full mosaic 4x4, 5e-5, only bkgmap)
cat1=fits.open(wd+'cdwfs_'+band+'_cat1_sim.fits')

ra_d=cat1[1].data['RA']
dec_d=cat1[1].data['DEC']
cts_full=cat1[1].data['TOT']
prob=cat1[1].data['PROB']
r90=cat1[1].data['AV_R90']
flux_d=cat1[1].data['FLUX']
print(len(ra_d),'total sample')

#define cut at 1% of spurious fraction
if band=='broad':
	cut=6e-5
elif band=='soft':
	cut=1.4e-4
elif band=='hard':
	cut=6e-5

#cut detected sample at a given probability threshold 
ra_d=ra_d[prob<=cut]
dec_d=dec_d[prob<=cut]
cts_full=cts_full[prob<=cut]
r90=r90[prob<=cut]
flux_d=flux_d[prob<=cut]
prob=prob[prob<=cut]

pool=build_struct(ra_d,dec_d,r90,flux_d,prob)
print(len(pool),'above rel threshold')

#take list of sources in input to simulation
#(flux_cdwfs,ra_cdwfs,dec_cdwfs)=np.genfromtxt(wd+'poiss_rand_'+band+'_lehmer.dat',unpack=True,skip_header=1)
### NEED TO FILTER THESE SOURCES WITH THE TOTAL FOV OF THE CDWFS, SOME OF THEM ARE OUTSIDE 
### AND CANNOT BE MATCHED BY DEFINITION
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
(flux_cdwfs,ra_cdwfs,dec_cdwfs)=np.genfromtxt(wd+'poiss_rand_lehmerx20_filtered.dat',unpack=True,skip_header=1,usecols=[0,1,2])
if band=='soft':
	flux_cdwfs=0.33*flux_cdwfs
elif band=='hard':
	flux_cdwfs=0.67*flux_cdwfs

# Sort them to start from the bright ones
ra_cdwfs=ra_cdwfs[::-1]
dec_cdwfs=dec_cdwfs[::-1]
flux_cdwfs=flux_cdwfs[::-1]

print('Starting to match...')
t_in=time.time()
flux_inp,flux_out=[],[]
unmatched,blendings=0,0
newpool=pool

#print(newpool[0],len(newpool[:,0]))
#myra=newpool[:,0][0]
#newpool=newpool[:,0][(newpool[:,0]>=myra-0.014) & (newpool[:,0]<=myra+0.014)]
#print(newpool,len(newpool))
#newpool2=np.delete(newpool,np.where(newpool[:,0]==myra),0)
#print(newpool2[0],len(newpool2))
#sys.exit()

r90eff,deff=[],[]
for i in range(len(ra_cdwfs)):
	input_source=[ra_cdwfs[i],dec_cdwfs[i]]
	
	found=0
	counterparts=[]
	count=0
	
	print(len(newpool[:,0]))

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
			r90eff.append(mat_rad[0])
			deff.append(distances[0])
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
			r90eff.append(mat_rad[probabilities==np.min(probabilities)])
			deff.append(distances[probabilities==np.min(probabilities)])
			#print(len(newpool[newpool[:,3]==counterparts[probabilities==np.min(probabilities)]]),newpool[newpool[:,3]==counterparts[probabilities==np.min(probabilities)]])
			#print(len(newpool[:,0]))
			newpool=np.delete(newpool,np.where(newpool[:,3]==counterparts[probabilities==np.min(probabilities)]),0)
			#print(len(newpool[:,0]))
			#sys.exit()
	
	else:
		unmatched=unmatched+1

t_out=time.time()
print(len(ra_cdwfs), 'in input')
print(len(flux_inp), 'matched')
print(unmatched, 'unmatched')
print(blendings,' blendings')
print(float(t_out-t_in),' seconds for the match.')

n,be=np.histogram(deff,bins=25)
be2=list((be[i+1]+be[i])/2. for i in range(len(be)-1))
nsum=np.cumsum(n)/float(np.sum(n))

f,(ax1,ax2)=plt.subplots(2,1,sharex=True)
ax1.hist(deff,bins=25)
ax1.set_ylabel('N')
ax2.plot(be2,nsum,'b-')
ax2.set_ylabel('Cumulative fraction')
ax2.set_xlabel('Distance [arcsec]')
plt.savefig(wd+'cdwfs_dist-cum.pdf',format='pdf')


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
plt.savefig(wd+'cdwfs_dist-r90.pdf',format='pdf')
#plt.show()

sys.exit()
'''
w=open(wd+'cdwfs_'+band+'_cat1_sim_unmatched.reg','w')
w2=open(wd+'cdwfs_'+band+'_cat1_sim_matched.dat','w')
unmatched=0
blendings=0
cts=[]
for i in range(len(ra_d)):
	kenter_source=[ra_d[i],dec_d[i]]
	match_rad=1.1*r90[i]
	found=0
	counterparts=[]
	at_least_one=False
	count=0
	
	delta = 0.014 #(0.028 ~100")
	ra_cdwfs_filt=ra_cdwfs[(ra_cdwfs>=ra_d[i]-delta) & (ra_cdwfs<=ra_d[i]+delta) & (dec_cdwfs>=dec_d[i]-delta) & (dec_cdwfs<=dec_d[i]+delta)]
	dec_cdwfs_filt=dec_cdwfs[(ra_cdwfs>=ra_d[i]-delta) & (ra_cdwfs<=ra_d[i]+delta) & (dec_cdwfs>=dec_d[i]-delta) & (dec_cdwfs<=dec_d[i]+delta)]
	flux_cdwfs_filt=flux_cdwfs[(ra_cdwfs>=ra_d[i]-delta) & (ra_cdwfs<=ra_d[i]+delta) & (dec_cdwfs>=dec_d[i]-delta) & (dec_cdwfs<=dec_d[i]+delta)]

	counterparts=[]
	for j in range(len(ra_cdwfs_filt)):
		cdwfs_source=[ra_cdwfs_filt[j],dec_cdwfs_filt[j]]
		d=distance(kenter_source,cdwfs_source)
		if d <= match_rad: #found a match
			if found==0: #it's the first
				if flux_d[i]/flux_cdwfs_filt[j] < 10.: #and the flux ratio is less than a factor of 10
					found=1
					#flux_inp.append(flux_cdwfs_filt[j])
					#flux_out.append(flux_d[i])
					counterparts.append(flux_cdwfs_filt[j])
					w2.write(str(ra_d[i])+' \t '+str(dec_d[i])+' \t '+str(prob[i])+' \t '+str(flux_cdwfs_filt[j])+' \t '+str(flux_d[i])+' \t '+str(d)+' \n')
					at_least_one=True
			else: #it's not the first match	
				if flux_d[i]/flux_cdwfs_filt[j] < 10.: #and the flux ratio is less than a factor of 10
					count=count+2
					blendings=blendings+1
					counterparts.append(flux_cdwfs_filt[j])
					#print('Found the '+str(count)+'nd/rd/th counterpart to '+str(kenter_source))
	if found == 1:
		if len(counterparts) == 1:
			flux_inp.append(counterparts[0])
			flux_out.append(flux_d[i])
		else:
			flux_inp.append(np.max(counterparts))
			flux_out.append(flux_d[i])
		
	else:
		unmatched=unmatched+1
		w.write('circle('+str(ra_d[i])+'d, '+str(dec_d[i])+'d, 10\") #width=2 color=cyan \n')
		cts.append(cts_full[i])
		#print(ra_k[i],dec_k[i],prob[i])
		#if cts_full[i] > 30.0:
			#print(ra_k[i],dec_k[i])

w.close()
w2.close()
t_out=time.time()
print(len(ra_d), 'detected')
print(len(flux_inp), 'matched')
print(unmatched, 'unmatched')
print(blendings,' blendings')
print(float(t_out-t_in),' seconds for the match.')
'''


plt.figure()
plt.plot(flux_inp,flux_out,'k.')
plt.xscale('log')
plt.yscale('log')
plt.show()

bins=np.logspace(np.log10(5e-17),np.log10(1e-12),50)
a,b=np.histogram(flux_cdwfs,bins=bins)
bincenters=list((b[i+1]+b[i])/2. for i in range(len(b)-1))
#plt.figure()
#plt.plot(bincenters,a,'go')
#plt.xscale('log')
#plt.yscale('log')
#plt.show()
#sys.exit()
a2,b2=np.histogram(flux_inp,bins=bins)

a=a/1.
a2=a2/1.

cum1=np.array(list(reversed(np.cumsum(list(reversed(a))))))
cum2=np.array(list(reversed(np.cumsum(list(reversed(a2))))))

cum=cum2/cum1
#cum[cum>1.]=1.

w=open(wd+'cdwfs_completeness.dat','w')
for i in range(len(bincenters)):
	w.write(str(bincenters[i])+' \t '+str(cum[i])+'\n')
w.close()

plt.figure()
plt.plot(bincenters,cum,'-',color='green')
plt.xscale('log')
plt.xlabel(r'0.5-7 keV flux [erg cm$^{-2}$ s$^(-1)$]')
plt.ylabel(r'Completeness')
plt.grid()
#plt.savefig(wd+'cdwfs_completeness.pdf',format='pdf')
plt.show()

#sys.exit()

#rat=[]
#for i in range(len(a)):
#	if float(a[i]) != 0.:
#		rat.append(float(a2[i])/float(a[i]))
#	else:
#		rat.append(1)
#rat=np.array(rat)*9.306
#cum=np.cumsum(rat)/np.sum(rat)
#cum=np.array(list(reversed(np.cumsum(list(reversed(rat))))))

#(f,k)=np.genfromtxt(wd+'kenter05_sens_05-7keV.dat',unpack=True)
#k=k/np.max(k)



rebin_factor=4.
#ecf=1.825E-11 # gamma=1.8 
ecf=1.847E-11 # gamma=1.4
scale=(0.492/3600.)*rebin_factor #pixel size in deg
arcsec2pix=scale*3600

sensb=fits.open(wd+'cdwfs_broad_sens_1.4Cy18.fits')
sens2b=sensb[0].data
sens2b=sens2b[sens2b>0]
print('Here.')

binsb=np.logspace(np.log10(min(sens2b)),np.log10(1e-13),100)
ab,bb=np.histogram(sens2b,bins=binsb)
areab=ab*scale**2
centersb=list((binsb[i+1]+binsb[i])/2. for i in range(len(binsb)-1))
cumb=np.cumsum(areab)

sensc=fits.open(wd+'cdwfs_broad_sens_r90.fits')
sens2c=sensc[0].data
sens2c=sens2c[sens2c>0]
print('Here.')

binsc=np.logspace(np.log10(min(sens2c)),np.log10(1e-13),100)
ac,bc=np.histogram(sens2c,bins=binsc)
areac=ac*scale**2
centersc=list((binsc[i+1]+binsc[i])/2. for i in range(len(binsc)-1))
cumc=np.cumsum(areac)

cum2=cum*np.max(cumb)

w=open(wd+'skycov.dat','w')
for k in range(len(bincenters)):
	w.write(str(bincenters[k])+' \t '+str(cum2[k])+'\n')
w.close()

plt.figure()
plt.plot(bincenters,cum2,'-',color='green',label='Completeness rescaled')
plt.plot(centersb,cumb,'b-',label='Sensitivity, R50')
plt.plot(centersc,cumc,'b--',label='Sensitivity, R90')
#plt.plot(f,k,'k--')
plt.xscale('log')
plt.xlabel(r'0.5-7 keV flux [erg cm$^{-2}$ s$^(-1)$]')
plt.ylabel(r'Area [deg$^2$]')
plt.legend()
plt.tight_layout()
#plt.yscale('log')
plt.show()
#plt.savefig(wd+'cdwfs_sensitivity.pdf',format='pdf')
sys.exit()

'''
------------
(flux_sim,ra_sim,dec_sim)=np.genfromtxt(wd+'poiss_rand_lehmerx20_filtered.dat',unpack=True,skip_header=1)
plt.figure()
for band in ['broad','soft','hard']:
	if band=='soft':
		flux_sim=0.33*flux_sim
	elif band=='hard':
		flux_sim=0.67*flux_sim

	(flux_inp,flux_out)=np.genfromtxt(wd+'sim_all/cdwfs_'+band+'_sim_poiss_matched.dat',unpack=True,usecols=[3,4])

	bins=np.logspace(np.log10(5e-17),np.log10(9e-13),20)
	a,b=np.histogram(flux_sim,bins=bins)
	bincenters=list((b[i+1]+b[i])/2. for i in range(len(b)-1))

	a2,b2=np.histogram(flux_inp,bins=bins)
	rat=[]
	#print(a)
	#print(a2)
	for i in range(len(a)):
		if float(a[i]) != 0.:
			rat.append(float(a2[i])/float(a[i]))
		else:
			rat.append(1)
	print(band,rat)
	
	plt.plot(bincenters,rat,'k-',label=band)
	if band=='broad':
		plt.xlabel('0.5-7 keV Flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=13)
	elif band=='soft':
		plt.xlabel('0.5-2 keV Flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=13)
	else:
		plt.xlabel('2-7 keV Flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=13)

plt.xscale('log')
plt.ylabel('Completeness',fontsize=13)
#plt.savefig(wd+'cdwfs_'+band+'_completeness.pdf',format='pdf',dpi=1000)
plt.show()
sys.exit()

band='hard'

(flux_sim,ra_sim,dec_sim)=np.genfromtxt(wd+'poiss_rand_lehmerx20_filtered.dat',unpack=True,skip_header=1)
if band=='soft':
	flux_sim=0.33*flux_sim
elif band=='hard':
	flux_sim=0.67*flux_sim

(ra,dec,flux_inp,flux_out)=np.genfromtxt(wd+'sim_all/cdwfs_'+band+'_sim_poiss_matched.dat',unpack=True,usecols=[0,1,3,4])

bins=np.logspace(np.log10(5e-17),np.log10(9e-13),20)
a,b=np.histogram(flux_sim,bins=bins)
bincenters=list((b[i+1]+b[i])/2. for i in range(len(b)-1))

a2,b2=np.histogram(flux_inp,bins=bins)

print(a)
print(a2)
a2[-4]=40
a2[-2]=10
for j in range(len(flux_inp)):
	if (flux_inp[j] <= 6e-13 and flux_inp[j] >=4e-13):
		print(ra[j],dec[j],flux_inp[j],flux_out[j],flux_inp[j]/flux_out[j])

rat=list(float(a2[i])/float(a[i]) for i in range(len(a)))
#print(rat)

plt.figure()
plt.plot(bincenters,rat,'k-')
plt.xscale('log')
if band=='broad':
	plt.xlabel('0.5-7 keV Flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=13)
elif band=='soft':
	plt.xlabel('0.5-2 keV Flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=13)
else:
	plt.xlabel('2-7 keV Flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=13)
plt.ylabel('Completeness',fontsize=13)
#plt.savefig(wd+'cdwfs_'+band+'_completeness.pdf',format='pdf',dpi=1000)
plt.show()
'''

color='black'
plt.figure()
for band in ['broad','soft','hard']:
	(flux_sim,ra_sim,dec_sim)=np.genfromtxt(wd+'poiss_rand_lehmerx20_filtered.dat',unpack=True,skip_header=1)
	if band=='soft':
		color='blue'
		flux_sim=0.33*flux_sim
	elif band=='hard':
		color='red'
		flux_sim=0.67*flux_sim

	(flux_inp,flux_out)=np.genfromtxt(wd+'sim_all/cdwfs_'+band+'_sim_poiss_matched.dat',unpack=True,usecols=[3,4])

	bins=np.logspace(np.log10(5e-17),np.log10(9e-13),20)
	a,b=np.histogram(flux_sim,bins=bins)
	bincenters=list((b[i+1]+b[i])/2. for i in range(len(b)-1))

	a2,b2=np.histogram(flux_inp,bins=bins)
	if band=='hard':
		a2[-4]=40
		a2[-2]=10
	rat=[]
	#print(a)
	#print(a2)
	for i in range(len(a)):
		if float(a[i]) != 0.:
			rat.append(float(a2[i])/float(a[i]))
		else:
			rat.append(1)
	#print(band,rat)
	rat=np.array(rat)*9.306
	#cum=np.cumsum(rat)/np.sum(rat)*9.306
	plt.plot(bincenters,rat,'-',color=color,label=band)

(f,k)=np.genfromtxt(wd+'kenter05_sens_05-7keV.dat',unpack=True)

plt.plot(f,k,'k--')
plt.xlabel('Flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=13)
plt.xscale('log')
plt.ylabel('Area [deg2]',fontsize=13)
plt.legend()
plt.tight_layout()
plt.savefig(wd+'cdwfs_all_completeness.pdf',format='pdf',dpi=1000)
#plt.show()
sys.exit()

x=np.linspace(1e-16,1e-12,100)
plt.figure()
plt.plot(flux_inp,flux_out,'r.')
plt.plot(x,x,'k--')
plt.xscale('log')
plt.yscale('log')
if band=='broad':
	plt.xlabel('0.5-7 keV Input Flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=13)
	plt.ylabel('0.5-7 keV Output Flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=13)
elif band=='soft':
	plt.xlabel('0.5-2 keV Input Flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=13)
	plt.ylabel('0.5-2 keV Output Flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=13)
else:
	plt.xlabel('2-7 keV Input Flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=13)
	plt.ylabel('2-7 keV Output Flux [erg cm$^{-2}$ s$^{-1}$]',fontsize=13)
#plt.hist(cts,bins=100)
plt.show()
#plt.savefig(wd+'cdwfs_'+band+'_Fin-Fout.pdf',format='pdf',dpi=1000)
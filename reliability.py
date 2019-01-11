# reliability
import numpy as np
import sys
from astropy.io import fits
import matplotlib.pyplot as plt
import time
from ciao_contrib.region.check_fov import FOVFiles

def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*xx*3600)**2 +((pointa[1]-pointb[1])*3600)**2)

wd='/Users/alberto/Desktop/XBOOTES/'

band='broad'

#take catalog of detected sources (wavdetect, full mosaic 4x4, 5e-5, only bkgmap)
cat1=fits.open(wd+'sim_all/mosaic_'+band+'_sim_poiss_cat0_3.fits')

ra_k=cat1[1].data['RA']
dec_k=cat1[1].data['DEC']
cts_full=cat1[1].data['TOT']
prob=cat1[1].data['PROB']
r90=cat1[1].data['AV_R90']
flux_k=cat1[1].data['FLUX']
#detml=cat1[1].data['DETML']
#prob=np.e**(-detml)
print(len(ra_k))

#Loop on the cuts
#cut detected sample at a given probability threshold 
cut=np.logspace(np.log10(1e-2),np.log10(3e-2),3)

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

	'''
	#write out the cut sample
	w=open(wd+'cdwfs_'+band+'_sim_poiss_'+str(cut)+'.reg','w')
	for i in range(len(ra_k)):
		w.write('circle('+str(ra_k[i])+'d, '+str(dec_k[i])+'d, '+str(r90[i])+'\") #width=2 color=red \n')
	w.close()
	'''
	'''
	#take list of sources in input to simulation
	(flux_cdwfs,ra_cdwfs,dec_cdwfs)=np.genfromtxt(wd+'poiss_rand_'+band+'_lehmer.dat',unpack=True,skip_header=1)
	### NEED TO FILTER THESE SOURCES WITH THE TOTAL FOV OF THE CDWFS, SOME OF THEM ARE OUTSIDE 
	### AND CANNOT BE MATCHED BY DEFINITION
	w=open(wd+'poiss_rand_'+band+'_lehmer_filtered.dat','w')
	w.write('Flux \t RA \t DEC \n')
	my_obs = FOVFiles('@'+wd+'fov.lis')
	for i in range(len(ra_cdwfs)):
		myobs = my_obs.inside(ra_cdwfs[i], dec_cdwfs[i])
		if len(myobs) > 0:
			w.write(str(flux_cdwfs[i])+' \t '+str(ra_cdwfs[i])+' \t '+str(dec_cdwfs[i])+' \n')
	w.close()
	print(len(ra_cdwfs))
	'''
	#take filtered list of sources in input to simulation
	(flux_cdwfs,ra_cdwfs,dec_cdwfs)=np.genfromtxt(wd+'poiss_rand_lehmerx20_filtered.dat',unpack=True,skip_header=1)
	if band=='soft':
		flux_cdwfs=0.33*flux_cdwfs
	elif band=='hard':
		flux_cdwfs=0.67*flux_cdwfs
	'''
	if band != 'broad':
		(flux_cdwfs,ra_cdwfs,dec_cdwfs)=np.genfromtxt(wd+'poiss_rand_'+band+'_lehmer_filtered.dat',unpack=True,skip_header=1)
	else:
		(flux_cdwfs,ra_cdwfs,dec_cdwfs)=np.genfromtxt(wd+'poiss_rand_lehmerx20_filtered.dat',unpack=True,skip_header=1)
	'''
	
	#print('Starting to match...')
	t_in=time.time()
	#match_rad=5.0
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

sp_frac=[]
for k in range(len(spurious)):
	sp_frac.append(float(spurious[k])/float(detected[k])*100.)

'''
# Write out sp fraction as a funciton of probability cut
w=open(wd+'sim_all/sp_frac_'+band+'.dat','w')
for i in range(len(cut)):
	w.write(str(cut[i])+' \t '+str(sp_frac[i])+'\n')
w.close()
'''

plt.figure()
plt.plot(cut,sp_frac,'k-')
plt.plot(cut,sp_frac,'go')
plt.axhline(y=1.0)
plt.axhline(y=3.0)
plt.xscale('log')
plt.xlabel('Prob cut')
plt.ylabel('Sp fraction (%)')
plt.tight_layout()
plt.show()
#sp_frac in the hard band
#sp_frac=np.array([41./3999.,58./4145.,68./4288.,83./4437.,98./4571.,113./4677.,131./4780.,146./4878.,169./4974.,191./5039.])*100.
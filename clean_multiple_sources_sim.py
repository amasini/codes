#match simulated sources detected to input ones, computing the reliability and completeness
# of the sample - WAY FASTER THAN clean_multiple_sources.py... NOTE THAT COULD BE A FACTOR OF 4 FASTER
# WITH SMARTER CODING
import numpy as np
import sys
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import time
#from ciao_contrib.region.check_fov import FOVFiles

def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*xx*3600)**2 +((pointa[1]-pointb[1])*3600)**2)

def build_struct(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p):
	s=[]
	for jj in range(len(a)):
		s.append([a[jj],b[jj],c[jj],d[jj],e[jj],f[jj],g[jj],h[jj],i[jj],j[jj],k[jj],l[jj],m[jj],n[jj],o[jj],p[jj]])
	s=np.array(s)
	return s

def choose_worst(pointa, pointb):
	"""
	Function that chooses the best between two points
	This version is based on the lower prob of being spurious
	"""
	if pointa[2] >= pointb[2]:
		return pointa, pointb
	else:
		return pointb, pointa 

wd='/Users/alberto/Desktop/XBOOTES/'

simfolder = 'sim_indep/'

band=str(sys.argv[1])

#take catalog of detected sources (wavdetect, full mosaic 4x4, 5e-5, only bkgmap)
cat1=fits.open(wd+simfolder+str(sys.argv[2])+'cdwfs_'+band+'_sim_cat0_exp-psf.fits')

ra_d=cat1[1].data['RA']
dec_d=cat1[1].data['DEC']
cts_full=cat1[1].data['TOT']
prob=cat1[1].data['PROB']
r90=cat1[1].data['AV_R90']
tot=cat1[1].data['TOT']
bkg=cat1[1].data['BKG']
net=cat1[1].data['NET']
enetp=cat1[1].data['E_NET_+']
enetn=cat1[1].data['E_NET_-']
exp=cat1[1].data['EXP']
cr=cat1[1].data['CR']
ecrp=cat1[1].data['E_CR_+']
ecrn=cat1[1].data['E_CR_-']
flux_d=cat1[1].data['FLUX']
efluxp=cat1[1].data['E_FLUX_+']
efluxn=cat1[1].data['E_FLUX_-']


pool=build_struct(ra_d,dec_d,prob,r90,tot,bkg,net,enetp,enetn,exp,cr,ecrp,ecrn,flux_d,efluxp,efluxn)
print(len(pool),'total sample')
'''
#take list of sources in input to simulation
(flux_cdwfs,ra_cdwfs,dec_cdwfs,gamma)=np.genfromtxt(wd+'poiss_rand_lehmer_newgamma.dat',unpack=True,skip_header=1)
### NEED TO FILTER THESE SOURCES WITH THE TOTAL FOV OF THE CDWFS, SOME OF THEM ARE OUTSIDE 
### AND CANNOT BE MATCHED BY DEFINITION
w=open(wd+'poiss_rand_lehmer_newgamma_filtered.dat','w')
w.write('Flux \t RA \t DEC \t Gamma\n')
my_obs = FOVFiles('@'+wd+'fov.lis')
for i in range(len(ra_cdwfs)):
	myobs = my_obs.inside(ra_cdwfs[i], dec_cdwfs[i])
	if len(myobs) > 0:
		w.write(str(flux_cdwfs[i])+' \t '+str(ra_cdwfs[i])+' \t '+str(dec_cdwfs[i])+' \t '+str(gamma[i])+'\n')
w.close()
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
'''
print('Starting to match...')
t_in=time.time()

newpool=pool

for i in range(len(ra_d)):
	input_source=[ra_d[i],dec_d[i]]
	
	multiple=0
	counterparts=[]
	
	#print(len(newpool[:,0]))

	delta = 0.007 #(0.028 ~100")
	ra_d_filt=newpool[:,0][(newpool[:,0]>=ra_d[i]-delta) & (newpool[:,0]<=ra_d[i]+delta) & (newpool[:,1]>=dec_d[i]-delta) & (newpool[:,1]<=dec_d[i]+delta)]
	dec_d_filt=newpool[:,1][(newpool[:,0]>=ra_d[i]-delta) & (newpool[:,0]<=ra_d[i]+delta) & (newpool[:,1]>=dec_d[i]-delta) & (newpool[:,1]<=dec_d[i]+delta)]
	prob_d_filt=newpool[:,2][(newpool[:,0]>=ra_d[i]-delta) & (newpool[:,0]<=ra_d[i]+delta) & (newpool[:,1]>=dec_d[i]-delta) & (newpool[:,1]<=dec_d[i]+delta)]
	r90_d_filt=newpool[:,3][(newpool[:,0]>=ra_d[i]-delta) & (newpool[:,0]<=ra_d[i]+delta) & (newpool[:,1]>=dec_d[i]-delta) & (newpool[:,1]<=dec_d[i]+delta)]
	tot_d_filt=newpool[:,4][(newpool[:,0]>=ra_d[i]-delta) & (newpool[:,0]<=ra_d[i]+delta) & (newpool[:,1]>=dec_d[i]-delta) & (newpool[:,1]<=dec_d[i]+delta)]
	bkg_d_filt=newpool[:,5][(newpool[:,0]>=ra_d[i]-delta) & (newpool[:,0]<=ra_d[i]+delta) & (newpool[:,1]>=dec_d[i]-delta) & (newpool[:,1]<=dec_d[i]+delta)]
	net_d_filt=newpool[:,6][(newpool[:,0]>=ra_d[i]-delta) & (newpool[:,0]<=ra_d[i]+delta) & (newpool[:,1]>=dec_d[i]-delta) & (newpool[:,1]<=dec_d[i]+delta)]
	enetp_d_filt=newpool[:,7][(newpool[:,0]>=ra_d[i]-delta) & (newpool[:,0]<=ra_d[i]+delta) & (newpool[:,1]>=dec_d[i]-delta) & (newpool[:,1]<=dec_d[i]+delta)]
	enetn_d_filt=newpool[:,8][(newpool[:,0]>=ra_d[i]-delta) & (newpool[:,0]<=ra_d[i]+delta) & (newpool[:,1]>=dec_d[i]-delta) & (newpool[:,1]<=dec_d[i]+delta)]
	exp_d_filt=newpool[:,9][(newpool[:,0]>=ra_d[i]-delta) & (newpool[:,0]<=ra_d[i]+delta) & (newpool[:,1]>=dec_d[i]-delta) & (newpool[:,1]<=dec_d[i]+delta)]
	cr_d_filt=newpool[:,10][(newpool[:,0]>=ra_d[i]-delta) & (newpool[:,0]<=ra_d[i]+delta) & (newpool[:,1]>=dec_d[i]-delta) & (newpool[:,1]<=dec_d[i]+delta)]
	ecrp_d_filt=newpool[:,11][(newpool[:,0]>=ra_d[i]-delta) & (newpool[:,0]<=ra_d[i]+delta) & (newpool[:,1]>=dec_d[i]-delta) & (newpool[:,1]<=dec_d[i]+delta)]
	ecrn_d_filt=newpool[:,12][(newpool[:,0]>=ra_d[i]-delta) & (newpool[:,0]<=ra_d[i]+delta) & (newpool[:,1]>=dec_d[i]-delta) & (newpool[:,1]<=dec_d[i]+delta)]
	flux_d_filt=newpool[:,13][(newpool[:,0]>=ra_d[i]-delta) & (newpool[:,0]<=ra_d[i]+delta) & (newpool[:,1]>=dec_d[i]-delta) & (newpool[:,1]<=dec_d[i]+delta)]
	efluxp_d_filt=newpool[:,14][(newpool[:,0]>=ra_d[i]-delta) & (newpool[:,0]<=ra_d[i]+delta) & (newpool[:,1]>=dec_d[i]-delta) & (newpool[:,1]<=dec_d[i]+delta)]
	efluxn_d_filt=newpool[:,15][(newpool[:,0]>=ra_d[i]-delta) & (newpool[:,0]<=ra_d[i]+delta) & (newpool[:,1]>=dec_d[i]-delta) & (newpool[:,1]<=dec_d[i]+delta)]
	
	
	counterparts,probabilities,distances,mat_rad=[],[],[],[]
	for j in range(len(ra_d_filt)):
		cdwfs_source=[ra_d_filt[j],dec_d_filt[j]]
		match_rad=1.1*r90_d_filt[j]
		d=distance(input_source,cdwfs_source)
		if (d <= match_rad) and (d != 0.0): #found a match
			
			multiple=1
			
			counterparts.append(np.array([ra_d_filt[j],dec_d_filt[j],prob_d_filt[j],r90_d_filt[j],tot_d_filt[j],bkg_d_filt[j],net_d_filt[j],enetp_d_filt[j],enetn_d_filt[j],exp_d_filt[j],cr_d_filt[j],ecrp_d_filt[j],ecrn_d_filt[j],flux_d_filt[j],efluxp_d_filt[j],efluxn_d_filt[j]]))
			print(pool[i])
			print('+'*10)
			print(counterparts)
			
	if multiple == 1:
	
		if len(counterparts) == 1:
			to_delete=choose_worst(pool[i],counterparts[-1])
			#print(len(newpool[:,0]))
			#print(len(newpool[newpool[:,3]==counterparts[0]]),newpool[newpool[:,3]==counterparts[0]])
			#newpool=np.delete(newpool,np.where(newpool[:,3]==to_delete[0]),0)
			newpool=np.delete(newpool,np.where(newpool[:,0]==to_delete[0][0]),0)

		else:

			if len(counterparts) == 2:
				to_delete_count=choose_worst(counterparts[0],counterparts[1])
				newpool=np.delete(newpool,np.where(newpool[:,0]==to_delete_count[0][0]),0)
				to_delete=choose_worst(pool[i],to_delete_count[1])
				newpool=np.delete(newpool,np.where(newpool[:,0]==to_delete[0][0]),0)
			else:
				print('im here')
				sys.exit()


t_out=time.time()

'''
sys.exit()

ra,dec,prob,av_r90,tot,bkg,net,e_net_up,e_net_lo,exptime,cr,e_cr_up,e_cr_lo,flux,e_flux_up,e_flux_lo=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
for i in range(len(final_list)):
	ra.append(final_list[i][0])
	dec.append(final_list[i][1])
	
	prob.append(final_list[i][2])
	av_r90.append(final_list[i][3])
	tot.append(final_list[i][4])
	bkg.append(final_list[i][5])
	net.append(final_list[i][6])
	e_net_up.append(final_list[i][7])
	e_net_lo.append(final_list[i][8])
	exptime.append(final_list[i][9])
	cr.append(final_list[i][10])
	e_cr_up.append(final_list[i][11])
	e_cr_lo.append(final_list[i][12])
	flux.append(final_list[i][13])
	e_flux_up.append(final_list[i][14])
	e_flux_lo.append(final_list[i][15])	
'''
#w=open(wd+'new_mosaics_detection/cdwfs_'+band+'_nomultiple.reg','w')
#for i in range(len(newpool)):
#	w.write('circle('+str(newpool[:,0][i])+'d,'+str(newpool[:,1][i])+'d,4\") # color=red\n')
#w.close()

#write catalog
cat=Table([newpool[:,0],newpool[:,1],newpool[:,2],newpool[:,3],newpool[:,4],newpool[:,5],newpool[:,6],newpool[:,7],newpool[:,8],newpool[:,9],newpool[:,10],newpool[:,11],newpool[:,12],newpool[:,13],newpool[:,14],newpool[:,15]],names=('RA','DEC','PROB','AV_R90','TOT','BKG','NET','E_NET_+','E_NET_-','EXP','CR','E_CR_+','E_CR_-','FLUX','E_FLUX_+','E_FLUX_-'))

#cat.write(wd+'new_mosaics_detection/cdwfs_'+band+'_cat1.fits',format='fits',overwrite=True)
cat.write(wd+simfolder+str(sys.argv[2])+'cdwfs_'+band+'_sim_cat1_exp-psf.fits',format='fits',overwrite=True)

print(len(pool), 'in input')
print(len(newpool), 'in output')
print(float(t_out-t_in),' seconds for the match.')


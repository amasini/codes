#match simulated sources detected to input ones, computing the reliability and completeness
# of the sample - WAY FASTER THAN clean_multiple_sources.py... NOTE THAT COULD BE A FACTOR OF 4 FASTER
# WITH SMARTER CODING
import numpy as np
import sys
from astropy.io import fits
from astropy.table import Table
import time

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

for band in ['broad','soft','hard']:

	#take catalog of wavdetect sources (wavdetect, full mosaic 4x4, 5e-5, only bkgmap)
	cat1=fits.open(wd+'new_mosaics_detection/cdwfs_'+band+'_cat0_200113.fits')

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
	print(len(pool),'total sample of the '+band+' band.')

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


	#write catalog
	cat=Table([newpool[:,0],newpool[:,1],newpool[:,2],newpool[:,3],newpool[:,4],newpool[:,5],newpool[:,6],newpool[:,7],newpool[:,8],newpool[:,9],newpool[:,10],newpool[:,11],newpool[:,12],newpool[:,13],newpool[:,14],newpool[:,15]],names=('RA','DEC','PROB','AV_R90','TOT','BKG','NET','E_NET_+','E_NET_-','EXP','CR','E_CR_+','E_CR_-','FLUX','E_FLUX_+','E_FLUX_-'))

	cat.write(wd+'new_mosaics_detection/cdwfs_'+band+'_cat1_200113.fits',format='fits',overwrite=True)

	print(len(pool), 'in input')
	print(len(newpool), 'in output')
	print(float(t_out-t_in),' seconds for the match.')

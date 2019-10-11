# Code to cross-correlate and merge F, S and H catalogs together
import numpy as np
from astropy.io import fits
from astropy.table import Table
import sys
import time
import matplotlib.pyplot as plt

t_start=time.time()

def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*3600*xx)**2 +((pointa[1]-pointb[1])*3600)**2)

wd='/Users/alberto/Desktop/XBOOTES/'

# Open full band catalog (no multiple, output from clean_multiple_sources.py)
fcat=fits.open(wd+'new_mosaics_detection/cdwfs_broad_cat1_exp-psf.fits')
raf=fcat[1].data['RA']
decf=fcat[1].data['DEC']
r90f=fcat[1].data['AV_R90']

filef=fcat[1].data

fcat.close()

w=open(wd+'new_mosaics_detection/cdwfs_broad_cat1_exp-psf.reg','w')
for i in range(len(raf)):
	w.write('circle('+str(raf[i])+'d,'+str(decf[i])+'d,'+str(r90f[i])+'\")\n')
w.close()

# Open soft band catalog (no multiple)
scat=fits.open(wd+'new_mosaics_detection/cdwfs_soft_cat1_exp-psf.fits')
ras=scat[1].data['RA']
decs=scat[1].data['DEC']
r90s=scat[1].data['AV_R90']

files=scat[1].data
scat.close()

w=open(wd+'new_mosaics_detection/cdwfs_soft_cat1_exp-psf.reg','w')
for i in range(len(ras)):
	w.write('circle('+str(ras[i])+'d,'+str(decs[i])+'d,'+str(r90s[i])+'\") #color=red \n')
w.close()

files=np.array(files)

tin=time.time()
fints=[] # List containing the intersection between F and S
delta = 0.007 #(0.028 ~100")
multi_counterpart=0
for linef in range(len(filef)):
	sourcef=[raf[linef],decf[linef]]
	found=0
	
	filt_files=files[(files['RA'] >= raf[linef]-delta) & (files['RA'] <= raf[linef]+delta) & (files['DEC'] >= decf[linef]-delta) & (files['DEC'] <= decf[linef]+delta)]
		
	counterparts,distances=[],[]
	for lines in range(len(filt_files)):
		sources=[filt_files['RA'][lines],filt_files['DEC'][lines]]
		dist = distance(sourcef,sources)
		if dist <= r90f[linef]:
			counterparts.append([filef[linef],filt_files[lines]])
			distances.append(dist)
			found=1
			
	if found == 0:
		fints.append([filef[linef],np.zeros(len(files[0]))])
		
	else:
		if len(counterparts) == 1:
			fints.append(counterparts[0])
			files=np.delete(files,np.where(files['RA']==counterparts[0][1][0]),0)
		else:
			print('Here for source',raf[linef],decf[linef])
			multi_counterpart += 1
			'''
			counterparts=np.array(counterparts)
			distances=np.array(distances)
			sys.exit()
			fints.append(counterparts[distances==np.min(distances)][0])
			files=np.delete(files,np.where(files['RA']==counterparts[0][1][0]),0)
			'''
for lines in range(len(files)):
	fints.append([np.zeros(len(filef[0])),files[lines]])

# Open hard band catalog (no multiple)
hcat=fits.open(wd+'new_mosaics_detection/cdwfs_hard_cat1_exp-psf.fits')
rah=hcat[1].data['RA']
dech=hcat[1].data['DEC']
r90h=hcat[1].data['AV_R90']

fileh=hcat[1].data
fileh=np.array(fileh)

fintsinth=[]
for i in range(len(fints)):
	if fints[i][0][0] != 0: # The source has been detected in the F band
		source_fors=[fints[i][0][0],fints[i][0][1]]
	else: # Else, use S band detection
		source_fors=[fints[i][1][0],fints[i][1][1]]
	
	found=0
	
	filt_fileh=fileh[(fileh['RA'] >= source_fors[0]-delta) & (fileh['RA'] <= source_fors[0]+delta) & (fileh['DEC'] >= source_fors[1]-delta) & (fileh['DEC'] <= source_fors[1]+delta)]
	
	counterparts,distances=[],[]
	for lineh in range(len(filt_fileh)):
		sourceh=[filt_fileh['RA'][lineh],filt_fileh['DEC'][lineh]]
		dist = distance(sourceh,source_fors)
		if dist <= r90h[lineh]:
			counterparts.append([fints[i],filt_fileh[lineh]])
			distances.append(dist)
			found=1
			
	if found==0:
		fintsinth.append([fints[i],np.zeros(len(fileh[0]))])
		
	else:
		if len(counterparts) == 1:
			fintsinth.append(counterparts[0])
			fileh=np.delete(fileh,np.where(fileh['RA']==counterparts[0][1][0]),0)
		else:
			print('Here for source',raf[linef],decf[linef])
			multi_counterpart += 1
			'''
			counterparts=np.array(counterparts)
			distances=np.array(distances)
			sys.exit()
			fints.append(counterparts[distances==np.min(distances)][0])
			files=np.delete(files,np.where(files['RA']==counterparts[0][1][0]),0)
			'''
			
for lineh in range(len(fileh)):
	fintsinth.append([[np.zeros(len(filef[0])),np.zeros(len(files[0]))],fileh[lineh]])

hcat.close()

print('Finished in',round((time.time()-t_start),1),'seconds, with',len(fintsinth),'sources.')
print('I found',multi_counterpart,'multiple counterparts')

out_raf,out_decf,out_probf,out_r90f,out_totf,out_bkgf,out_netf,out_enetf_up,out_enetf_lo,out_expf,out_crf,out_ecrf_up,out_ecrf_lo,out_fluxf,out_efluxf_up,out_efluxf_lo=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
out_ras,out_decs,out_probs,out_r90s,out_tots,out_bkgs,out_nets,out_enets_up,out_enets_lo,out_exps,out_crs,out_ecrs_up,out_ecrs_lo,out_fluxs,out_efluxs_up,out_efluxs_lo=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
out_rah,out_dech,out_probh,out_r90h,out_toth,out_bkgh,out_neth,out_eneth_up,out_eneth_lo,out_exph,out_crh,out_ecrh_up,out_ecrh_lo,out_fluxh,out_efluxh_up,out_efluxh_lo=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
for i in range(len(fintsinth)):
	out_raf.append(fintsinth[i][0][0][0])
	out_decf.append(fintsinth[i][0][0][1])
	out_probf.append(fintsinth[i][0][0][2])
	out_r90f.append(fintsinth[i][0][0][3])
	out_totf.append(fintsinth[i][0][0][4])
	out_bkgf.append(fintsinth[i][0][0][5])
	out_netf.append(fintsinth[i][0][0][6])
	out_enetf_up.append(fintsinth[i][0][0][7])
	out_enetf_lo.append(fintsinth[i][0][0][8])
	out_expf.append(fintsinth[i][0][0][9])
	out_crf.append(fintsinth[i][0][0][10])
	out_ecrf_up.append(fintsinth[i][0][0][11])
	out_ecrf_lo.append(fintsinth[i][0][0][12])
	out_fluxf.append(fintsinth[i][0][0][13])
	out_efluxf_up.append(fintsinth[i][0][0][14])
	out_efluxf_lo.append(fintsinth[i][0][0][15])
	
	out_ras.append(fintsinth[i][0][1][0])
	out_decs.append(fintsinth[i][0][1][1])
	out_probs.append(fintsinth[i][0][1][2])
	out_r90s.append(fintsinth[i][0][1][3])
	out_tots.append(fintsinth[i][0][1][4])
	out_bkgs.append(fintsinth[i][0][1][5])
	out_nets.append(fintsinth[i][0][1][6])
	out_enets_up.append(fintsinth[i][0][1][7])
	out_enets_lo.append(fintsinth[i][0][1][8])
	out_exps.append(fintsinth[i][0][1][9])
	out_crs.append(fintsinth[i][0][1][10])
	out_ecrs_up.append(fintsinth[i][0][1][11])
	out_ecrs_lo.append(fintsinth[i][0][1][12])
	out_fluxs.append(fintsinth[i][0][1][13])
	out_efluxs_up.append(fintsinth[i][0][1][14])
	out_efluxs_lo.append(fintsinth[i][0][1][15])
	
	out_rah.append(fintsinth[i][1][0])
	out_dech.append(fintsinth[i][1][1])
	out_probh.append(fintsinth[i][1][2])
	out_r90h.append(fintsinth[i][1][3])
	out_toth.append(fintsinth[i][1][4])
	out_bkgh.append(fintsinth[i][1][5])
	out_neth.append(fintsinth[i][1][6])
	out_eneth_up.append(fintsinth[i][1][7])
	out_eneth_lo.append(fintsinth[i][1][8])
	out_exph.append(fintsinth[i][1][9])
	out_crh.append(fintsinth[i][1][10])
	out_ecrh_up.append(fintsinth[i][1][11])
	out_ecrh_lo.append(fintsinth[i][1][12])
	out_fluxh.append(fintsinth[i][1][13])
	out_efluxh_up.append(fintsinth[i][1][14])
	out_efluxh_lo.append(fintsinth[i][1][15])

#write catalog
cat=Table([out_raf,out_decf,out_probf,out_r90f,out_totf,out_bkgf,out_netf,out_enetf_up,out_enetf_lo,out_expf,out_crf,out_ecrf_up,out_ecrf_lo,out_fluxf,out_efluxf_up,out_efluxf_lo,out_ras,out_decs,out_probs,out_r90s,out_tots,out_bkgs,out_nets,out_enets_up,out_enets_lo,out_exps,out_crs,out_ecrs_up,out_ecrs_lo,out_fluxs,out_efluxs_up,out_efluxs_lo,out_rah,out_dech,out_probh,out_r90h,out_toth,out_bkgh,out_neth,out_eneth_up,out_eneth_lo,out_exph,out_crh,out_ecrh_up,out_ecrh_lo,out_fluxh,out_efluxh_up,out_efluxh_lo],names=('RA_F','DEC_F','PROB_F','R90_F','TOT_F','BKG_F','NET_F','E_NET_F_+','E_NET_F_-','EXP_F','CR_F','E_CR_F_+','E_CR_F_-','FLUX_F','E_FLUX_F_+','E_FLUX_F_-', 'RA_S','DEC_S','PROB_S','R90_S','TOT_S','BKG_S','NET_S','E_NET_S_+','E_NET_S_-','EXP_S','CR_S','E_CR_S_+','E_CR_S_-','FLUX_S','E_FLUX_S_+','E_FLUX_S_-', 'RA_H','DEC_H','PROB_H','R90_H','TOT_H','BKG_H','NET_H','E_NET_H_+','E_NET_H_-','EXP_H','CR_H','E_CR_H_+','E_CR_H_-','FLUX_H','E_FLUX_H_+','E_FLUX_H_-'))
cat.write(wd+'new_mosaics_detection/cdwfs_merged_cat0_exp-psf.fits',format='fits',overwrite=True)

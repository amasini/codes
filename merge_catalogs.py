# Code to cross-correlate and merge F, S and H catalogs together
import numpy as np
from astropy.io import fits
from astropy.table import Table
import sys
import time

t_start=time.time()

def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*3600*xx)**2 +((pointa[1]-pointb[1])*3600)**2)

wd='/Users/alberto/Desktop/XBOOTES/'

cutf,cuts,cuth=1.4e-2,1e-2,3.5e-3 # These are the probability cuts in F,S,H bands at 97% rel

# Open full band catalog (no multiple, output from clean_multiple_sources.py)
fcat=fits.open(wd+'mosaic_broad_cat1_3.fits')
raf=fcat[1].data['RA']
decf=fcat[1].data['DEC']
r90f=fcat[1].data['AV_R90']

filef=fcat[1].data
fcat.close()

# Open soft band catalog (no multiple)
scat=fits.open(wd+'mosaic_soft_cat1_3.fits')
ras=scat[1].data['RA']
decs=scat[1].data['DEC']
r90s=scat[1].data['AV_R90']

files=scat[1].data

fints=[] # List containing the output file
for linef in range(len(filef)):
	sourcef=[raf[linef],decf[linef]]
	found=0
	
	for lines in range(len(files)):
		sources=[ras[lines],decs[lines]]
		if distance(sourcef,sources) <= r90f[linef]:
			if found==1:
				print('Something wrong with',raf[linef],decf[linef])
				
			else:
				fints.append([filef[linef],files[lines]])
				found=1
			
	if found==0:
		fints.append([filef[linef],np.zeros(len(files[0]))])


for lines in range(len(files)):
	sources=[ras[lines],decs[lines]]
	found=0

	for i in range(len(fints)):
		if fints[i][1][0] == ras[lines]: # The source has been already matched
			found=1
			
	if found==0:
		fints.append([np.zeros(len(filef[0])),files[lines]])

# Open hard band catalog (no multiple)
hcat=fits.open(wd+'mosaic_hard_cat1_3.fits')
rah=hcat[1].data['RA']
dech=hcat[1].data['DEC']
r90h=hcat[1].data['AV_R90']

fileh=hcat[1].data

fintsinth=[]
for i in range(len(fints)):
	if fints[i][0][0] != 0: # The source has been detected in the F band
		source_fors=[fints[i][0][0],fints[i][0][1]]
	else: # Else, use S band detection
		source_fors=[fints[i][1][0],fints[i][1][1]]
	
	found=0
	
	for lineh in range(len(fileh)):
		sourceh=[rah[lineh],dech[lineh]]
		if distance(sourceh,source_fors) <= r90h[lineh]:
			fintsinth.append([fints[i],fileh[lineh]])
			found=1
		
	if found==0:
		fintsinth.append([fints[i],np.zeros(len(fileh[0]))])

for lineh in range(len(fileh)):
	sourceh=[rah[lineh],dech[lineh]]
	found=0

	for i in range(len(fintsinth)):
		if fintsinth[i][1][0] == rah[lineh]: # The source has been already matched
			found=1
			
	if found==0:
		fintsinth.append([[np.zeros(len(filef[0])),np.zeros(len(files[0]))],fileh[lineh]])

print((time.time()-t_start)/60.,'minutes')
print(len(fintsinth))
print(fintsinth[0])
print(fintsinth[9669])
print(fintsinth[9670])
#for i in range(len(fintsinth)):
#	print(i, fintsinth[i][0][0][0])
sys.exit()
print(fintsinth[0][0][0][0])
'''
new=[]
a,b,c=[],[],[]
for i in range(len(fintsinth)):
	a.append(fintsinth[i][0][0])
	b.append(fintsinth[i][0][1])
	c.append(fintsinth[i][1])
print(a[0],b[0],c[0])
'''
out_raf,out_decf,out_probf,out_r90f,out_totf,out_bkgf,out_netf,out_expf,out_crf,out_fluxf=[],[],[],[],[],[],[],[],[],[]
out_ras,out_decs,out_probs,out_r90s,out_tots,out_bkgs,out_nets,out_exps,out_crs,out_fluxs=[],[],[],[],[],[],[],[],[],[]
out_rah,out_dech,out_probh,out_r90h,out_toth,out_bkgh,out_neth,out_exph,out_crh,out_fluxh=[],[],[],[],[],[],[],[],[],[]
for i in range(len(fintsinth)):
	out_raf.append(fintsinth[i][0][0][0])
	out_decf.append(fintsinth[i][0][0][1])
	out_probf.append(fintsinth[i][0][0][2])
	out_r90f.append(fintsinth[i][0][0][3])
	out_totf.append(fintsinth[i][0][0][4])
	out_bkgf.append(fintsinth[i][0][0][5])
	out_netf.append(fintsinth[i][0][0][6])
	out_expf.append(fintsinth[i][0][0][7])
	out_crf.append(fintsinth[i][0][0][8])
	out_fluxf.append(fintsinth[i][0][0][9])
	
	out_ras.append(fintsinth[i][0][1][0])
	out_decs.append(fintsinth[i][0][1][1])
	out_probs.append(fintsinth[i][0][1][2])
	out_r90s.append(fintsinth[i][0][1][3])
	out_tots.append(fintsinth[i][0][1][4])
	out_bkgs.append(fintsinth[i][0][1][5])
	out_nets.append(fintsinth[i][0][1][6])
	out_exps.append(fintsinth[i][0][1][7])
	out_crs.append(fintsinth[i][0][1][8])
	out_fluxs.append(fintsinth[i][0][1][9])
	
	out_rah.append(fintsinth[i][1][0])
	out_dech.append(fintsinth[i][1][1])
	out_probh.append(fintsinth[i][1][2])
	out_r90h.append(fintsinth[i][1][3])
	out_toth.append(fintsinth[i][1][4])
	out_bkgh.append(fintsinth[i][1][5])
	out_neth.append(fintsinth[i][1][6])
	out_exph.append(fintsinth[i][1][7])
	out_crh.append(fintsinth[i][1][8])
	out_fluxh.append(fintsinth[i][1][9])
#write catalog
cat=Table([out_raf,out_decf,out_probf,out_r90f,out_totf,out_bkgf,out_netf,out_expf,out_crf,out_fluxf,out_ras,out_decs,out_probs,out_r90s,out_tots,out_bkgs,out_nets,out_exps,out_crs,out_fluxs,out_rah,out_dech,out_probh,out_r90h,out_toth,out_bkgh,out_neth,out_exph,out_crh,out_fluxh],names=('RA_F','DEC_F','PROB_F','R90_F','TOT_F','BKG_F','NET_F','EXP_F','CR_F','FLUX_F', 'RA_S','DEC_S','PROB_S','R90_S','TOT_S','BKG_S','NET_S','EXP_S','CR_S','FLUX_S', 'RA_H','DEC_H','PROB_H','R90_H','TOT_H','BKG_H','NET_H','EXP_H','CR_H','FLUX_H'))
cat.write(wd+'prova_cdwfs_merged_cat0.fits',format='fits',overwrite=True)

print((time.time()-t_start)/60.,'minutes')
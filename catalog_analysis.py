import numpy as np
from astropy.io import fits
import sys
import matplotlib.pyplot as plt

wd='/Users/alberto/Desktop/XBOOTES/'

cat=fits.open(wd+'prova_cdwfs_merged_cat0.fits')
raf=cat[1].data['RA_F']
fluxf=cat[1].data['FLUX_F']
probf=cat[1].data['PROB_F']
ras=cat[1].data['RA_S']
fluxs=cat[1].data['CR_S']
probs=cat[1].data['PROB_S']
rah=cat[1].data['RA_H']
fluxh=cat[1].data['CR_H']
probh=cat[1].data['PROB_H']

#cutf,cuts,cuth=1.4e-2,1e-2,3.5e-3 # These are the probability cuts in F,S,H bands at 97% rel -> 9240 srcs
cutf,cuts,cuth=2e-4,2e-4,1e-4 # These are the probability cuts in F,S,H bands at 99% rel -> 7666 srcs

probf[raf==0.0]=9999
probs[ras==0.0]=9999
probh[rah==0.0]=9999

count=0
prob,flux=[],[]
prob_gray,flux_gray=[],[]
for j in range(len(probf)):
    if (probf[j] <= cutf or probs[j] <= cuts or probh[j] <= cuth):
        count=count+1

for i in range(len(probf)):
	if probf[i] <= cutf:
		flux.append(fluxf[i])
		prob.append(probf[i])
	else:
		if probf[i] <=1:
			flux_gray.append(fluxf[i])
			prob_gray.append(probf[i])
		

print(count,len(raf))
print(len(prob))
plt.figure()
plt.plot(prob,flux,'k.')
plt.plot(prob_gray,flux_gray,'.',color='gray')
plt.axvline(x=cutf,linestyle='dashed')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Prob of being spurious')
plt.ylabel('0.5-7 keV flux [erg cm-2 s-1]')
plt.axis([1e-30,1,1e-16,1e-12])
plt.tight_layout()
plt.savefig(wd+'cdwfs_full_flux.pdf',format='pdf',dpi=1000)
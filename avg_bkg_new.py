import numpy as np
import sys
import subprocess as s
import os
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from astropy.stats import sigma_clip

wd='/Users/alberto/Desktop/XBOOTES/'

#(old_bkg,old_e_bkg)=np.genfromtxt(wd+'avg_bkg.dat',unpack=True,skip_header=1,usecols=[2,3])
#old_obs=np.genfromtxt(wd+'avg_bkg.dat',unpack=True,skip_header=1,usecols=1,dtype='str')
(mjd,bkg,e_bkg)=np.genfromtxt(wd+'avg_bkg_new.dat',unpack=True,skip_header=1,usecols=[1,2,3])
obs=np.genfromtxt(wd+'avg_bkg_new.dat',unpack=True,skip_header=1,usecols=0,dtype='str')
'''
old_mjd,mjd=[],[]

for i in range(len(old_obs)):
	if len(old_obs[i]) == 4:
		stem='0'+old_obs[i]
	elif len(old_obs[i]) == 3:
		stem='00'+old_obs[i]
	elif len(old_obs[i]) == 5:
		stem=old_obs[i]
		
	aux_mjd=s.check_output('dmkeypar '+wd+'data/'+old_obs[i]+'/repro_new_asol/acisf'+stem+'_repro_evt2.fits MJD_OBS echo=yes',shell=True)
	old_mjd.append(float(aux_mjd))

for j in range(len(obs)):
	if len(obs[j]) == 4:
		stem='0'+obs[j]
	elif len(obs[j]) == 3:
		stem='00'+obs[j]
	elif len(obs[j]) == 5:
		stem=obs[j]
		
	aux_mjd=s.check_output('dmkeypar '+wd+'data/'+obs[j]+'/repro_new_asol/acisf'+stem+'_repro_evt2.fits MJD_OBS echo=yes',shell=True)
	mjd.append(float(aux_mjd))


mjd=np.array(mjd)
#bkg=bkg[mjd>53000]
#mjd=mjd[mjd>53000]
'''

T=11.0
omega=(2*np.pi)/(365*T)
peakday=54800
phi=(np.pi/2)-omega*peakday

fine_t = np.arange(min(mjd),max(mjd),0.1)
my_fit=3e-08*np.sin(omega*fine_t+phi)+1e-7

fig,ax=plt.subplots(1)
#plt.plot(old_mjd,old_bkg,'b.',label='Before')
ax.errorbar(mjd,bkg,yerr=e_bkg,color='red',fmt='.',label='After')
ax.plot(fine_t, my_fit, 'k--',label='Sinusoidal T=11 yrs')
ax.axvline(x=54800,linestyle='dotted',color='k')
ax.annotate('Cycle 24 starts',xy=(54550,2e-7),rotation=90)
ax.axvline(x=56748,linestyle='dotted',color='k')
ax.annotate('Cycle 24 maximum',xy=(56500,2e-7),rotation=90)
ax.set_yscale('log')
ax.set_xlabel('MJD')
ax.set_ylabel(r'Bkg Surf Brightness (cts pixel$^{-2}$ s$^{-1}$)')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.tick_params(axis='both',which='major',direction='in',length=8)
ax.tick_params(axis='y',which='minor',direction='in',length=6)

ax2ticks=[53005,53736,54466,55197,55927,56658,57388,58119]
ax2tickslabels=['2004','2006','2008','2010','2012','2014','2016','2018']
ax2=ax.twiny()
ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(ax2ticks)
ax2.set_xticklabels(ax2tickslabels)
ax2.set_xlabel('Year')
ax2.tick_params(axis='both',which='major',direction='in',length=8)
ax2.tick_params(axis='y',which='minor',direction='in',length=6)
#plt.legend()
plt.tight_layout()
#plt.show()
plt.savefig(wd+'avg_bkg_new.pdf',format='pdf')
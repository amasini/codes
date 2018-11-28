import numpy as np
import subprocess as s
import sys

wd='/Users/alberto/Desktop/XBOOTES/'
(dirs,obsid)=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=[0,1],dtype=str)


for i in range(len(dirs)):
    if len(obsid[i]) == 4:
        stem='0'+obsid[i]
    elif len(obsid[i]) == 3:
        stem='00'+obsid[i]
    elif len(obsid[i]) == 5:
        stem=obsid[i]
    s.call('ds9 '+wd+dirs[i]+'/'+obsid[i]+'/repro_new_asol/expo_broad_flux_4rebinned.img -scale log -cmap sls -zoom 0.6',shell=True)

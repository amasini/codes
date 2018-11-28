### Simulation script for just one obsid of the CDWFS. choose the dir and obsid below

import matplotlib.pyplot as plt
import sys
import subprocess as s
import os
import time
import numpy as np
from astropy.io import fits
from ciao_contrib.region.check_fov import FOVFiles

start_time=time.time()

################
dirs='CDWFS'   #
obsid='19684'  #
################
if len(obsid) == 4:
    stem='0'+obsid
elif len(obsid) == 3:
    stem='00'+obsid
elif len(obsid) == 5:
    stem=obsid

wd='/Users/alberto/Desktop/XBOOTES/'

(f_s,ra_s,dec_s)=np.genfromtxt(wd+'poiss_rand.dat',skip_header=1,unpack=True,usecols=[0,1,2])
ra_s=ra_s[f_s>1e-17]
dec_s=dec_s[f_s>1e-17]
f_s=f_s[f_s>1e-17]
#print(len(ra_s))
dat=fits.open('psf1_rebinned.fits')
data_cts=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=2)
(dirs0,obsid0)=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=[0,1],dtype=str)

datafull=data_cts[obsid0==obsid]
sources_ra,sources_dec,sources_flux=[],[],[]
my_obs = FOVFiles(wd+'data/'+obsid+'/repro_new_asol/fov_acisI.fits')
for k in range(len(ra_s)):
    myobs = my_obs.inside(ra_s[k], dec_s[k])
    if myobs !=[]:
        sources_ra.append(ra_s[k])
        sources_dec.append(dec_s[k])
        sources_flux.append(f_s[k])
print('In obsid '+obsid+' there are '+str(len(sources_ra))+' sources')
                    
bkg=fits.open(wd+'data/'+obsid+'/repro_new_asol/out/acisf'+stem+'_bkgmap.fits')
back=bkg[0].data
backheader=bkg[0].header
bkg.close()

exp=fits.open(wd+'data/'+obsid+'/repro_new_asol/out/acisf'+stem+'_expomap.fits')
expomap=exp[0].data
exp.close()
    
noback=np.zeros_like(back)
    
path=wd+'data/'+obsid+'/repro_new_asol/out/acisf'+stem+'_expomap.fits'
roll=s.check_output('dmkeypar '+path+' ROLL_PNT echo=yes',shell=True) #arcmin
roll=float(roll)
for i in range(len(sources_ra)):
    #dmcoords from cel to msc (ra,dec)->(theta,phi,logicalx,logicaly)
    s.call('punlearn dmcoords',shell=True)
    s.call('dmcoords '+path+' asol=none opt=cel celfmt=deg ra='+str(sources_ra[i])+' dec='+str(sources_dec[i])+'',shell=True)
    res=s.check_output('pget dmcoords theta phi logicalx logicaly',shell=True)
    theta,phi,logicalx,logicaly=res.splitlines()
    theta,phi,logicalx,logicaly=float(theta),float(phi),float(logicalx),float(logicaly)
    elev=theta*np.sin((phi+roll)*np.pi/180.)
    azim=theta*np.cos((phi+roll)*np.pi/180.)
    index0=10+int(round(elev))
    index1=10+int(round(azim))
    if index0 > 20:
        index0=20
    if index1 > 20:
        index1=20
    roundedx=round(logicalx)
    roundedy=round(logicaly)
    #21elevations 21azimuths 5energies 1defocus 64x64 pixels
    psf=dat[0].data[index0][index1][1][0]
    expo=expomap[int(roundedy-1)][int(roundedx-1)]
    #compute total counts and rescale psf; PIMMS predicts 5.438E+10 cps with CHANDRA ACIS-I (0.5-7 keV, Gamma=1.8; would be 5.392E+10 w/ Gamma=1.4)
    counts=sources_flux[i]*expo*5.438E+10
    newpsf=psf*(counts/np.sum(psf))
    index2=int(roundedy-1-31.5)
    index3=int(roundedx-1-31.5)
    index4=int(roundedx-1+31.5)
    index5=int(roundedy-1+31.5)
    c1=[index2,index3]
    c2=[index2,index4]
    c3=[index5,index3]
    c4=[index5,index4]
    #check if the psf image is inside the FoV and adjust accordingly
    nlines=len(noback)
    ncol=len(noback[0])
    if c1[0] < 0:
        c1[0]=0
        c2[0]=0
    if c1[1] < 0:
        c1[1]=0
        c3[1]=0
    if c4[0] > nlines-1:
        c4[0]=nlines-1
        c3[0]=nlines-1
    if c4[1] > ncol-1:
        c4[1]=ncol-1
        c2[1]=ncol-1
                        
    n=0       
    for jj in range(c1[0],c3[0]+1):
        m=0
        for ii in range(c1[1],c2[1]+1):
            noback[jj][ii]=noback[jj][ii]+newpsf[n][m]
            m=m+1
        n=n+1

sim_sources=np.sum(noback)
newback=back*((datafull-sim_sources)/np.sum(back))
simulation=newback+noback
                    
hdu0 = fits.PrimaryHDU(simulation,header=backheader)              
hdu0.writeto(wd+'/sim_full/acisf'+stem+'_sim.fits',overwrite=True) 
                    
poiss_sim=np.random.poisson(simulation)
hdu1 = fits.PrimaryHDU(poiss_sim,header=backheader)
hdu1.writeto(wd+'/sim_full/acisf'+stem+'_sim_poiss.fits',overwrite=True)  
                    
elapsed_time=time.time()-start_time
print(float(elapsed_time)/60.)
dat.close()

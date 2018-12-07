# THIS IS THE SIMULATION SCRIPT - THE STEPS ARE:
# TAKE THE INPUT LIST OF SOURCES AND FOR EACH SOURCE:
# TAkE THE CORRECT PSF IMAGE (TAKING INTO ACCOUNT ROLL ANGLE OF OBSID) AND SUM IT TO BKG
# CREATE A TOTAL IMAGE OF SOURCES + BACKGROUND - A SIMULATION IS THEN THE POISSONIAN 
# REALIZATION OF SUCH IMAGE
# NOTE - THE SIMULATIONS ARE FULL RESOLUTION (0.492"/PIXEL)
import matplotlib.pyplot as plt
import sys
import subprocess as s
import os
import time
import numpy as np
from astropy.io import fits
from ciao_contrib.region.check_fov import FOVFiles
from scipy.ndimage import rotate

def rot(image, xy, angle):
    im_rot = rotate(image,angle) 
    org_center = (np.array(image.shape[:2][::-1])-1)/2.
    rot_center = (np.array(im_rot.shape[:2][::-1])-1)/2.
    org = xy-org_center
    a = np.deg2rad(angle)
    new = np.array([org[0]*np.cos(a) + org[1]*np.sin(a),
            -org[0]*np.sin(a) + org[1]*np.cos(a) ])
    return im_rot, new+rot_center

start_time=time.time()

wd='/Users/alberto/Desktop/XBOOTES/'

#################
# take input sources in full band
(f_s,ra_s,dec_s)=np.genfromtxt(wd+'poiss_rand_lehmerx20.dat',skip_header=1,unpack=True,usecols=[0,1,2])
#################

#################
# take psffile, remember that the energies are 0.277 keV, 1.49 keV, 4.51 keV, 6.4 keV, 8.6 keV
dat=fits.open(wd+'psf1_rebinned.fits')
#################

#################
# take names of directories
obsid0=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=1,dtype=str)
#################

# START THE LOOP ON THE OBSIDS
for j in range(len(obsid0)):
	if len(obsid0[j]) == 4:
		stem='0'+obsid0[j]
	elif len(obsid0[j]) == 3:
		stem='00'+obsid0[j]
	elif len(obsid0[j]) == 5:
		stem=obsid0[j]
	else:
		print('Something\'s wrong with '+obsid0[j]+'/')
		sys.exit()

	if obsid0[j]!='AAA': # this if is here to filter out some obsids, if needed
		sources_ra,sources_dec,sources_flux=[],[],[]
		my_obs = FOVFiles(wd+'data/'+obsid0[j]+'/repro_new_asol/fov_acisI.fits')
		for k in range(len(ra_s)):
			myobs = my_obs.inside(ra_s[k], dec_s[k])
			if myobs !=[]:
				sources_ra.append(ra_s[k])
				sources_dec.append(dec_s[k])
				sources_flux.append(f_s[k])
		print('In obsid '+obsid0[j]+' there are '+str(len(sources_ra))+' sources')
	
		path=wd+'data/'+obsid0[j]+'/repro_new_asol/out/acisf'+stem+'_broad_expomap.fits'
		roll=s.check_output('dmkeypar '+path+' ROLL_PNT echo=yes',shell=True) #arcmin
		roll=float(roll)
		
		#Start loop on bands
		for band in ['broad','soft','hard']:
	
			if band=='broad':
				band2='broad'
				cf=5.392E+10 #conversion factor from flux to count rate from PIMMS (Gamma=1.4)
				flux_ratio=1.
			elif band=='soft':
				band2='0.5-2'
				flux_ratio=0.33 #(Gamma=1.4)
				cf=6.758E+10 #conversion factor from flux to count rate from PIMMS (Gamma=1.4)
			elif band=='hard':
				band2='hard'
				flux_ratio=0.67 #(Gamma=1.4)
				cf=4.721E+10 #conversion factor from flux to count rate from PIMMS (Gamma=1.4)
			
			bkg=fits.open(wd+'data/'+obsid0[j]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap.fits')
			back=bkg[0].data
			backheader=bkg[0].header
			bkg.close()
	
			exp=fits.open(wd+'data/'+obsid0[j]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_expomap.fits')
			expomap=exp[0].data
			exp.close()
	
			noback=np.zeros_like(back)

			print('Throwing ('+band+') sources on bkg in '+obsid0[j]+'...')
		
			for i in range(len(sources_ra)):
				#dmcoords from cel to msc (ra,dec)->(theta,phi,logicalx,logicaly)
				s.call('punlearn dmcoords',shell=True)
				s.call('dmcoords '+path+' asol=none opt=cel celfmt=deg ra='+str(sources_ra[i])+' dec='+str(sources_dec[i])+'',shell=True)
				res=s.check_output('pget dmcoords theta phi logicalx logicaly',shell=True)
				theta,phi,logicalx,logicaly=res.splitlines()
				theta,phi,logicalx,logicaly=float(theta),float(phi),float(logicalx),float(logicaly)
				elev=theta*np.sin(phi*np.pi/180.)
				azim=theta*np.cos(phi*np.pi/180.)
				index0=10+int(round(elev))
				index1=10+int(round(azim))
				if index0 > 20:
					index0=20
				if index1 > 20:
					index1=20
				roundedx=round(logicalx)
				roundedy=round(logicaly)
				#21elevations 21azimuths 5energies 1defocus 64x64 pixels
				if band != 'hard':
					psf=dat[0].data[index0][index1][1][0]
				else:
					psf=dat[0].data[index0][index1][2][0]
				expo=expomap[int(roundedy-1)][int(roundedx-1)]
				#compute total counts and rescale psf; PIMMS predicts 5.438E+10 cps with CHANDRA ACIS-I (0.5-7 keV, Gamma=1.8; would be 5.392E+10 w/ Gamma=1.4)
				counts=sources_flux[i]*expo*cf*flux_ratio
				newpsf=psf*(counts/np.sum(psf))
				
				x0,y0 = 128,128 # (xrot,yrot) should point there
				
				psf_rot = rotate(newpsf,roll) 
				
				index2=int(roundedy-1-((psf_rot.shape[0]/2.)-0.5))
				index3=int(roundedx-1-((psf_rot.shape[0]/2.)-0.5))
				index4=int(roundedx-1+((psf_rot.shape[0]/2.)-0.5))
				index5=int(roundedy-1+((psf_rot.shape[0]/2.)-0.5))
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
						noback[jj][ii]=noback[jj][ii]+psf_rot[n][m]
						m=m+1
					n=n+1
			
			#sim_sources=np.sum(noback)
			#newback=back*((datafull-sim_sources)/np.sum(back)) #rescale background to get same total counts
			
			#print('Original background is',np.sum(back),'Rescaled one is',np.sum(newback))
		
			#simulation=newback+noback
			simulation=back+noback
			
			simulation[simulation<0]=0      
			hdu0 = fits.PrimaryHDU(simulation,header=backheader)                
			hdu0.writeto(wd+'/sim_all/acisf'+stem+'_'+band+'_sim.fits',overwrite=True) 
	
			#old method
			poiss_sim=np.random.poisson(simulation)
			hdu1 = fits.PrimaryHDU(poiss_sim,header=backheader)
			hdu1.writeto(wd+'/sim_all/acisf'+stem+'_'+band+'_sim_poiss.fits',overwrite=True)
			poiss_sim2=poiss_sim.astype(float) #this prevents the simulation to have bitpix=64 causing issues with dmextract
			hdu2 = fits.PrimaryHDU(poiss_sim2,header=backheader)
			hdu2.writeto(wd+'/sim_all/acisf'+stem+'_'+band+'_sim_poiss_bitpix-64.fits',overwrite=True)
	    	            
elapsed_time=time.time()-start_time
print(float(elapsed_time)/3600.,'hours for the whole simulation.')
dat.close()

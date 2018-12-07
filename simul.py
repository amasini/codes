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

################
# define the band
band='soft'
band2='0.5-2' #different for soft band only, used 0.5-2 for bkgmap and expomap, while broad and hard are ok
#################

#################
# take input sources; remember that full band requires lehmerx20.dat
#(f_s,ra_s,dec_s)=np.genfromtxt(wd+'poiss_rand_'+band+'_lehmer.dat',skip_header=1,unpack=True,usecols=[0,1,2])
(f_s,ra_s,dec_s)=np.genfromtxt(wd+'poiss_rand_lehmerx20.dat',skip_header=1,unpack=True,usecols=[0,1,2])
#################

#################
# take psffile, remember that the energies are 0.277 keV, 1.49 keV, 4.51 keV, 6.4 keV, 8.6 keV
dat=fits.open(wd+'psf1_rebinned.fits')
#################

#################
# take counts in the images, based on the band adopted
if band=='broad':
	data_cts=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=2) #dirs obsid cts_full cts_soft cts_hard
	cf=5.438E+10 #conversion factor from flux to count rate from PIMMS (Gamma=1.8)
elif band=='soft':
	data_cts=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=3) #dirs obsid cts_full cts_soft cts_hard
	#cf=6.127E+10 #conversion factor from flux to count rate from PIMMS (Gamma=1.8)
	cf=6.758E+10 #conversion factor from flux to count rate from PIMMS (Gamma=1.4)
elif band=='hard':
	data_cts=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=4) #dirs obsid cts_full cts_soft cts_hard
	cf=4.960E+10 #conversion factor from flux to count rate from PIMMS (Gamma=1.8)
#################

#################
# take names of directories
(dirs0,obsid0)=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=[0,1],dtype=str)
#################

# START THE LOOP ON THE OBSIDS
kk=0
for j in range(len(obsid0)):
	if len(obsid0[j]) == 4:
		stem='0'+obsid0[j]
	elif len(obsid0[j]) == 3:
		stem='00'+obsid0[j]
	elif len(obsid0[j]) == 5:
		stem=obsid0[j]
	else:
		print('Something\'s wrong with '+dirs0[j]+'/'+obsid0[j]+'/')
		sys.exit()

	#produce fov file for ACIS-I
	#filename='acisf'+stem+'_repro_fov1.fits'

	#fovfits=fits.open(wd+'data/'+obsid0[j]+'/repro_new_asol/'+filename)

	#fov=fovfits[1].data
	#fovfits[1].data=fov[fov["CCD_ID"]<4]
	#fovfits.writeto(wd+'data/'+obsid0[j]+'/repro_new_asol/fov_acisI.fits',overwrite=True)
	#fovfits.close()

	if obsid0[j]=='6998': # this if is here to filter out some obsids, if needed
		sources_ra,sources_dec,sources_flux=[],[],[]
		my_obs = FOVFiles(wd+'data/'+obsid0[j]+'/repro_new_asol/fov_acisI.fits')
		for k in range(len(ra_s)):
			myobs = my_obs.inside(ra_s[k], dec_s[k])
			if myobs !=[]:
				sources_ra.append(ra_s[k])
				sources_dec.append(dec_s[k])
				sources_flux.append(f_s[k])
		print('In obsid '+obsid0[j]+' there are '+str(len(sources_ra))+' sources')
	
		### WATCH OUT HERE!!!
		datafull=data_cts[j]
	
		bkg=fits.open(wd+'data/'+obsid0[j]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap.fits')
		back=bkg[0].data
		backheader=bkg[0].header
		bkg.close()
	
		exp=fits.open(wd+'data/'+obsid0[j]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_expomap.fits')
		expomap=exp[0].data
		exp.close()
	
		noback=np.zeros_like(back)
	
		path=wd+'data/'+obsid0[j]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_expomap.fits'
		roll=s.check_output('dmkeypar '+path+' ROLL_PNT echo=yes',shell=True) #arcmin
		roll=float(roll)
		print('Throwing sources on bkg in '+obsid0[j]+'...')
		for i in range(len(sources_ra)):
			#print(i)
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
			counts=sources_flux[i]*expo*cf*2
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
	
		
		sim_sources=np.sum(noback)
		#newsources=noback*((datafull-np.sum(back))/sim_sources) #rescale sources instead of background
		newback=back*((datafull-sim_sources)/np.sum(back)) #rescale background to get same total counts
		
		print('Original background is',np.sum(back),'Rescaled one is',np.sum(newback))
		
		#print(datafull, sim_sources,np.sum(back),np.sum(newback))
		#print(len(newback[newback<0]))
		simulation=newback+noback
		#simulation=back+noback
		#simulation=newsources+back
		
		simulation[simulation<0]=0      
		hdu0 = fits.PrimaryHDU(simulation,header=backheader)                
		hdu0.writeto(wd+'/sim_'+band+'/acisf'+stem+'_sim.fits',overwrite=True) 

		#old method
		poiss_sim=np.random.poisson(simulation)
		hdu1 = fits.PrimaryHDU(poiss_sim,header=backheader)
		hdu1.writeto(wd+'/sim_'+band+'/acisf'+stem+'_sim_poiss.fits',overwrite=True)
		poiss_sim2=poiss_sim.astype(float) #this prevents the simulation to have bitpix=64 causing issues with dmextract
		hdu2 = fits.PrimaryHDU(poiss_sim2,header=backheader)
		hdu2.writeto(wd+'/sim_'+band+'/acisf'+stem+'_sim_poiss_bitpix-64.fits',overwrite=True)
	                
elapsed_time=time.time()-start_time
print(float(elapsed_time)/3600.,'hours for the whole simulation.')
dat.close()
'''
my_obs = FOVFiles(wd+'*/*/repro/fov_acisI.fits')
for i in range(len(ra_s)):
    myobs = my_obs.inside(ra_s[i], dec_s[i])
    if myobs !=[]:
        print(ra_s[i],dec_s[i])
        for j in range(len(myobs)):
            path=myobs[j][:-14]+'out2/*_expomap.fits'
            print(path)
            #dmcoords from cel to msc (ra,dec)->(theta,phi)
            s.call('punlearn dmcoords',shell=True)
            s.call('dmcoords '+path+' asol=none opt=cel celfmt=deg ra='+str(ra_s[i])+' dec='+str(dec_s[i])+'',shell=True)
            theta=s.check_output('pget dmcoords theta',shell=True) #arcmin
            phi=s.check_output('pget dmcoords phi',shell=True) #degrees
            logicalx=s.check_output('pget dmcoords logicalx',shell=True)
            logicaly=s.check_output('pget dmcoords logicaly',shell=True)
            roll=s.check_output('dmkeypar '+path+' ROLL_PNT echo=yes',shell=True) #arcmin
            print(float(theta),float(phi),float(roll))
            elev=float(theta)*np.sin((float(phi)+float(roll))*np.pi/180.)
            azim=float(theta)*np.cos((float(phi)+float(roll))*np.pi/180.)
            index0=10+int(round(elev))
            index1=10+int(round(azim))
            #print(index0,index1)
            if index0 > 20 or index0 < 0:
                print('problem with source in '+path+' at coordinates'+str(ra_s[i])+', '+str(dec_s[i])+'\n')
                sys.exit()
            if index1 > 20 or index1 < 0:
                print('problem with source in '+path+' at coordinates'+str(ra_s[i])+', '+str(dec_s[i])+'\n')
                sys.exit()

            dat=fits.open('psf1.fits')
            #21elevations 21azimuths 5energies 1defocus 256x256 pixels
            psf=dat[0].data[index0][index1][1][0]
            dat.close()
            
            #extract avg exposure in 2"
            s.call('punlearn dmextract',shell=True)
            s.call('pset dmextract infile="'+path+'[bin pos=circle('+str(ra_s[i])+'d,'+str(dec_s[i])+'d,0.0005556d)]" mode=h outfile=expo.fits opt=generic mode=h clobber=yes',shell=True)
            s.call('dmextract',shell=True)
            s.call('dmlist "expo.fits[cols COUNTS, AREA]" data,clean | grep -v COUNTS > expo.dat',shell=True)
            (totexpo,npix)=np.genfromtxt('expo.dat',unpack=True)
            expo=totexpo/npix
            #print(expo,totexpo,npix)
            s.call('rm -f expo.fits',shell=True)
            s.call('rm -f expo.dat',shell=True)

            #compute total counts and rescale psf
            #PIMMS predicts 5.438E+10 cps with CHANDRA ACIS-I (0.5-7 keV, Gamma=1.8-5.392E+10 w/ Gamma=1.4)
            counts=f_s[i]*expo*5.438E+10
            psf=psf*(counts/np.sum(psf))
            hdu = fits.PrimaryHDU(psf)
            hdu.writeto('mypsf.fits',overwrite=True)
            s.call('dmcopy "mypsf.fits[bin x=::4,y=::4]" mypsf_rebinned.fits clobber=yes',shell=True)
            dat2=fits.open('mypsf_rebinned.fits')
            newpsf=dat2[0].data    
            dat2.close()

            name=s.check_output('ls '+myobs[j][:-14]+'out2/*_bkgmap.fits',shell=True)
            bkg=fits.open(name[:-1])
            back=bkg[0].data
            bkg.close()
            #print(logicalx,logicaly,back[int(round(float(logicaly)))-1][int(round(float(logicalx)))-1])
            c1=[int(round(float(logicaly))-1-31.5),int(round(float(logicalx))-1-31.5)]
            c2=[int(round(float(logicaly))-1-31.5),int(round(float(logicalx))-1+31.5)]
            c3=[int(round(float(logicaly))-1+31.5),int(round(float(logicalx))-1-31.5)]
            c4=[int(round(float(logicaly))-1+31.5),int(round(float(logicalx))-1+31.5)]
            n=0       
            for jj in range(c1[0],c3[0]+1):
                m=0
                for ii in range(c1[1],c2[1]+1):
                    back[jj][ii]=back[jj][ii]+newpsf[n][m]
                    #back[jj][ii]=1.0
                    m=m+1
                n=n+1
            #print(logicalx,logicaly)        
            #print(counts,np.sum(newpsf))
            plt.imshow(back,origin='lower')
            plt.show()
            #s.call('rm -f mypsf*.fits',shell=True)
             
        sys.exit()


for dirs in os.listdir(wd):
    if os.path.isdir(dirs) == True:
        for obsid in os.listdir(wd+dirs+'/'):
            if os.path.isdir(wd+dirs+'/'+obsid+'/') == True:
                if len(obsid) == 4:
                    stem='0'+obsid
                elif len(obsid) == 3:
                    stem='00'+obsid
                elif len(obsid) == 5:
                    stem=obsid
                else:
                    print('Something\'s wrong with '+dirs+'/'+obsid+'/')
                    sys.exit()
                sources_ra,sources_dec,sources_flux=[],[],[]
                my_obs = FOVFiles(wd+dirs+'/'+obsid+'/repro/fov_acisI.fits')
                for k in range(len(ra_s)):
                    myobs = my_obs.inside(ra_s[k], dec_s[k])
                    if myobs !=[]:
                        sources_ra.append(ra_s[k])
                        sources_dec.append(dec_s[k])
                        sources_flux.append(f_s[k])
                print('In obsid '+dirs+'/'+obsid+' there are '+str(len(sources_ra))+' sources')

                bkg=fits.open(wd+dirs+'/'+obsid+'/repro/out2/acisf'+stem+'_bkgmap.fits')
                back=bkg[0].data
                back=0.8*back
                backheader=bkg[0].header
                bkg.close()

                exp=fits.open(wd+dirs+'/'+obsid+'/repro/out2/acisf'+stem+'_expomap.fits')
                expomap=exp[0].data
                exp.close()

                if noback == True:
                    back=np.zeros_like(back)

                path=wd+dirs+'/'+obsid+'/repro/out2/acisf'+stem+'_expomap.fits'
                roll=s.check_output('dmkeypar '+path+' ROLL_PNT echo=yes',shell=True) #arcmin

                for i in range(len(sources_ra)):
                    print(i)
                    
                    #dmcoords from cel to msc (ra,dec)->(theta,phi,logicalx,logicaly)
                    s.call('punlearn dmcoords',shell=True)
                    s.call('dmcoords '+path+' asol=none opt=cel celfmt=deg ra='+str(sources_ra[i])+' dec='+str(sources_dec[i])+'',shell=True)
                    res=s.check_output('pget dmcoords theta phi logicalx logicaly',shell=True)
                    theta,phi,logicalx,logicaly=res.splitlines()
                    #theta=s.check_output('pget dmcoords theta',shell=True) #arcmin
                    #phi=s.check_output('pget dmcoords phi',shell=True) #degrees
                    #logicalx=s.check_output('pget dmcoords logicalx',shell=True)
                    #logicaly=s.check_output('pget dmcoords logicaly',shell=True)
                    elev=float(theta)*np.sin((float(phi)+float(roll))*np.pi/180.)
                    azim=float(theta)*np.cos((float(phi)+float(roll))*np.pi/180.)
                    index0=10+int(round(elev))
                    index1=10+int(round(azim))
                    #print(index0,index1)
                                     
                    if index0 > 20 or index0 < 0:
                        print('problem with source in '+path+' at coordinates'+str(sources_ra[i])+', '+str(sources_dec[i])+'\n')
                        sys.exit()
                    if index1 > 20 or index1 < 0:
                        print('problem with source in '+path+' at coordinates'+str(sources_ra[i])+', '+str(sources_dec[i])+'\n')
                        sys.exit()

                    
                    
                    dat=fits.open('psf1.fits')
                    #21elevations 21azimuths 5energies 1defocus 256x256 pixels
                    psf=dat[0].data[index0][index1][1][0]
                    dat.close()
                    
                    #extract avg exposure in 2"
                    s.call('punlearn dmextract',shell=True)
                    s.call('pset dmextract infile="'+path+'[bin pos=circle('+str(sources_ra[i])+'d,'+str(sources_dec[i])+'d,0.0005556d)]" mode=h outfile=expo.fits opt=generic mode=h clobber=yes',shell=True)
                    s.call('dmextract',shell=True)
                    s.call('dmlist "expo.fits[cols COUNTS, AREA]" data,clean | grep -v COUNTS > expo.dat',shell=True)
                    (totexpo,npix)=np.genfromtxt('expo.dat',unpack=True)
                    expo=totexpo/npix
                    #print(expo,totexpo,npix)
                    s.call('rm -f expo.fits',shell=True)
                    s.call('rm -f expo.dat',shell=True)
                    
                    #21elevations 21azimuths 5energies 1defocus 64x64 pixels
                    psf=dat[0].data[index0][index1][1][0]
                    
                    expo=expomap[int(round(float(logicaly))-1)][int(round(float(logicalx))-1)]
                    #compute total counts and rescale psf
                    #PIMMS predicts 5.438E+10 cps with CHANDRA ACIS-I (0.5-7 keV, Gamma=1.8-5.392E+10 w/ Gamma=1.4)
                    counts=sources_flux[i]*expo*5.438E+10
                    newpsf=psf*(counts/np.sum(psf))
                    #psf=psf*(counts/np.sum(psf))
                    #hdu = fits.PrimaryHDU(psf)
                    #hdu.writeto('mypsf.fits',overwrite=True)
                    #s.call('dmcopy "mypsf.fits[bin x=::4,y=::4]" mypsf_rebinned.fits clobber=yes',shell=True)
                    #dat2=fits.open('mypsf_rebinned.fits')
                    #newpsf=dat2[0].data    
                    #dat2.close()

                    #print(logicalx,logicaly,back[int(round(float(logicaly)))-1][int(round(float(logicalx)))-1])
                    c1=[int(round(float(logicaly))-1-31.5),int(round(float(logicalx))-1-31.5)]
                    c2=[int(round(float(logicaly))-1-31.5),int(round(float(logicalx))-1+31.5)]
                    c3=[int(round(float(logicaly))-1+31.5),int(round(float(logicalx))-1-31.5)]
                    c4=[int(round(float(logicaly))-1+31.5),int(round(float(logicalx))-1+31.5)]
                    #check if the psf image is inside the FoV and adjust accordingly
                    if c1[0] < 0:
                        c1[0]=0
                        c2[0]=0
                    if c1[1] < 0:
                        c1[1]=0
                        c3[1]=0
                    if c4[0] > len(back)-1:
                        c4[0]=len(back)-1
                        c3[0]=len(back)-1
                    if c4[1] > len(back[0])-1:
                        c4[1]=len(back[0])-1
                        c2[1]=len(back[0])-1
                    
                    n=0       
                    for jj in range(c1[0],c3[0]+1):
                        m=0
                        for ii in range(c1[1],c2[1]+1):
                            back[jj][ii]=back[jj][ii]+newpsf[n][m]
                            m=m+1
                        n=n+1

                hdu0 = fits.PrimaryHDU(back,header=backheader)                
                poiss_back=np.random.poisson(back)
                hdu1 = fits.PrimaryHDU(poiss_back,header=backheader)
                if noback == True:
                    hdu0.writeto(wd+'/sim_full/acisf'+stem+'_sim_noback.fits',overwrite=True)
                    hdu1.writeto(wd+'/sim_full/acisf'+stem+'_sim_poiss_noback.fits',overwrite=True)  
                else:
                    hdu0.writeto(wd+'/sim_full/acisf'+stem+'_sim.fits',overwrite=True) 
                    hdu1.writeto(wd+'/sim_full/acisf'+stem+'_sim_poiss.fits',overwrite=True)  
                
                elapsed_time=time.time()-start_time
                print(float(elapsed_time)/60.)
                dat.close()
                sys.exit()
'''


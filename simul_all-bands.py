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
import smtplib
from astropy.io import fits
from ciao_contrib.region.check_fov import FOVFiles
from scipy.ndimage import rotate
from email.mime.text import MIMEText

####### FUNCTIONS #######

# A simple Gaussian
def gauss(x,mu,sigma):
	g=np.exp(-(x-mu)**2/(2*sigma**2))
	return g

# Function to rotate an image (for the PSF image)	
def rot(image, xy, angle):
    im_rot = rotate(image,angle) 
    org_center = (np.array(image.shape[:2][::-1])-1)/2.
    rot_center = (np.array(im_rot.shape[:2][::-1])-1)/2.
    org = xy-org_center
    a = np.deg2rad(angle)
    new = np.array([org[0]*np.cos(a) + org[1]*np.sin(a),
            -org[0]*np.sin(a) + org[1]*np.cos(a) ])
    return im_rot, new+rot_center

# Lehmer+12 CDFS dn/ds for AGN
def dnds(fx,band):
	if band == 'broad':
		k,b1,b2,fb,fref = 562.2e14,-1.34,-2.35,8.1e-15,1e-14
	elif band == 'soft':
		k,b1,b2,fb,fref = 169.56e14,-1.49,-2.48,6.0e-15,1e-14
	elif band == 'hard':
		k,b1,b2,fb,fref = 573.13e14,-1.32,-2.55,6.4e-15,1e-14
	else:
		print('Band not recognized. Exit')
		sys.exit()
		
	k1 = k*(fb/fref)**(b1-b2)
	
	if type(fx) == float:
		if fx <= fb:
			return k*(fx/fref)**b1
		else:
			return k1*(fx/fref)**b2
	elif (type(fx) == list)	or (type(fx) == np.ndarray):
		if type(fx) == list:
			fx = np.array(fx)
		aux = k*(fx/fref)**b1
		aux[fx > fb] = k1*(fx[fx > fb]/fref)**b2
		return aux

# Function to generate a new list of input sources
def generate_new_list(band, srcsfilename):

	# Define range of fluxes
	flux=np.logspace(np.log10(5e-17),np.log10(1e-12),101)
	centers0=list((flux[i+1]+flux[i])/2 for i in range(0,len(flux)-1))

	# AGNs per unit flux per squared degree, in soft and hard bands
	dnds0=dnds(centers0,band) # See the relative function
	
	# Write file with sources
	w=open(wd+srcsfilename[:-4]+'.reg','w')
	w2=open(wd+srcsfilename,'w')
	w2.write(band+' flux \t RA \t DEC \t Index \t Gamma\n')
	#choose rectangular area of 4x3.5 deg2 centered on the center of the field
	(minra,maxra)=(215.82,220.1)
	(minde,maxde)=(32.2,36.2)
	area=((maxra-minra)/57.29*(np.sin(maxde/57.29)-np.sin(minde/57.29)))*57.29**2
	
	dn=[]
	for i in range(len(dnds0)):
	
		#sources per square degree in each flux bin
		dn.append(dnds0[i]*(flux[i+1]-flux[i]))

		#sources in total in each flux bin
		n=dnds0[i]*(flux[i+1]-flux[i])*area

		#Poissonian realization of the total number of sources
		N=np.random.poisson(n)

		j=0
		while j < N:
			randdec=np.random.uniform(minde,maxde)
			prob=np.random.uniform(0,1)
			if np.cos(randdec*np.pi/180.) > prob:
				j=j+1
				randra=np.random.uniform(minra,maxra)
				random_index = np.random.choice(xvals, p=new) # Pull a random Gamma for the source
				random_gamma = (random_index+9)/10.
				w.write('circle('+str(randra)+'d,'+str(randdec)+'d,1\")\n')
				w2.write(str(flux[i])+'\t'+str(randra)+'\t'+str(randdec)+'\t'+str(random_index)+'\t'+str(round(random_gamma,1))+'\n')
	w.close()
	w2.close()

# Check the result of the generate_new_list function (put this in the code if needed)
'''
##########################################################################
#recover logn-logs of input sources written in the file
(f_s,ra_s,dec_s)=np.genfromtxt(wd+srcsfilename,skip_header=1,unpack=True,usecols=[0,1,2])
bins00=np.logspace(np.log10(5e-17),np.log10(1e-12),101)
centers00=list((bins00[i+1]+bins00[i])/2 for i in range(0,len(bins00)-1))

quante,bincenters_s=np.histogram(f_s,bins=bins00)

#choose rectangular area of 4x3.5 deg2 centered on the center of the field
(minra,maxra)=(215.82,220.1)
(minde,maxde)=(32.2,36.2)
area=((maxra-minra)/57.29*(np.sin(maxde/57.29)-np.sin(minde/57.29)))*57.29**2
quante_perarea=quante/area # Here, area is larger than 9.3 deg2!

ncum_in=list(reversed(np.cumsum(list(reversed(quante_perarea)))))


#check the cumulatives
dnds0 = dnds(centers00,band)
dn=[]
for kk in range(len(dnds0)):
	dn.append(dnds0[kk]*(bins00[kk+1]-bins00[kk]))
ncum=list(reversed(np.cumsum(list(reversed(dn)))))

#make the plots
f,ax1=plt.subplots(1,1)
ax1.plot(centers00,ncum,'r-',linewidth=2,label='Lehmer+12 0.5-7 keV')
ax1.plot(centers00,ncum_in,'kD',linewidth=2, label='In input file')
ax1.set_xlabel('S [cgs]')
ax1.set_ylabel(r'N(>S) [deg$^-2$]')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.axis([5e-17,1e-13,5,1e5])
ax1.legend()
plt.tight_layout()
plt.show()
sys.exit()
##########################################################################
'''
#########################

start_time=time.time()

wd='/Users/alberto/Desktop/XBOOTES/'

# Define the Gamma PDF (i.e., a random line [0,10] in the ECF files)
xvals = np.arange(0,11,1)
mu=5 # This corresponds to Gamma=1.4
sigma=2
p=[]
for ii in range(len(xvals)):
	p.append(gauss(xvals[ii],mu,sigma))
p=np.array(p)
new=p/np.sum(p)

######## INPUTS #########
# Name of output folder
simfolder='sim_indep_12-Dec-19/'

# Bands to adopt
bands = ['broad','soft','hard']

# PSF file; remember that the energies are 0.277 keV, 1.49 keV, 4.51 keV, 6.4 keV, 8.6 keV
dat=fits.open(wd+'psf1_rebinned.fits')

# Names of directories and their cycles
(obsid0,cy)=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=[1,6],dtype=str)

# Create new list of input sources?
new_list_srcs = True
#########################


#### START THE LOOP ON THE OBSIDS ####

for band in bands:
	if band=='broad':
		band2='broad'
	elif band=='soft':
		band2='0.5-2'
	elif band=='hard':
		band2='hard'
	
	# Generate/Take input sources	
	if new_list_srcs == True:
		srcsfilename = 'poiss_rand_'+band+'_'+simfolder[:-1]+'.dat'
		generate_new_list(band, srcsfilename)		
		(f_s,ra_s,dec_s,ind,gam)=np.genfromtxt(wd+srcsfilename,skip_header=1,unpack=True,usecols=[0,1,2,3,4])
	else:
		srcsfilename = 'poiss_rand_'+band+'_filtered_new.dat'
		(f_s,ra_s,dec_s)=np.genfromtxt(wd+srcsfilename,skip_header=1,unpack=True,usecols=[0,1,2])
	
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
		
		if obsid0[j] != 'AA': # this if is here to filter out some obsids, if needed
			
			# Take the ECFs as a function of Gamma and band (for a given cycle)
			(cf_f,cf_s,cf_h)=np.genfromtxt(wd+'cdwfs_ecf_flux-to-cr_CY'+cy[j]+'.dat',unpack=True,skip_header=1,usecols=[1,2,3])
			cf_f=cf_f*1e10
			cf_s=cf_s*1e10
			cf_h=cf_h*1e10
			
			# Consider only the input sources that fall on the FoV of the considered OBSID
			sources_ra,sources_dec,sources_flux,sources_index,sources_gamma=[],[],[],[],[]
			my_obs = FOVFiles(wd+'data/'+obsid0[j]+'/repro_new_asol/fov_acisI.fits')
			for k in range(len(ra_s)):
				myobs = my_obs.inside(ra_s[k], dec_s[k])
				if myobs !=[]:
					sources_ra.append(ra_s[k])
					sources_dec.append(dec_s[k])
					sources_flux.append(f_s[k])
					sources_index.append(ind[k])
					sources_gamma.append(gam[k])
			print('In obsid '+obsid0[j]+' there are '+str(len(sources_ra))+' sources')
			
			# Extract roll angle (unit degrees)
			path=wd+'data/'+obsid0[j]+'/repro_new_asol/out/acisf'+stem+'_broad_expomap.fits'
			roll=s.check_output('dmkeypar '+path+' ROLL_PNT echo=yes',shell=True) 
			roll=float(roll)
			
			# Take the instrumental background
			bkg=fits.open(wd+'data/'+obsid0[j]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_instr.fits')
			back=bkg[0].data
			backheader=bkg[0].header
			bkg.close()
			# Create an empty space where to sum up sources (same shape of instr bkg)
			noback=np.zeros_like(back)

			print('Throwing ('+band+') sources on bkg in '+obsid0[j]+'...')
			
			# Start putting sources on the background
			for i in range(len(sources_ra)):
				if band=='broad':
					cf=cf_f[int(sources_index[i])] # ECF from flux to count rate given the index pulled
					flux_ratio=1.
				elif band=='soft':
					flux_ratio=1.
					cf=cf_s[int(sources_index[i])]
				elif band=='hard':
					flux_ratio=1.
					cf=cf_h[int(sources_index[i])]
				
				if not cf:
					print("Wrong CF. Exit.")
					print(sources_ra[i],sources_dec[i],sources_flux[i],sources_gamma[i],xvals)
					sys.exit()
				
				# Call dmcoords to transform from cel to msc (ra,dec)->(theta,phi,logicalx,logicaly)
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
				# 21elevations 21azimuths 5energies 1defocus [64x64 pixels]
				if band != 'hard':
					psf=dat[0].data[index0][index1][1][0]
					r90=1.+10.*(theta/10.)**2
				else:
					psf=dat[0].data[index0][index1][2][0]
					r90=1.8+10.*(theta/10.)**2
				
				# Extract the average exposure from the expomap
				s.call('dmextract infile="'+wd+'data/'+obsid0[j]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_expomap.fits[bin pos=circle('+str(sources_ra[i])+'d,'+str(sources_dec[i])+'d,'+str(r90*0.000277778)+'d)]" mode=h outfile=expo.fits opt=generic mode=h clobber=yes',shell=True)
				s.call('dmlist "expo.fits[cols COUNTS, AREA]" data,clean | grep -v COUNTS > expo.dat',shell=True)
				(totexpo,npix)=np.genfromtxt('expo.dat',unpack=True)
				expo=totexpo/npix
				
				# Compute total counts and rescale psf
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
				
				# Check if the psf image is inside the FoV and cut the edges accordingly
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
		
			simulation=back+noback
			
			# If there are negative pixels, put them to zero
			simulation[simulation<0]=0      
			
			# Save the simulation
			hdu0 = fits.PrimaryHDU(simulation,header=backheader)                
			hdu0.writeto(wd+'/'+simfolder+'/acisf'+stem+'_'+band+'_sim.fits',overwrite=True) 
			
			# Rebin the simulation and remove the full-resolution one
			s.call('dmcopy \"'+wd+'/'+simfolder+'/acisf'+stem+'_'+band+'_sim.fits[bin (x,y)=::4]\" outfile='+wd+'/'+simfolder+'/acisf'+stem+'_'+band+'_sim_4reb.fits clobber=yes',shell=True)
			s.call('rm -f '+wd+'/'+simfolder+'/acisf'+stem+'_'+band+'_sim.fits',shell=True)
	    	            
elapsed_time=time.time()-start_time
print(float(elapsed_time)/3600.,'hours for the whole simulation.')
dat.close()

#### SEND AN EMAIL WHEN DONE ####

msg = MIMEText(str((time.time()-start_time)/3600.)+' hours for the simulation.')

me = 'alberto311290@gmail.com'
you = 'alberto.masini@dartmouth.edu'
msg['Subject'] = 'Code is done!'
msg['From'] = me
msg['To'] = you

# Send the message via our own SMTP server, but don't include the
# envelope header.
s0 = smtplib.SMTP('localhost')
s0.sendmail(me, [you], msg.as_string())
s0.quit()
import numpy as np
import sys
import subprocess as s
from astropy.io import fits
from ciao_contrib.region.check_fov import FOVFiles
import matplotlib.pyplot as plt

def distance(pointa, pointb):
	xx = np.cos(pointa[1]/180*3.141592)
	return np.sqrt(((pointa[0]-pointb[0])*xx*3600)**2 +((pointa[1]-pointb[1])*3600)**2)

wd="/Users/alberto/Desktop/XBOOTES/"

band='broad'
band2='broad'

obs=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=1,dtype='str') # List of observations
(cts_f,cts_s,cts_h)=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=[2,3,4],dtype='float') # Total 0.5-7 keV counts

cutf,cuts,cuth=1.4e-2,1e-2,3.5e-3 # These are the probability cuts in F,S,H bands at 97% rel -> 9240 srcs
# Take the catalog of detected sources (cut at 97% reliability)
cat=fits.open('/Users/alberto/Desktop/prova_cdwfs_merged_cat0.fits')
raf=cat[1].data['RA_F']
decf=cat[1].data['DEC_F']
r90f=cat[1].data['R90_F']
probf=cat[1].data['PROB_F']

ras=cat[1].data['RA_S']
decs=cat[1].data['DEC_S']
r90s=cat[1].data['R90_S']
probs=cat[1].data['PROB_S']

rah=cat[1].data['RA_H']
dech=cat[1].data['DEC_H']
r90h=cat[1].data['R90_H']
probh=cat[1].data['PROB_H']

probf[raf==0.0]=9999
probs[ras==0.0]=9999
probh[rah==0.0]=9999

# Take Kenter+05 XBOOTES Clusters
cat0=fits.open(wd+'kenter_clusters.fits')
ra_cl=cat0[1].data['_RAJ2000']
dec_cl=cat0[1].data['_DEJ2000']
r90_cl=cat0[1].data['Size']

ra_src,dec_src,r90_src=[],[],[]
for j in range(len(probf)):
	#if (probf[j] <= cutf or probs[j] <= cuts or probh[j] <= cuth):
	if band=='broad':
		if probf[j] <= cutf: # Select only srcs detected above threshold in F band
			ra_src.append(raf[j])
			dec_src.append(decf[j])
			r90_src.append(r90f[j])
	elif band=='soft':
		if probs[j] <= cuts: # Select only srcs detected above threshold in F band
			ra_src.append(ras[j])
			dec_src.append(decs[j])
			r90_src.append(r90s[j])
	elif band=='hard':
		if probh[j] <= cuth: # Select only srcs detected above threshold in F band
			ra_src.append(rah[j])
			dec_src.append(dech[j])
			r90_src.append(r90h[j])

# Add clusters to total list of sources 
for j in range(len(ra_cl)):
	ra_src.append(ra_cl[j])
	dec_src.append(dec_cl[j])
	r90_src.append(r90_cl[j])

ra_src=np.array(ra_src)
dec_src=np.array(dec_src)
r90_src=np.array(r90_src)

dev=[]
for i in range(len(obs)):
	if len(obs[i]) == 5:
		stem=obs[i]
	elif len(obs[i]) == 4:
		stem='0'+obs[i]
	elif len(obs[i]) == 3:
		stem='00'+obs[i]

	sources_ra,sources_dec,sources_r90,distances=[],[],[],[]
	ra_aim=s.check_output('dmkeypar '+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img RA_PNT echo=yes',shell=True)
	dec_aim=s.check_output('dmkeypar '+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img DEC_PNT echo=yes',shell=True)
	exp=s.check_output('dmkeypar '+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img LIVETIME echo=yes',shell=True)
	ra_aim=float(ra_aim)
	dec_aim=float(dec_aim)
	
	# Create region file with center of field of radius 6'
	w=open(wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_'+band+'_src.reg','w')
	w.write('circle('+str(ra_aim)+'d,'+str(dec_aim)+'d,860\")\n')
	w.close()
	
	my_obs = FOVFiles(wd+'data/'+obs[i]+'/repro_new_asol/fov_acisI.fits')
	
	for k in range(len(ra_src)): # We can substitute this to the direct dmextract to account for extended sources outside of FoV
		myobs = my_obs.inside(ra_src[k], dec_src[k])
		if myobs !=[]:
			sources_ra.append(ra_src[k])
			sources_dec.append(dec_src[k])
			sources_r90.append(r90_src[k])
	
	print('In obsid '+obs[i]+' there are '+str(len(sources_ra))+' sources')
	
	# Create region file with list of sources
	w=open(wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_'+band+'_bkg.reg','w')
	for k in range(len(sources_ra)):
		src=[sources_ra[k],sources_dec[k]]
		aim=[ra_aim,dec_aim]
		if distance(aim,src) <= 3600.:
			radius=2.0*sources_r90[k]
			w.write('circle('+str(sources_ra[k])+'d,'+str(sources_dec[k])+'d,'+str(radius)+'\")\n')
	w.close()
	
	# DMCOPY
	#s.call('dmcopy infile="'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img[exclude sky=region('+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_'+band+'_bkg.reg)]" outfile=excised.fits clobber=yes',shell=True)
	#sys.exit()
	
	s.call('punlearn dmextract',shell=True)
	if band=='broad':
		s.call('dmextract infile="'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img[bin sky=region('+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_'+band+'_src.reg)]" bkg="'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img[bin sky=region('+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_'+band+'_bkg.reg)]" outfile=out2.fits opt=generic mode=h clobber=yes',shell=True)
		#s.call('dmextract infile="'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img[bin sky=region('+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_prova_todelete.reg)]" outfile=out2.fits opt=generic mode=h clobber=yes',shell=True)
		#s.call('dmextract "'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img" out2.fits bkg="'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img[bin sky=region('+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_prova_todelete.reg)]"',shell=True)
	elif band=='soft':
		#s.call('dmextract infile="'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to2keV.img[bin sky=region('+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_'+band+'_src.reg)]" outfile=out2.fits opt=generic mode=h clobber=yes',shell=True)
		s.call('dmextract infile="'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to2keV.img[bin sky=region('+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_'+band+'_src.reg)]" bkg="'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to2keV.img[bin sky=region('+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_'+band+'_bkg.reg)]" outfile=out2.fits opt=generic mode=h clobber=yes',shell=True)
	elif band=='hard':
		#s.call('dmextract infile="'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_2to7keV.img[bin sky=region('+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_'+band+'_src.reg)]" outfile=out2.fits opt=generic mode=h clobber=yes',shell=True)
		s.call('dmextract infile="'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_2to7keV.img[bin sky=region('+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_'+band+'_src.reg)]" bkg="'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_2to7keV.img[bin sky=region('+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_'+band+'_bkg.reg)]" outfile=out2.fits opt=generic mode=h clobber=yes',shell=True)
	src_cts=s.check_output('dmlist "out2.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
	src_area=s.check_output('dmlist "out2.fits[cols AREA]" data,clean | grep -v AREA',shell=True)
	bkg_cts=s.check_output('dmlist "out2.fits[cols BG_COUNTS]" data,clean | grep -v COUNTS',shell=True)
	bkg_area=s.check_output('dmlist "out2.fits[cols BG_AREA]" data,clean | grep -v AREA',shell=True)
	print(src_cts,src_area,bkg_cts,bkg_area)
	cts_new=(float(src_cts)-float(bkg_cts))*(2048**2/(2048**2-float(bkg_area)))
	val=cts_new/(2094**2*float(exp))
	dat=fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/'+band2+'_thresh.expmap')
	expo=dat[0].data
	header=dat[0].header
	new=expo*val
	hdu0 = fits.PrimaryHDU(new,header=header)
	hdu0.writeto('bkgmap_new.fits',overwrite=True)
	#val=(float(src_cts)-float(bkg_cts))/((float(src_area)-float(bkg_area)))
	
	#print(val)
	sys.exit()
	
	dat=fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/'+band2+'_thresh.expmap')
	expo=dat[0].data
	header=dat[0].header
	dat.close()
	#rescale expomap to build bkgmap
	#expo[expo==0]=np.nan
	#new=(expo*val)/expo
	new=expo*val
	#write bkgmap
	hdu0 = fits.PrimaryHDU(new,header=header)
	hdu0.writeto(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_new.fits',overwrite=True)
	
	bkgmap=fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap.fits')
	bkg=bkgmap[0].data
	bkgmap.close()
	oldback=np.sum(bkg)
	
	bkgmap2=fits.open(wd+'data/'+obs[i]+'/repro_new_asol/out/acisf'+stem+'_'+band2+'_bkgmap_new.fits')
	bkg2=bkgmap2[0].data
	bkgmap2.close()
	newback=np.sum(bkg2)
	
	x=(newback-oldback)/oldback
	'''
	if band=='broad':
		x=(newback-oldback)/oldback
	elif band=='soft':
		x=(cts_s[i]-float(src_cts))*(2048.**2/(2048.**2-float(src_area)))
	elif band=='hard':
		x=(cts_h[i]-float(src_cts))*(2048.**2/(2048.**2-float(src_area)))
	'''
	if (x*100.) >= -30:
		dev.append(x*100.)
	print(obs[i],oldback,newback,x*100.)

print(np.mean(dev),np.median(dev))
plt.figure()
plt.hist(dev,bins=20)
plt.show()

'''
for k in range(len(ra_src)):
    # Compute distance array for all sources
    
    distances.append(distance(src,aim))

# Mask sources within 15' from the aimpoint to loop on
distances=np.array(distances)
ra_src2=ra_src[distances <= 900.]
dec_src2=dec_src[distances <= 900.]
r90_src2=r90_src[distances <= 900.]
    
for k in range(len(ra_src2)):
        	#radius=2.0*r90_src2[k]
        	#s.call('dmextract infile="'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img[bin pos=circle('+str(ra_src2[k])+'d,'+str(dec_src2[k])+'d,'+str(radius*0.000277778)+'d)]" outfile=out2.fits opt=generic mode=h clobber=yes',shell=True)
        	#src_cts=s.check_output('dmlist "out2.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
        	#if float(src_cts) > 0.0:
'''

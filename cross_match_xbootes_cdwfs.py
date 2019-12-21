# Compare XBOOTES (Kenter+05) and CDWFS (Masini+20) catalogs - all is done (for now) in the broad band
import numpy as np
import sys
import subprocess as s
from astropy.io import fits
import matplotlib.pyplot as plt
import time
import scipy
import scipy.stats.distributions

def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*xx*3600)**2 +((pointa[1]-pointb[1])*3600)**2)

wd='/Users/alberto/Desktop/XBOOTES/'

# Open the CDWFS catalog
cat=fits.open(wd+'CDWFS_I-Ks-3.6_v2.fits')
ra=cat[1].data['CHA_RA']
dec=cat[1].data['CHA_DEC']
r90=cat[1].data['CHA_R90_F']
#tot=cat[1].data['TOT_F']

# Write out the region file of the CDWFS
#w=open(wd+'new_mosaics_detection/cdwfs_merged_cat1.reg','w')
#for j in range(len(ra)):
#	w.write('circle('+str(ra[j])+'d,'+str(dec[j])+'d,'+str(r90[j])+'\") #width=3 color=red\n')
#w.close()

# Open the XBOOTES catalog
ken=fits.open(wd+'xbootes_kenter+05.fits')
rak=ken[1].data['RAJ2000']
deck=ken[1].data['DEJ2000']
totk=ken[1].data['Fcts']
ctss=ken[1].data['Scts']
ctsh=ken[1].data['Hcts']

band='broad'
path=wd+'new_mosaics_detection/cdwfs_'+band+'_r90_4reb.fits'
ima=fits.open(path)
im=ima[0].data
ima.close()
imagemap=wd+'new_mosaics_detection/cdwfs_'+band+'_4reb.fits'
backmap=wd+'new_mosaics_detection/cdwfs_'+band+'_bkgmap_4reb.fits'

mat=0
cts_old,cts_new,prob,dist,hr,ra_un,dec_un=[],[],[],[],[],[],[]
w=open(wd+'xbootes_unmatched.reg','w')
w1=open(wd+'xbootes_unmatched.dat','w')
w1.write('RA \t DEC \t Kenter_Fcts \t CDWFS_Fcts \t CDWFS_Fbkg \t R90 \t prob \n')
tin=time.time()
for i in range(len(rak)):
	xb=[rak[i],deck[i]]
	found=0
	for j in range(len(ra)):
		my=[ra[j],dec[j]]
		d=distance(my,xb)
		if d <= 1.1*r90[j]:
			if found==0:
				dist.append(d)
				found=1
				mat=mat+1
	
	if found==0:
		w.write('circle('+str(rak[i])+'d,'+str(deck[i])+'d,'+str(totk[i])+'\") #width=4 text={'+str(totk[i])+'}\n')
		#print(rak[i],deck[i],totk[i])
		ra_un.append(rak[i])
		dec_un.append(deck[i])
		cts_old.append(totk[i])
		hr.append((ctsh[i]-ctss[i])/(ctsh[i]+ctss[i]))
		
		# Extract the counts in my CDWFS using R90 and compute prob to be compared with the threshold
		s.call('dmcoords '+path+' asol=none opt=cel celfmt=deg ra='+str(rak[i])+' dec='+str(deck[i])+'',shell=True)
		res=s.check_output('pget dmcoords logicalx logicaly',shell=True)
		logicalx,logicaly=res.splitlines()
		logicalx,logicaly=float(logicalx),float(logicaly)
		av_r90=im[int(round(logicaly)-1),int(round(logicalx)-1)]

		#extract exposure from vignetting-corrected expomap
		#expomap=wd+'new_mosaics_detection/cdwfs_'+band+'_expomap_4reb.fits'
		#s.call('dmextract infile="'+expomap+'[bin pos=circle('+str(src_ra[i])+'d,'+str(src_dec[i])+'d,'+str(av_r90*0.000277778)+'d)]" mode=h outfile=expo.fits opt=generic mode=h clobber=yes',shell=True)
		#s.call('dmlist "expo.fits[cols COUNTS, AREA]" data,clean | grep -v COUNTS > expo.dat',shell=True)
		#(totexpo,npix)=np.genfromtxt('expo.dat',unpack=True)
		#av_exp=16.0*totexpo/npix # npix is the number of NATIVE CHANDRA pixels, so need divide it by 16!

		#extract counts and background using r90
		s.call('dmextract infile="'+imagemap+'[bin pos=circle('+str(rak[i])+'d,'+str(deck[i])+'d,'+str(av_r90*0.000277778)+'d)]" mode=h bkg="'+backmap+'[bin pos=circle('+str(rak[i])+'d,'+str(deck[i])+'d,'+str(av_r90*0.000277778)+'d)]" outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
		cts=s.check_output('dmlist "counts.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
		bkg=s.check_output('dmlist "counts.fits[cols BG_COUNTS]" data,clean | grep -v BG_COUNTS',shell=True)

		cts,bkg=float(cts),float(bkg)
		cts_new.append(cts)
		prob.append(scipy.stats.distributions.poisson.pmf(cts,bkg))
		w1.write(str(rak[i])+' \t '+str(deck[i])+' \t '+str(totk[i])+' \t '+str(cts)+' \t '+str(bkg)+' \t '+str(av_r90)+' \t '+str(prob[-1])+' \n')
		
w.close()
w1.close()
prob=np.array(prob)
print('CDWFS:',len(ra),'XBOOTES:',len(rak))
print(mat,'matches')
print(len(ra_un),'XBOOTES missing in CDWFS.')
print(len(prob[prob<7e-5]),'XBOOTES sources which satisfy reliability cut.')
print((time.time()-tin)/60.,'minutes for the match.')

sys.exit()

(ra,dec,cts_old,cts_new,r90,prob)=np.genfromtxt(wd+'xbootes_unmatched.dat',skip_header=1,unpack=True,usecols=[0,1,2,3,5,6])

bins=np.logspace(np.log10(1e-10),np.log10(1),40)
plt.figure()
plt.hist(prob,bins=bins,histtype='step',color='green',linewidth=3)
plt.axvline(x=7e-5)
plt.xscale('log')
plt.show()

r=ra[prob<7e-5]
d=dec[prob<7e-5]
r9=r90[prob<7e-5]
cts_o=cts_old[prob<7e-5]
cts_n=cts_new[prob<7e-5]
p=prob[prob<7e-5]

#w=open(wd+'xbootes_unmatched_BUTsignificant.reg','w')
#for j in range(len(r)):
	#w.write('circle('+str(r[j])+'d,'+str(d[j])+'d,'+str(r9[j])+'\") #width=4 color=yellow\n')
	#s.call('ds9 '+wd+'new_mosaics_detection/cdwfs_broad_4reb.fits -region '+wd+'new_mosaics_detection/cdwfs_merged_cat1.reg -region '+wd+'xbootes_unmatched_BUTsignificant.reg -crop '+str(r[j])+' '+str(d[j])+' 70 70 wcs fk5 arcsec -zoom to fit -pan to '+str(r[j])+' '+str(d[j])+' wcs fk5 -scale log -scale limits 0 25',shell=True)
#w.close()

detml=-np.log10(p)
plt.figure()
plt.scatter(cts_o,cts_n,s=2*detml,color='green',marker='o')
plt.show()

# Create cutouts with ds9
'''
for i in range(len(ra_un)):
	s.call('ds9 '+wd+'new_mosaics_detection/cdwfs_broad_4reb.fits -region '+wd+'new_mosaics_detection/cdwfs_merged_cat1.reg -region '+wd+'xbootes_unmatched.reg -crop '+str(ra_un[i])+' '+str(dec_un[i])+' 70 70 wcs fk5 arcsec -zoom to fit -pan to '+str(ra_un[i])+' '+str(dec_un[i])+' wcs fk5 -scale log -view colorbar off -saveimage '+wd+'cutouts/xb_'+str(cts[i])+'_'+str(i)+'.png -exit',shell=True)
'''
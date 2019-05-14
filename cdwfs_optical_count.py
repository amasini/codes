# Script to work out the optical counterparts to CDWFS sources: good luck to me!
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import sys
import time
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*xx*3600)**2 +((pointa[1]-pointb[1])*3600)**2)
 
cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)
   
wd='/Users/alberto/Desktop/XBOOTES/'

doo=1.0 # optical-optical matching radius
dxo=5.0 # X-ray-optical matching radius
'''
cat4=fits.open('/Users/alberto/Downloads/NDWFS_I_33_34_cat_m.fits.gz')
ra_ndwfs=cat4[1].data['ALPHA_J2000']
dec_ndwfs=cat4[1].data['DELTA_J2000']
cat4.close()
w=open(wd+'NDWFS_I_33_34_src.reg','w')
for k in range(len(ra_ndwfs)):
	w.write('circle('+str(ra_ndwfs[k])+'d,'+str(dec_ndwfs[k])+'d,'+str(doo)+'\") #color=yellow\n')
w.close()

sys.exit()
'''
tin=time.time()

#OPTICAL MATCH
# Take pure X-ray catalog
cat=fits.open(wd+'cdwfs_merged_cat1.fits')
murray_id=cat[1].data['XB_ID']
fullf=cat[1].data['FLUX_F']
softflux=cat[1].data['FLUX_S']
hardflux=cat[1].data['FLUX_H']
cdwfs_ra=cat[1].data['RA']
cdwfs_dec=cat[1].data['DEC']
cdwfs_r90f=cat[1].data['R90_F']
cdwfs_netf=cat[1].data['NET_F']
cdwfs_r90s=cat[1].data['R90_S']
cdwfs_nets=cat[1].data['NET_S']
cdwfs_r90h=cat[1].data['R90_H']
cdwfs_neth=cat[1].data['NET_H']

# Make sure that each source has net counts and r90
s=cdwfs_netf
s[s<=0]=cdwfs_nets[s<=0]
s[s<=0]=cdwfs_neth[s<=0]
r90=cdwfs_r90f
r90[r90<=0]=cdwfs_r90s[r90<=0]
r90[r90<=0]=cdwfs_r90h[r90<=0]

cdwfs_poserr=r90/np.sqrt(s)


bins=np.logspace(np.log10(min(cdwfs_poserr)),np.log10(max(cdwfs_poserr)),100)
plt.figure()
plt.hist(cdwfs_poserr,bins=bins)
plt.xscale('log')
plt.xlabel('Positional error ["]')
plt.show()
sys.exit()

# Make sure that each source has a full band flux (detected or extrapolated from soft/hard bands)
fluxf=fullf
#print(len(fluxf[fluxf<=0]))

fluxf[fluxf<=0]=2.277*softflux[fluxf<=0]
#print(len(fluxf[fluxf<=0]))

fluxf[fluxf<=0]=1.783*hardflux[fluxf<=0]
#print(len(fluxf[fluxf<=0]))

cdwfs=[]
for i in range(len(cdwfs_ra)):
	cdwfs.append([cdwfs_ra[i],cdwfs_dec[i],fluxf[i]])

cdwfs_ar = np.array(cdwfs)
cdwfs_sorted = cdwfs_ar[cdwfs_ar[:,1].argsort()]

# Check sorting by DEC
#ra=cdwfs_sorted[:,0]
#dec=cdwfs_sorted[:,1]
#plt.figure()
#plt.plot(dec,'k.')
#for j in range(32,36):
#	plt.axhline(y=j)
#plt.show()
#sys.exit()

# Open the FLAMINGOS Ks band catalog - doesn't add anything wrt J band...
(ra_ks,dec_ks,mag_ks)=np.genfromtxt('/Users/alberto/Downloads/bootes_dr1/BOOTES_ks_V1.0.cat',unpack=True,skip_header=73,usecols=[7,8,17])

# Open the Spitzer 3.6 micron catalog (SDWFS)
(ra_36m,dec_36m,mag_36m)=np.genfromtxt('/Users/alberto/Downloads/SDWFS_ch1_stack.v34.txt',unpack=True,skip_header=21,usecols=[0,1,18])

# Open the I band NDWFS catalogs
ra_ndwfs,dec_ndwfs,mag_auto=[],[],[]
tin=time.time()
for j in range(32,36):
	print('Loading NDWFS_I_'+str(j)+'_'+str(j+1)+'_cat_m.fits.gz')
	cat4=fits.open('/Users/alberto/Downloads/NDWFS_I_'+str(j)+'_'+str(j+1)+'_cat_m.fits.gz')
	ra=cat4[1].data['ALPHA_J2000']
	dec=cat4[1].data['DELTA_J2000']
	mag=cat4[1].data['MAG_AUTO']
	for i in range(len(ra)):
		ra_ndwfs.append(ra[i])
		dec_ndwfs.append(dec[i])
		mag_auto.append(mag[i])
	cat4.close()

print('Full NDWFS catalog loaded. Starting the match...')
ra_ndwfs=np.array(ra_ndwfs)
dec_ndwfs=np.array(dec_ndwfs)
mag_auto=np.array(mag_auto)
nmI,nmK,nm36m,qmI,qmK,qm36m=[],[],[],[],[],[]
foundI=0
for i in range(len(cdwfs_sorted)):
	f=0
	delta = 0.028 #(0.028 ~100")
	ra_ndwfs_filt=ra_ndwfs[(ra_ndwfs>=cdwfs_sorted[i][0]-delta) & (ra_ndwfs<=cdwfs_sorted[i][0]+delta) & (dec_ndwfs>=cdwfs_sorted[i][1]-delta) & (dec_ndwfs<=cdwfs_sorted[i][1]+delta)]
	dec_ndwfs_filt=dec_ndwfs[(ra_ndwfs>=cdwfs_sorted[i][0]-delta) & (ra_ndwfs<=cdwfs_sorted[i][0]+delta) & (dec_ndwfs>=cdwfs_sorted[i][1]-delta) & (dec_ndwfs<=cdwfs_sorted[i][1]+delta)]
	mag_auto_filt=mag_auto[(ra_ndwfs>=cdwfs_sorted[i][0]-delta) & (ra_ndwfs<=cdwfs_sorted[i][0]+delta) & (dec_ndwfs>=cdwfs_sorted[i][1]-delta) & (dec_ndwfs<=cdwfs_sorted[i][1]+delta)]
	for k in range(len(ra_ndwfs_filt)):
		ndwfs=[ra_ndwfs_filt[k],dec_ndwfs_filt[k]]
		dist=distance(cdwfs_sorted[i,:-1],ndwfs)
		if (dist >= 5.) and (dist <= 30.):
			if mag_auto_filt[k] < 30.0:
				nmI.append(mag_auto_filt[k])
		elif dist <= 1.0:
			if mag_auto_filt[k] < 30.0:
				qmI.append(mag_auto_filt[k])
				f=1
	if f==1:
		foundI=foundI+1
'''
for j in range(32,36):
	print('Opening NDWFS_I_'+str(j)+'_'+str(j+1)+'_cat_m.fits.gz')
	cat4=fits.open('/Users/alberto/Downloads/NDWFS_I_'+str(j)+'_'+str(j+1)+'_cat_m.fits.gz')
	ra_ndwfs=cat4[1].data['ALPHA_J2000']
	dec_ndwfs=cat4[1].data['DELTA_J2000']
	mag_auto=cat4[1].data['MAG_AUTO']
	cat4.close()
	for i in range(len(cdwfs_sorted)):
		if (cdwfs_sorted[i][1] > j) and (cdwfs_sorted[i][1] < j+1):
			# Insert here the mask to restrict the distance computation?
			delta = 0.028 #(0.028 ~100")
			ra_ndwfs_filt=ra_ndwfs[(ra_ndwfs>=cdwfs_sorted[i][0]-delta) & (ra_ndwfs<=cdwfs_sorted[i][0]+delta) & (dec_ndwfs>=cdwfs_sorted[i][1]-delta) & (dec_ndwfs<=cdwfs_sorted[i][1]+delta)]
			dec_ndwfs_filt=dec_ndwfs[(ra_ndwfs>=cdwfs_sorted[i][0]-delta) & (ra_ndwfs<=cdwfs_sorted[i][0]+delta) & (dec_ndwfs>=cdwfs_sorted[i][1]-delta) & (dec_ndwfs<=cdwfs_sorted[i][1]+delta)]
			mag_auto_filt=mag_auto[(ra_ndwfs>=cdwfs_sorted[i][0]-delta) & (ra_ndwfs<=cdwfs_sorted[i][0]+delta) & (dec_ndwfs>=cdwfs_sorted[i][1]-delta) & (dec_ndwfs<=cdwfs_sorted[i][1]+delta)]
			for k in range(len(ra_ndwfs_filt)):
				ndwfs=[ra_ndwfs_filt[k],dec_ndwfs_filt[k]]
				dist=distance(cdwfs_sorted[i,:-1],ndwfs)
				if (dist >= 5.) and (dist <= 30.):
					if mag_auto_filt[k] < 30.0:
						nmI.append(mag_auto_filt[k])
'''
print('I band completed. Now moving to the NIR-MIR samples...')
foundK,found36mu=0,0
for i in range(len(cdwfs_sorted)):
	f,g=0,0				
	delta = 0.028 #(~100")
	ra_ks_filt=ra_ks[(ra_ks>=cdwfs_sorted[i][0]-delta) & (ra_ks<=cdwfs_sorted[i][0]+delta) & (dec_ks>=cdwfs_sorted[i][1]-delta) & (dec_ks<=cdwfs_sorted[i][1]+delta)]
	dec_ks_filt=dec_ks[(ra_ks>=cdwfs_sorted[i][0]-delta) & (ra_ks<=cdwfs_sorted[i][0]+delta) & (dec_ks>=cdwfs_sorted[i][1]-delta) & (dec_ks<=cdwfs_sorted[i][1]+delta)]
	mag_ks_filt=mag_ks[(ra_ks>=cdwfs_sorted[i][0]-delta) & (ra_ks<=cdwfs_sorted[i][0]+delta) & (dec_ks>=cdwfs_sorted[i][1]-delta) & (dec_ks<=cdwfs_sorted[i][1]+delta)]
	for k in range(len(ra_ks_filt)):
		ndwfs=[ra_ks_filt[k],dec_ks_filt[k]]
		dist=distance(cdwfs_sorted[i,:-1],ndwfs)
		if (dist >= 5.) and (dist <= 30.):
			if mag_ks_filt[k] < 30.0:
				nmK.append(mag_ks_filt[k])
		elif dist <= 1.0:
			if mag_ks_filt[k] < 30.0:
				qmK.append(mag_ks_filt[k])
				f=1
	if f==1:
		foundK=foundK+1
	
	ra_36m_filt=ra_36m[(ra_36m>=cdwfs_sorted[i][0]-delta) & (ra_36m<=cdwfs_sorted[i][0]+delta) & (dec_36m>=cdwfs_sorted[i][1]-delta) & (dec_36m<=cdwfs_sorted[i][1]+delta)]
	dec_36m_filt=dec_36m[(ra_36m>=cdwfs_sorted[i][0]-delta) & (ra_36m<=cdwfs_sorted[i][0]+delta) & (dec_36m>=cdwfs_sorted[i][1]-delta) & (dec_36m<=cdwfs_sorted[i][1]+delta)]
	mag_36m_filt=mag_36m[(ra_36m>=cdwfs_sorted[i][0]-delta) & (ra_36m<=cdwfs_sorted[i][0]+delta) & (dec_36m>=cdwfs_sorted[i][1]-delta) & (dec_36m<=cdwfs_sorted[i][1]+delta)]			
	for j in range(len(ra_36m_filt)):
		ndwfs=[ra_36m_filt[j],dec_36m_filt[j]]
		dist=distance(cdwfs_sorted[i,:-1],ndwfs)
		if (dist >= 5.) and (dist <= 30.):
			if mag_36m_filt[j] < 30.0:
				nm36m.append(mag_36m_filt[j])
		elif dist <= 1.0:
			if mag_36m_filt[j] < 30.0:
				qm36m.append(mag_36m_filt[j])
				g=1
	if g==1:
		found36mu=found36mu+1
cat.close()

print('+'*20)
print('Everything done in',(time.time()-tin)/60.,' minutes')
#print('Redshifts of Murray counterparts:',len(z))
#print('Redshifts of other CDWFS:',len(zcdwfs))
#print('Total redshifts found:',len(z)+len(zcdwfs),'of which '+str(len(redshift_flag[redshift_flag=='specz']))+' are specz and '+str(len(redshift_flag[redshift_flag=='photoz']))+' are photoz')
print('X-ray sources found:',foundI, 'I band,', foundK, 'K band,', found36mu,'3.6 micron.')
#print('Unmatched sources:',len(unmat))
#print('Min - max z:',min(z),max(z))
#print('Blendings:',blend)
print('+'*20)
sys.exit()

'''
#write out n(m) functions
w=open(wd+'LR_nm_I.dat','w')
for k in range(len(nmI)):
	w.write(str(nmI[k])+'\n')
w.close()
w=open(wd+'LR_qm_I.dat','w')
for k in range(len(qmI)):
	w.write(str(qmI[k])+'\n')
w.close()
w=open(wd+'LR_nm_K.dat','w')
for k in range(len(nmK)):
	w.write(str(nmK[k])+'\n')
w.close()
w=open(wd+'LR_qm_K.dat','w')
for k in range(len(qmK)):
	w.write(str(qmK[k])+'\n')
w.close()
w=open(wd+'LR_nm_3.6mu.dat','w')
for k in range(len(nm36m)):
	w.write(str(nm36m[k])+'\n')
w.close()
w=open(wd+'LR_qm_3.6mu.dat','w')
for k in range(len(qm36m)):
	w.write(str(qm36m[k])+'\n')
w.close()
'''
binwidth=0.5 # mag

nbinsI=int((max(nmI)-min(nmI))/binwidth)
qbinsI=int((max(qmI)-min(qmI))/binwidth)
plt.figure()
plt.hist(nmI,bins=nbinsI,histtype='step',linestyle='dashed',color='k')
plt.hist(qmI,bins=qbinsI,histtype='step',color='b')
plt.show()

nbinsK=int((max(nmK)-min(nmK))/binwidth)
qbinsK=int((max(qmK)-min(qmK))/binwidth)
plt.figure()
plt.hist(nmK,bins=nbinsK,histtype='step',linestyle='dashed',color='k')
plt.hist(qmK,bins=qbinsK,histtype='step',color='b')
plt.show()

nbins36m=int((max(nm36m)-min(nm36m))/binwidth)
qbins36m=int((max(qm36m)-min(qm36m))/binwidth)
plt.figure()
plt.hist(nm36m,bins=nbins36m,histtype='step',linestyle='dashed',color='k')
plt.hist(qm36m,bins=qbins36m,histtype='step',color='b')
plt.show()

# Write out region file of unmatched sources and their coordinates/flux
#w=open(wd+'cdwfs_no_optical_match.reg','w')
#w1=open(wd+'cdwfs_no_optical_match.dat','w')
#for k in range(len(ra_unmat)):
#	w.write('circle('+str(ra_unmat[k])+'d,'+str(dec_unmat[k])+'d,'+str(dxo)+'\")\n')
#	w1.write(str(ra_unmat[k])+' \t '+str(dec_unmat[k])+' \t '+str(unmat[k])+'\n')
#w.close()
#w1.close()

#plt.figure()
#plt.hist(distance_to_edges,bins=20)
#plt.show()
#END OF OPTICAL MATCH 
'''

# NIR MATCH
# open the catalog of the 380 sources not detected in NDWFS
(ra_opt_un,dec_opt_un,flux_opt_un)=np.genfromtxt(wd+'cdwfs_no_optical_match.dat',unpack=True)

# Open the FLAMINGOS J band catalog
(ra_ir,dec_ir)=np.genfromtxt('/Users/alberto/Downloads/bootes_dr1/BOOTES_j_V1.0.cat',unpack=True,skip_header=73,usecols=[7,8])

# Open the FLAMINGOS Ks band catalog - doesn't add anything
#(ra_ks,dec_ks)=np.genfromtxt('/Users/alberto/Downloads/bootes_dr1/BOOTES_ks_V1.0.cat',unpack=True,skip_header=73,usecols=[7,8])

# Open the Spitzer 3.6 mum catalog (SDWFS)
(ra_ks,dec_ks)=np.genfromtxt('/Users/alberto/Downloads/SDWFS_ch1_stack.v34.txt',unpack=True,skip_header=21,usecols=[0,1])

ra_un,dec_un,un=[],[],[]
blend=0
for i in range(len(ra_opt_un)):
	source=[ra_opt_un[i],dec_opt_un[i]]
	found=0
	for j in range(len(ra_ir)):
		obj=[ra_ir[j],dec_ir[j]]
		d=distance(source,obj)
		if d <= dxo:
			if found==0:
				found=1
			else:
				found=found+1
				blend=blend+1
				print(str(found)+' J-band counterparts for '+str(source))
	
	if found == 0:
		for k in range(len(ra_ks)):
			obj=[ra_ks[k],dec_ks[k]]
			d=distance(source,obj)
			if d <= dxo:
				if found==0:
					found=1
				else:
					found=found+1
					blend=blend+1
					print(str(found)+' Spitzer counterparts for '+str(source))
				
		if found==0:
			ra_un.append(ra_opt_un[i])
			dec_un.append(dec_opt_un[i])
			un.append(flux_opt_un[i])

print('+'*20)
print('Total sources with no optical counterpart:',len(ra_opt_un))
print('NIR Unmatched sources:',len(un))
print('Blendings:',blend)
print('+'*20)

# Write out region file of unmatched sources and their coordinates/flux
w=open(wd+'cdwfs_no_optical-NIR_match.reg','w')
w1=open(wd+'cdwfs_no_optical-NIR_match.dat','w')
for k in range(len(ra_un)):
	w.write('circle('+str(ra_un[k])+'d,'+str(dec_un[k])+'d,'+str(dxo)+'\")\n')
	w1.write(str(ra_un[k])+' \t '+str(dec_un[k])+' \t '+str(un[k])+'\n')
w.close()
w1.close()
'''

'''
# Take pure X-ray catalog
cat=fits.open(wd+'cdwfs_merged_cat1.fits')
murray_id=cat[1].data['XB_ID']
fullf=cat[1].data['FLUX_F']
softflux=cat[1].data['FLUX_S']
hardflux=cat[1].data['FLUX_H']
cdwfs_ra=cat[1].data['RA']
cdwfs_dec=cat[1].data['DEC']

# Make sure that each source has a full band flux (detected or extrapolated from soft/hard bands)
fluxf=fullf
#print(len(fluxf[fluxf<=0]))

fluxf[fluxf<=0]=2.277*softflux[fluxf<=0]
#print(len(fluxf[fluxf<=0]))

fluxf[fluxf<=0]=1.783*hardflux[fluxf<=0]
#print(len(fluxf[fluxf<=0]))


# Open brand+05 catalog to pull optical counterpart coordinates
cat1=fits.open(wd+'xbootes_brand+05.fits')
id=cat1[1].data['CXOXB']
ra_opt=cat1[1].data['RAo']
dec_opt=cat1[1].data['DEo']
rank=cat1[1].data['Ropt']

# Open kochanek+12 catalog to pull redshifts
cat2=fits.open(wd+'ages_kochanek+12.fits')
ra_koc=cat2[1].data['RAJ2000']
dec_koc=cat2[1].data['DEJ2000']
zspec1=cat2[1].data['z1']
snr1=cat2[1].data['S_N1']
zspec2=cat2[1].data['z2']
snr2=cat2[1].data['S_N2']
zspec3=cat2[1].data['z3']
snr3=cat2[1].data['S_N3']
zphot=cat2[1].data['zph']
qso=cat2[1].data['qso']
xcts=cat2[1].data['Xct']
agn=cat2[1].data['agn']

# Open Chung+14 catalog
cat3=fits.open(wd+'xbootes_chung+14.fits')
ra_chu=cat3[1].data['RAJ2000']
dec_chu=cat3[1].data['DEJ2000']
zsp=cat3[1].data['zsp']
zph=cat3[1].data['z_G_A_'] # Photometric redshift with Galaxy + AGN template

blend=0
tin=time.time()
flux,z,fluxcdwfs,zcdwfs=[],[],[],[]
unmat,ra_unmat,dec_unmat=[],[],[]
redshift_flag=[]
for i in range(len(murray_id)):
	if murray_id[i] != '0': # If CDWFS source has Murray counterpart
		for j in range(len(id)):
			if (murray_id[i] == id[j] and rank[j] ==1): # Find the primary optical counterpart in Brand+05 catalog and use its optical coordinates to match
				rao=ra_opt[j]
				deo=dec_opt[j]
				brand=[rao,deo]
				found=0
				d,qso2,xcts2,agn2,zcont=[],[],[],[],[]
				# Look for counterpart in Kochanek
				for k in range(len(ra_koc)):
					kochanek=[ra_koc[k],dec_koc[k]]
					if distance(brand,kochanek) <= doo:
						if found==0: # Is the first I find:
							found=1
							if (str(snr1[k]) != 'nan' or str(snr2[k]) != 'nan' or str(snr3[k]) != 'nan'): # If counterpart has at least one spec-z
								snr=np.array([snr1[k],snr2[k],snr3[k]])	
								zs=np.array([zspec1[k],zspec2[k],zspec3[k]])
								tryz=zs[snr==max(snr)]
								if tryz > 0:
									z.append(float(tryz))
									flux.append(fluxf[i])
									d.append(distance(brand,kochanek))
									qso2.append(qso[k])
									xcts2.append(xcts[k])
									agn2.append(agn[k])
									zcont.append(float(tryz))
									redshift_flag.append('specz')
								else:
									if str(zphot[k]) != 'nan':
										z.append(float(zphot[k]))
										flux.append(fluxf[i])
										d.append(distance(brand,kochanek))
										qso2.append(qso[k])
										xcts2.append(xcts[k])
										agn2.append(agn[k])
										zcont.append(float(zphot[k]))
										redshift_flag.append('photoz')
							else:
								if str(zphot[k]) != 'nan':
									z.append(float(zphot[k]))
									flux.append(fluxf[i])
									d.append(distance(brand,kochanek))
									qso2.append(qso[k])
									xcts2.append(xcts[k])
									agn2.append(agn[k])
									zcont.append(float(zphot[k]))
									redshift_flag.append('photoz')
						else:
							blend=blend+1
							if (str(snr1[k]) != 'nan' or str(snr2[k]) != 'nan' or str(snr3[k]) != 'nan'):
								snr=np.array([snr1[k],snr2[k],snr3[k]])	
								zs=np.array([zspec1[k],zspec2[k],zspec3[k]])
								tryz=zs[snr==max(snr)]
								if tryz > 0:
									d.append(distance(brand,kochanek))
									qso2.append(qso[k])
									xcts2.append(xcts[k])
									agn2.append(agn[k])
									zcont.append(float(tryz))
								else:
									if str(zphot[k]) != 'nan':
										d.append(distance(brand,kochanek))
										qso2.append(qso[k])
										xcts2.append(xcts[k])
										agn2.append(agn[k])
										zcont.append(float(zphot[k]))
							else:
								if str(zphot[k]) != 'nan':
									d.append(distance(brand,kochanek))
									qso2.append(qso[k])
									xcts2.append(xcts[k])
									agn2.append(agn[k])
									zcont.append(float(zphot[k]))
							
							print('Double optical counterpart for Brand '+str(brand))
							print(d)
							print(qso2, '1 for QSO, 0 for galaxy')
							print(agn2, '1 for Extended AGN')
							print(xcts2,' X-ray counts')
							print(zcont, ' Redshifts')
							
				# If no counterpart is found in Kochanek, try with Chung+14
				if found == 0: 
					for m in range(len(ra_chu)):
						chung=[ra_chu[m],dec_chu[m]]
						if distance(brand,chung) <= doo:
							if found==0: # Is the first I find:
								found=1
								if str(zsp[m]) != 'nan': # use spec-z
									z.append(float(zsp[m]))
									flux.append(fluxf[i])
									redshift_flag.append('specz')

								else: # use photo-z
									z.append(float(zph[m]))
									flux.append(fluxf[i])
									redshift_flag.append('photoz')
							else:
								blend=blend+1
								print('Double optical counterpart in Chung for Brand '+str(brand))

	else: # If CDWFS has NOT a Murray counterpart, match using X-ray position
		rao=cdwfs_ra[i]
		deo=cdwfs_dec[i]
		cdwfs=[rao,deo]
		found=0
		d,qso2,xcts2,agn2,zcont=[],[],[],[],[]
		# Look for counterpart in Kochanek
		for k in range(len(ra_koc)):
			kochanek=[ra_koc[k],dec_koc[k]]
			if distance(cdwfs,kochanek) <= dxo:
				if found==0: # Is the first I find:
					found=1
					if (str(snr1[k]) != 'nan' or str(snr2[k]) != 'nan' or str(snr3[k]) != 'nan'):
						snr=np.array([snr1[k],snr2[k],snr3[k]])	
						zs=np.array([zspec1[k],zspec2[k],zspec3[k]])
						tryz=zs[snr==max(snr)]
						if len(tryz) > 1:
							tryz=tryz[0]
						if tryz > 0:
							zcdwfs.append(float(tryz))
							fluxcdwfs.append(fluxf[i])
							d.append(distance(cdwfs,kochanek))
							qso2.append(qso[k])
							xcts2.append(xcts[k])
							agn2.append(agn[k])
							zcont.append(float(tryz))
							redshift_flag.append('specz')
						else:
							if str(zphot[k]) != 'nan':
								zcdwfs.append(float(zphot[k]))
								fluxcdwfs.append(fluxf[i])
								d.append(distance(cdwfs,kochanek))
								qso2.append(qso[k])
								xcts2.append(xcts[k])
								agn2.append(agn[k])
								zcont.append(float(zphot[k]))
								redshift_flag.append('photoz')
					else:
						if str(zphot[k]) != 'nan':
							zcdwfs.append(float(zphot[k]))
							fluxcdwfs.append(fluxf[i])
							d.append(distance(cdwfs,kochanek))
							qso2.append(qso[k])
							xcts2.append(xcts[k])
							agn2.append(agn[k])
							zcont.append(float(zphot[k]))
							redshift_flag.append('photoz')
				else:
					blend=blend+1
					if (str(snr1[k]) != 'nan' or str(snr2[k]) != 'nan' or str(snr3[k]) != 'nan'):
						snr=np.array([snr1[k],snr2[k],snr3[k]])	
						zs=np.array([zspec1[k],zspec2[k],zspec3[k]])
						tryz=zs[snr==max(snr)]
						if len(tryz) > 1:
							tryz=tryz[0]
						if tryz > 0:
							d.append(distance(cdwfs,kochanek))
							qso2.append(qso[k])
							xcts2.append(xcts[k])
							agn2.append(agn[k])
							zcont.append(float(tryz))
						else:
							if str(zphot[k]) != 'nan':
								d.append(distance(cdwfs,kochanek))
								qso2.append(qso[k])
								xcts2.append(xcts[k])
								agn2.append(agn[k])
								zcont.append(float(zphot[k]))
					else:
						if str(zphot[k]) != 'nan':
							d.append(distance(cdwfs,kochanek))
							qso2.append(qso[k])
							xcts2.append(xcts[k])
							agn2.append(agn[k])
							zcont.append(float(zphot[k]))
					
					print('Double optical counterpart for CDWFS '+str(cdwfs))
					print(d)
					print(qso2, '1 for QSO, 0 for galaxy')
					print(agn2, '1 for Extended AGN')
					print(xcts2,' X-ray counts')
					print(zcont, ' Redshifts')
					
		# If no counterpart in Kochanek, try with Chung
		if found == 0: 
			for m in range(len(ra_chu)):
				chung=[ra_chu[m],dec_chu[m]]
				if distance(cdwfs,chung) <= dxo:
					if found==0: # Is the first I find:
						found=1
						if str(zsp[m]) != 'nan': # use spec-z
							zcdwfs.append(float(zsp[m]))
							fluxcdwfs.append(fluxf[i])
							redshift_flag.append('specz')

						else: # use photo-z
							zcdwfs.append(float(zph[m]))
							fluxcdwfs.append(fluxf[i])
							redshift_flag.append('photoz')
					else:
						blend=blend+1
						print('Double optical counterpart in Chung for CDWFS '+str(cdwfs))
	
	if found == 0: # if still no source, try with the full NDWFS catalog
		cat4=fits.open('/Users/alberto/Downloads/NDWFS_I_'+str(int(cdwfs_dec[i]))+'_'+str(int(cdwfs_dec[i])+1)+'_cat_m.fits.gz')
		ra_ndwfs=cat4[1].data['ALPHA_J2000']
		dec_ndwfs=cat4[1].data['DELTA_J2000']
		cat4.close()
		for k in range(len(ra_ndwfs)):
			ndwfs=[ra_ndwfs[k],dec_ndwfs[k]]
			if distance(cdwfs,ndwfs) <= dxo:
				if found==0: # Is the first I find:
					found=1
					counter=1
				else:
					blend=blend+1
					counter=counter+1
					print(str(counter)+' optical counterpart for NDWFS '+str(cdwfs))
	
	if found == 0:
		unmat.append(fluxf[i])
		ra_unmat.append(cdwfs_ra[i])
		dec_unmat.append(cdwfs_dec[i])
		edge=[cdwfs_ra[i],round(cdwfs_dec[i])]
		distance_to_edges.append(distance(cdwfs,edge))

cat.close()
cat1.close()
cat2.close()
cat3.close()

flux=np.array(flux)
z=np.array(z)

zcdwfs=np.array(zcdwfs)
fluxcdwfs=np.array(fluxcdwfs)
redshift_flag=np.array(redshift_flag)

# Write out region file of unmatched sources and their coordinates/flux
w=open(wd+'cdwfs_no_optical_match.reg','w')
w1=open(wd+'cdwfs_no_optical_match.dat','w')
for k in range(len(ra_unmat)):
	w.write('circle('+str(ra_unmat[k])+'d,'+str(dec_unmat[k])+'d,'+str(dxo)+'\")\n')
	w1.write(str(ra_unmat[k])+' \t '+str(dec_unmat[k])+' \t '+str(unmat[k])+'\n')
w.close()
w1.close()

print('+'*20)
print((time.time()-tin)/60.,' minutes')
print('Redshifts of Murray counterparts:',len(z))
print('Redshifts of other CDWFS:',len(zcdwfs))
print('Total redshifts found:',len(z)+len(zcdwfs),'of which '+str(len(redshift_flag[redshift_flag=='specz']))+' are specz and '+str(len(redshift_flag[redshift_flag=='photoz']))+' are photoz')
print('Total sources in catalog:',len(murray_id[murray_id!='0'])+len(murray_id[murray_id=='0']))
print('Unmatched sources:',len(unmat))
print('Min - max z:',min(z),max(z))
print('Blendings:',blend)
print('+'*20)


plt.figure()
plt.hist(distance_to_edges,bins=20)
plt.show()

fluxrestframe=flux*(1+z)**(-0.2)
dl=cosmo.luminosity_distance(z)
dl2=dl.value*3.086e24
lfull=4*3.141592*fluxrestframe*dl2**2

fluxrestframe_cdwfs=fluxcdwfs*(1+zcdwfs)**(-0.2)
dl=cosmo.luminosity_distance(zcdwfs)
dl2=dl.value*3.086e24
lfullcdwfs=4*3.141592*fluxrestframe_cdwfs*dl2**2

# Write out data so plots are faster to be made
w=open(wd+'cdwfs_lxz_0.dat','w')
for i in range(len(z)):
	w.write(''+str(z[i])+' \t '+str(lfull[i])+'\n')
w.close()
w=open(wd+'cdwfs_lxz_1.dat','w')
for i in range(len(zcdwfs)):
	w.write(''+str(zcdwfs[i])+' \t '+str(lfullcdwfs[i])+'\n')
w.close()

zz=np.logspace(np.log10(1e-4),np.log10(5),100)
flim=1e-15
flimrestframe=flim*((1+zz)**(-0.2))
dl=cosmo.luminosity_distance(zz)
dl2=dl.value*3.086e24
l=flimrestframe*4*3.141592*dl2**2

(z,lfull)=np.genfromtxt(wd+'cdwfs_lxz_0.dat',unpack=True)
(zcdwfs,lfullcdwfs)=np.genfromtxt(wd+'cdwfs_lxz_1.dat',unpack=True)

totz=list(z)+list(zcdwfs)
totflux=list(flux)+list(fluxcdwfs)

bins0=np.logspace(np.log10(np.min(totflux)),np.log10(np.max(totflux)),15)
bins=np.logspace(np.log10(np.min(unmat)),np.log10(np.max(unmat)),15)
plt.figure()
plt.hist(totflux,bins=bins0,label='Total CDWFS sources')
plt.hist(unmat,bins=bins,label='No optical counterpart')
plt.xscale('log')
plt.xlabel('Full band Chandra flux [cgs]')
plt.ylabel('N')
plt.legend()
plt.tight_layout()
plt.savefig(wd+'cdwfs_unmatched-flux.pdf',format='pdf')
#plt.show()

#bins=np.logspace(np.log10(min(z)),np.log10(max(z)),30)
plt.figure()
plt.hist(totz,bins=30)
#plt.xscale('log')
plt.xlabel('z',fontsize=20)
plt.ylabel('N',fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=15)
plt.tight_layout()
plt.savefig(wd+'cdwfs_zdistr.pdf',format='pdf')
#plt.show()

plt.figure()
plt.plot(z,lfull,'b.',label='XBOOTES')
plt.plot(zcdwfs,lfullcdwfs,'r.',label='CDWFS')
plt.plot(zz,l,'k--')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('z',fontsize=20)
plt.ylabel(r'$L_{0.5-7}$ [erg/s]',fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=15)
plt.axis([0.01,5,5e39,1e46])
plt.annotate(r'$N=$'+str(len(totz)),xy=(0.02,5e44))
plt.legend()
plt.tight_layout()
plt.savefig(wd+'cdwfs_lxz.pdf',format='pdf')
#plt.show()
'''
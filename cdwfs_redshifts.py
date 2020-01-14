# Code to work out redshifts after the optical-NIR counterparts have been found with NWAY
# Match with Chung+14
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table,Column
import sys
import time

def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*xx*3600)**2 +((pointa[1]-pointb[1])*3600)**2)
    
wd='/Users/alberto/Desktop/XBOOTES/'

'''
imag=[]
for i in range(32,36):
	cat=fits.open('/Users/alberto/Downloads/nway-master/nway_cdwfs-ndwfs/cdwfs_NDWFS_I_'+str(i)+'_'+str(i+1)+'.fits')
	matchcol=cat[1].data['match_flag']
	pany=cat[1].data['p_any']
	ra=cat[1].data['CDWFS_RA']
	dec=cat[1].data['CDWFS_DEC']
	fluxs=cat[1].data['CDWFS_FLUX_S']
	fluxh=cat[1].data['CDWFS_FLUX_H']
	imag0=cat[1].data['NDWFS_MAG_AUTO']
	#zsp=cat[1].data['AGES_z1']
	#zph=cat[1].data['AGES_zph']
	sep=cat[1].data['Separation_NDWFS_CDWFS']
	cat.close()
	
	if i ==32:
		p_any_cut=0.00
	elif i==33:	
		p_any_cut=0.05
	elif i==34:
		p_any_cut=0.12
	elif i==35:
		p_any_cut=0.02
	imag1=imag0[(matchcol==1) & (pany>p_any_cut)]
	for k in range(len(imag1)):
		imag.append(imag1[k])
'''

# Open the matched master catalog (should contain only one CDWFS source per line)  
cat=fits.open(wd+'/nway-master/200109_trial.fits')

p_any_cut=0.31 # 0.54 (0.12 as of 4-Sep-19, 0.20 as of 15-Oct-19, 0.16 as of 30-Nov-19, 0.13 as of 01-Dec-19, 0.31 as of 09-Jan-20) needed to have <5% false associations

#cat=fits.open(wd+'CDWFS_I-Ks-3.6_v2.fits')
data=cat[1].data
cols=cat[1].columns
names=cols.names
pany=data['p_any']
fluxs=data['CDWFS_FLUX_S']
fluxh=data['CDWFS_FLUX_H']
hr0=data['CDWFS_HR']
ehr0p=data['CDWFS_E_HR_+']
ehr0n=data['CDWFS_E_HR_-']
imag0=data['NDWFS_MAG_AUTO']
kmag0=data['IBIS_MAG_BEST']
spmag0=data['SDWFS_ch1_ma']
sep0=data['Separation_NDWFS_CDWFS']
sep0k=data['Separation_IBIS_CDWFS']
sep0sp=data['Separation_SDWFS_CDWFS']
#z=data['z']
#zflag=data['zflag']
cat.close()
print('='*30)
print('The full CDWFS catalog has',len(data),'sources.')

#outcat=[]
#for i in range(len(data)):
#	outcat.append(data[i])
#outcat=np.array(outcat)

imag=imag0[pany>p_any_cut]
kmag=kmag0[pany>p_any_cut]
spmag=spmag0[pany>p_any_cut]
sep=sep0[pany>p_any_cut]
sepk=sep0k[pany>p_any_cut]
sepsp=sep0sp[pany>p_any_cut]
#z=z[pany>p_any_cut]
#zflag=zflag[pany>p_any_cut]
hr1=hr0[pany>p_any_cut]
ehr1p=ehr0p[pany>p_any_cut]
ehr1n=ehr0n[pany>p_any_cut]
#Imag_chu=imag[zflag!=-99]
#Imag_und=imag[zflag==-99]
#Imag_chu_zsp=imag[zflag==1]
#Imag_chu_zph=imag[zflag==0]
print('A total of',len(pany[pany>0]),' (p_any>0) opt-NIR associations found.')
print('A total of',len(imag),'ROBUST (p_any>'+str(p_any_cut)+') opt-NIR associations found.')

sep=sep[(imag<99.0) & (imag!=-99)]
imag=imag[(imag<99.0) & (imag!=-99)]
sepk=sepk[(kmag<99.0) & (kmag!=-99)]
kmag=kmag[(kmag<99.0) & (kmag!=-99)]
sepsp=sepsp[(spmag<99.0) & (spmag!=-99)]
spmag=spmag[(spmag<99.0) & (spmag!=-99)]
print('I band good photometry:',len(imag))
print('Ks band good photometry:',len(kmag))
print('[3.6] band good photometry:',len(spmag))

'''
hr=hr1[zflag!=-99]
ehrp=ehr1p[zflag!=-99]
ehrn=ehr1n[zflag!=-99]
zsp=z[zflag==1]
zph=z[zflag==0]
ztot=list(zsp)+list(zph)
print(len(ztot),' redshifts:',len(zsp),'spec-z,',len(zph),'photo-z in the robust subsample.')
print('='*30)
bins=np.linspace(0,5,30)
plt.figure()
plt.hist(ztot,bins=bins,histtype='step',linewidth=3,color='k')
plt.hist(zsp,bins=bins,histtype='step',linewidth=2,color='cyan',label='Spec-z')
plt.hist(zph,bins=bins,histtype='step',linewidth=2,color='red',label='Photo-z')
plt.legend()
plt.show()

plt.figure()
plt.scatter(ztot,hr,marker='.',color='black')
plt.xlabel(r'$z$')
plt.ylabel(r'HR=(H-S)/(H+S)')
plt.show()
sys.exit()
'''
#bins=np.linspace(np.min(Imag[Imag>0.]),np.max(Imag[Imag>0.]),31)
bins=30
plt.figure()
plt.hist(imag,bins=bins,histtype='step',color='green',linewidth=3,label='I band')
plt.hist(kmag,bins=bins,histtype='step',color='orange',linewidth=3,label='Ks band')
plt.hist(spmag,bins=bins,histtype='step',color='red',linewidth=3,label='[3.6] band')
plt.xlabel(r'Magnitude')
plt.legend(loc='upper left')
plt.tight_layout()
plt.show()

plt.figure()
plt.hist(sep,bins=bins,histtype='step',color='g',linewidth=3,label='I band')
plt.hist(sepk,bins=bins,histtype='step',color='orange',linewidth=3,label='Ks band')
plt.hist(sepsp,bins=bins,histtype='step',color='r',linewidth=3,label='[3.6] band')
plt.xlabel(r'Separation ["]')
plt.legend()
plt.tight_layout()
plt.show()

f,ax=plt.subplots(1,3,figsize=[11,5],sharey=True)
ax[0].scatter(sep,imag,color='green',marker='.')
ax[0].set_xlabel('Separation ["]')
ax[0].set_ylabel('I-band mag')
ax[0].axis([0,15,8,30])

ax[1].scatter(sepk,kmag,color='orange',marker='.')
ax[1].set_xlabel('Separation ["]')
ax[1].set_ylabel('K-band mag')

ax[2].scatter(sepsp,spmag,color='red',marker='.')
ax[2].set_xlabel('Separation ["]')
ax[2].set_ylabel('[3.6]-band mag')

plt.show()

'''++++++++++++++++++++++ REDSHIFT PART ++++++++++++++++++++++'''
'''+++++++++++++++++++ MATCH with CHUNG+14 +++++++++++++++++++'''
opt_ra=data['NDWFS_RA_J2000']
opt_dec=data['NDWFS_DEC_J2000']
Imag=data['NDWFS_MAG_AUTO']

nir_ra=data['IBIS_RA_PEAK_J2000']
nir_dec=data['IBIS_DEC_PEAK_J2000']

ir_ra=data['SDWFS_ra']
ir_dec=data['SDWFS_dec']

ra=opt_ra
ra[opt_ra==-99.0]=ir_ra[opt_ra==-99.0]
ra[ra==-99.0]=nir_ra[ra==-99.0]

dec=opt_dec
dec[opt_dec==-99.0]=ir_dec[opt_dec==-99.0]
dec[dec==-99.0]=nir_dec[dec==-99.0]

# Open chung catalog to get redshifts
chucat=fits.open(wd+'xbootes_chung+14.fits')
chu_ra=chucat[1].data['RAJ2000']
chu_dec=chucat[1].data['DEJ2000']
chu_zsp=chucat[1].data['zsp']
chu_zph_g=chucat[1].data['z_G_']
chu_chi_g=chucat[1].data['chi2_G_']
chu_zph_ga=chucat[1].data['z_G_A_']
chu_chi_ga=chucat[1].data['chi2_G_A_']
chu_ebv=chucat[1].data['E_B-V_']
chu_chi_s=chucat[1].data['chi2_S_']
chucat.close()

#z=chu_zsp
spectro=np.full_like(chu_zsp,True,dtype='bool')
spectro[np.isnan(chu_zsp)]=False
#z[np.isnan(chu_zsp)]=chu_zph[np.isnan(chu_zsp)]

tin=time.time()
mat=0
Imag_chu,Imag_und,sp,z_out_sp,z_out_ga,z_out_g,z_out_s,zflag,ebv,chi_ga,chi_g,chi_s=[],[],[],[],[],[],[],[],[],[],[],[]
for i in range(len(ra)):
	found=0
	opt=[ra[i],dec[i]]
	delta = 0.00056 #(0.028 ~100")
	
	cut_chu_ra=chu_ra[(chu_ra>=ra[i]-delta) & (chu_ra<=ra[i]+delta) & (chu_dec>=dec[i]-delta) & (chu_dec<=dec[i]+delta)]
	cut_chu_dec=chu_dec[(chu_ra>=ra[i]-delta) & (chu_ra<=ra[i]+delta) & (chu_dec>=dec[i]-delta) & (chu_dec<=dec[i]+delta)]
	
	cut_spectro=spectro[(chu_ra>=ra[i]-delta) & (chu_ra<=ra[i]+delta) & (chu_dec>=dec[i]-delta) & (chu_dec<=dec[i]+delta)]
	cut_zsp=chu_zsp[(chu_ra>=ra[i]-delta) & (chu_ra<=ra[i]+delta) & (chu_dec>=dec[i]-delta) & (chu_dec<=dec[i]+delta)]
	cut_zph_ga=chu_zph_ga[(chu_ra>=ra[i]-delta) & (chu_ra<=ra[i]+delta) & (chu_dec>=dec[i]-delta) & (chu_dec<=dec[i]+delta)]
	cut_zph_g=chu_zph_g[(chu_ra>=ra[i]-delta) & (chu_ra<=ra[i]+delta) & (chu_dec>=dec[i]-delta) & (chu_dec<=dec[i]+delta)]
	
	cut_chi_ga=chu_chi_ga[(chu_ra>=ra[i]-delta) & (chu_ra<=ra[i]+delta) & (chu_dec>=dec[i]-delta) & (chu_dec<=dec[i]+delta)]
	cut_chi_g=chu_chi_g[(chu_ra>=ra[i]-delta) & (chu_ra<=ra[i]+delta) & (chu_dec>=dec[i]-delta) & (chu_dec<=dec[i]+delta)]
	cut_chi_s=chu_chi_s[(chu_ra>=ra[i]-delta) & (chu_ra<=ra[i]+delta) & (chu_dec>=dec[i]-delta) & (chu_dec<=dec[i]+delta)]
	cut_ebv=chu_ebv[(chu_ra>=ra[i]-delta) & (chu_ra<=ra[i]+delta) & (chu_dec>=dec[i]-delta) & (chu_dec<=dec[i]+delta)]
	
	for j in range(len(cut_chu_ra)):
		chu=[cut_chu_ra[j],cut_chu_dec[j]]
		d=distance(opt,chu)
		if (d < 0.5) and (found==0):
			found=1
			mat=mat+1
			Imag_chu.append(Imag[i])
			
			if cut_spectro[j] == True:
				sp.append(True)
				#zflag.append(1)
				z_out_sp.append(cut_zsp[j])
				z_out_ga.append(cut_zph_ga[j])
				z_out_g.append(cut_zph_g[j])
				ebv.append(cut_ebv[j])
				chi_ga.append(cut_chi_ga[j])
				chi_g.append(cut_chi_g[j])
				chi_s.append(cut_chi_s[j])
			else:
				sp.append(False)
				#zflag.append(0)
				z_out_sp.append(-99.0)
				z_out_ga.append(cut_zph_ga[j])
				ebv.append(cut_ebv[j])
				chi_ga.append(cut_chi_ga[j])
				z_out_g.append(cut_zph_g[j])
				chi_g.append(cut_chi_g[j])
				chi_s.append(cut_chi_s[j])
	if found==0:
		Imag_und.append(Imag[i])
		#zflag.append(-99.0)
		z_out_sp.append(-99.0)
		z_out_ga.append(-99.0)
		ebv.append(-99.0)
		chi_ga.append(-99.0)
		z_out_g.append(-99.0)
		chi_g.append(-99.0)
		chi_s.append(-99.0)
print(mat,'matched with Chung in',(time.time()-tin)/60.,'minutes.')

# Write out catalog
table = Table.read(wd+'nway-master/200109_trial.fits', format='fits')

t1 = Column(name='zsp', data=z_out_sp)
t2 = Column(name='zph_G+A', data=z_out_ga)
t3 = Column(name='E_B-V', data=ebv)
t4 = Column(name='chi2_G+A', data=chi_ga)
t5 = Column(name='zph_G', data=z_out_g)
t6 = Column(name='chi2_G', data=chi_g)
t7 = Column(name='chi2_S', data=chi_s)
table.add_column(t1)
table.add_column(t2)
table.add_column(t3)
table.add_column(t4)
table.add_column(t5)
table.add_column(t6)
table.add_column(t7)
table.write(wd+'CDWFS_I-Ks-3.6_v200109.fits', format='fits', overwrite='True')
sys.exit()
'''
vet_col=[]
for i in range(len(names)):
	vet_col.append(outcat[:,i])

names.append('z')
names.append('zflag')

vet_col.append(z_out)
vet_col.append(zflag)

cat=Table(vet_col,names=names)
cat.write(wd+'CDWFS_I-Ks-3.6_v2.fits',format='fits',overwrite=True)
'''

bins=np.linspace(12,max(Imag_chu[(Imag_chu<99.) & (Imag_chu!=-99.)]),30)
plt.figure()
plt.hist(Imag_chu[(Imag_chu<99.) & (Imag_chu!=-99.)],bins=bins,histtype='step',color='blue',linewidth=3,label='With redshift')
plt.hist(Imag_chu_zsp[(Imag_chu_zsp<99.) & (Imag_chu_zsp!=-99.)],bins=bins,histtype='step',color='cyan',linewidth=2,label='Spec-z')
plt.hist(Imag_chu_zph[(Imag_chu_zph<99.) & (Imag_chu_zph!=-99.)],bins=bins,histtype='step',color='red',linewidth=2,label='Photo-z')
plt.hist(Imag_und[(Imag_und<99.) & (Imag_und!=-99.)],bins=30,histtype='step',color='orange',linewidth=3,linestyle='dashed',label='No redshift')
plt.axvline(x=22.5,color='k',label='AGES limit')
plt.xlabel('I-band magnitude')
plt.ylabel('#')
plt.legend()
plt.show()

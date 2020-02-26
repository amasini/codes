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
date = '200226'

# Clean the useless columns from NDWFS and SDWFS
table = Table.read(wd+'nway-master/'+date+'.fits', format='fits')

table.remove_columns(['NDWFS_X_IMAGE', 'NDWFS_Y_IMAGE', 'NDWFS_XPEAK_IMAGE', 'NDWFS_YPEAK_IMAGE', 'NDWFS_FLUX_ISO', 'NDWFS_FLUXERR_ISO', 'NDWFS_MAG_ISO', 'NDWFS_MAGERR_ISO', 'NDWFS_FLUX_APER_01', 'NDWFS_FLUX_APER_02', 'NDWFS_FLUX_APER_03', 'NDWFS_FLUX_APER_04', 'NDWFS_FLUX_APER_05', 'NDWFS_FLUX_APER_06', 'NDWFS_FLUX_APER_07', 'NDWFS_FLUX_APER_08','NDWFS_FLUX_APER_09', 'NDWFS_FLUX_APER_10', 'NDWFS_FLUX_APER_15', 'NDWFS_FLUX_APER_20','NDWFS_FLUXERR_APER_01', 'NDWFS_FLUXERR_APER_02', 'NDWFS_FLUXERR_APER_03', 'NDWFS_FLUXERR_APER_04', 'NDWFS_FLUXERR_APER_05', 'NDWFS_FLUXERR_APER_06', 'NDWFS_FLUXERR_APER_07', 'NDWFS_FLUXERR_APER_08','NDWFS_FLUXERR_APER_09', 'NDWFS_FLUXERR_APER_10', 'NDWFS_FLUXERR_APER_15', 'NDWFS_FLUXERR_APER_20','NDWFS_MAG_APER_01', 'NDWFS_MAG_APER_02', 'NDWFS_MAG_APER_03', 'NDWFS_MAG_APER_04', 'NDWFS_MAG_APER_05', 'NDWFS_MAG_APER_06', 'NDWFS_MAG_APER_07', 'NDWFS_MAG_APER_08','NDWFS_MAG_APER_09', 'NDWFS_MAG_APER_10', 'NDWFS_MAG_APER_15', 'NDWFS_MAG_APER_20','NDWFS_MAGERR_APER_01', 'NDWFS_MAGERR_APER_02', 'NDWFS_MAGERR_APER_03', 'NDWFS_MAGERR_APER_04', 'NDWFS_MAGERR_APER_05', 'NDWFS_MAGERR_APER_06', 'NDWFS_MAGERR_APER_07', 'NDWFS_MAGERR_APER_08','NDWFS_MAGERR_APER_09', 'NDWFS_MAGERR_APER_10', 'NDWFS_MAGERR_APER_15', 'NDWFS_MAGERR_APER_20','NDWFS_KRON_RADIUS', 'NDWFS_BACKGROUND', 'NDWFS_THRESHOLD', 'NDWFS_FLUX_MAX', 'NDWFS_ISOAREA_IMAGE', 'NDWFS_ALPHAPEAK_J2000', 'NDWFS_DELTAPEAK_J2000', 'NDWFS_X2_IMAGE', 'NDWFS_Y2_IMAGE', 'NDWFS_XY_IMAGE', 'NDWFS_CXX_IMAGE', 'NDWFS_CYY_IMAGE', 'NDWFS_CXY_IMAGE', 'NDWFS_CXX_WORLD', 'NDWFS_CYY_WORLD', 'NDWFS_CXY_WORLD', 'NDWFS_A_IMAGE', 'NDWFS_B_IMAGE', 'NDWFS_A_WORLD', 'NDWFS_B_WORLD', 'NDWFS_THETA_IMAGE', 'NDWFS_THETA_WORLD', 'NDWFS_ELONGATION', 'NDWFS_ELLIPTICITY', 'NDWFS_ERRX2_IMAGE', 'NDWFS_ERRY2_IMAGE', 'NDWFS_ERRXY_IMAGE', 'NDWFS_ERRA_IMAGE', 'NDWFS_ERRB_IMAGE', 'NDWFS_ERRTHETA_IMAGE', 'NDWFS_FWHM_IMAGE', 'NDWFS_IMAFLAGS_ISO', 'NDWFS_FLAG_SPLITMATCH','SDWFS_ch1_4', 'SDWFS_ch2_4', 'SDWFS_ch3_4','SDWFS_ch4_4', 'SDWFS_err1_4', 'SDWFS_err2_4', 'SDWFS_err3_4','SDWFS_err4_4','SDWFS_ch1_6', 'SDWFS_ch2_6', 'SDWFS_ch3_6','SDWFS_ch4_6', 'SDWFS_err1_6', 'SDWFS_err2_6', 'SDWFS_err3_6','SDWFS_err4_6'])
table.write(wd+'nway-master/'+date+'-cl.fits', format='fits', overwrite='True')

# Open the cleaned master catalog (should contain only one CDWFS source per line)  
cat=fits.open(wd+'/nway-master/'+date+'-cl.fits')

data=cat[1].data
cols=cat[1].columns
names=cols.names
cat.close()

print('The full CDWFS catalog has',len(data),'sources.')

'''++++++++++++++++++++++ REDSHIFT PART ++++++++++++++++++++++'''
'''+++++++++++++++++++ MATCH with CHUNG+14 +++++++++++++++++++'''
opt_ra=data['NDWFS_RA_J2000']
opt_dec=data['NDWFS_DEC_J2000']
Imag=data['NDWFS_MAG_AUTO']

ir_ra=data['SDWFS_ra']
ir_dec=data['SDWFS_dec']

ra=opt_ra
ra[opt_ra==-99.0]=ir_ra[opt_ra==-99.0]

dec=opt_dec
dec[opt_dec==-99.0]=ir_dec[opt_dec==-99.0]

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

spectro=np.full_like(chu_zsp,True,dtype='bool')
spectro[np.isnan(chu_zsp)]=False

print('Now matching with Chung+14 for redshifts...')

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
				
				z_out_sp.append(cut_zsp[j])
				z_out_ga.append(cut_zph_ga[j])
				z_out_g.append(cut_zph_g[j])
				ebv.append(cut_ebv[j])
				chi_ga.append(cut_chi_ga[j])
				chi_g.append(cut_chi_g[j])
				chi_s.append(cut_chi_s[j])
			else:
				sp.append(False)
				
				z_out_sp.append(-99.0)
				z_out_ga.append(cut_zph_ga[j])
				ebv.append(cut_ebv[j])
				chi_ga.append(cut_chi_ga[j])
				z_out_g.append(cut_zph_g[j])
				chi_g.append(cut_chi_g[j])
				chi_s.append(cut_chi_s[j])
	if found==0:
		Imag_und.append(Imag[i])
		
		z_out_sp.append(-99.0)
		z_out_ga.append(-99.0)
		ebv.append(-99.0)
		chi_ga.append(-99.0)
		z_out_g.append(-99.0)
		chi_g.append(-99.0)
		chi_s.append(-99.0)
print(mat,'matched with Chung in',(time.time()-tin)/60.,'minutes.')

# Write out catalog
table = Table.read(wd+'nway-master/'+date+'-cl.fits', format='fits')

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
table.write(wd+'CDWFS_I-3.6_v'+date+'.fits', format='fits', overwrite='True')
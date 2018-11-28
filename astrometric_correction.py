#following Civano+16, run of wavdetect and match with optical counterparts to correct with reproject_aspect
import numpy as np
import sys
from astropy.io import fits
from astropy.table import Table
import subprocess as s
from ciao_contrib.region.check_fov import FOVFiles
from coords.format import sex2deg

def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*3600*xx)**2 +((pointa[1]-pointb[1])*3600)**2)

wd='/Users/alberto/Desktop/XBOOTES/'

#(dirs,obsid)=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=[0,1],dtype='str')

dirs=['CDWFS','CDWFS']
obsid=['19676','19777','19778','19779']

########### USING METHOD=TRANS FOR THESE ONES
#dirs=['bootes2','CDWFS','CDFWS','murray05','murray05','murray05']
#obsid=['7009','19651','19673','4248','4249','4263']
###########

for i in range(len(obsid)):
	if len(obsid[i]) == 4:
		stem='0'+obsid[i]
	elif len(obsid[i]) == 3:
		stem='00'+obsid[i]
	elif len(obsid[i]) == 5:
		stem=obsid[i]
	ra_pnt=s.check_output('dmkeypar '+wd+'data/'+obsid[i]+'/repro/acisf'+stem+'_repro_evt2.fits RA_PNT echo=yes',shell=True)
	dec_pnt=s.check_output('dmkeypar '+wd+'data/'+obsid[i]+'/repro/acisf'+stem+'_repro_evt2.fits DEC_PNT echo=yes',shell=True)
	pnt=[float(ra_pnt),float(dec_pnt)]
	dat=fits.open(wd+'data/junk_output/'+stem+'_src_new_asol.fits')
	data = dat[1].data
	goodvalues=[]
	for j in range(len(dat[1].data)):
		net=data['NET_COUNTS'][j]
		bkg=data['BKG_COUNTS'][j]
		#significance cut (3.5sigma - hope the formula is correct!)
		if net >= 3.5*np.sqrt(bkg):
			ra_src=data['RA'][j]
			dec_src=data['DEC'][j]
			src=[ra_src,dec_src]
			if distance(src,pnt) < 360.0:
				goodvalues.append(j)
	dat.close()
	newdata=data[goodvalues]
	hdu = fits.BinTableHDU(data=newdata)
	hdu.writeto(wd+'/data/'+obsid[i]+'/src_new_asol.fits',overwrite=True)
	
	'''
	dat2=fits.open(wd+'xbootes_brand+05.fits')
	data2=dat2[1].data
	data2=data2[data2['St']==1]
	goodvalues2=[]
	my_obs = FOVFiles(wd+'data/'+obsid[i]+'/repro/fov_acisI.fits')
	for k in range(len(data2)):
		ra=data2['RAo'][k]
		dec=data2['DEo'][k]
		#ra,dec = sex2deg(str(RAo), str(DEo))
		myobs = my_obs.inside(ra,dec)
		if myobs !=[]:
			goodvalues2.append(k)
		
		newdata2=data2[goodvalues2]
		dat2.close()
	hdu2 = fits.BinTableHDU(data=newdata2)
	hdu2.writeto(wd+'/data/'+obsid[i]+'/acisf'+stem+'_brand+05.fits',overwrite=True)
	
	t=Table(newdata2)
	t['RAo'].name='RA'
	t['DEo'].name='DEC'
	hdu3 = fits.BinTableHDU(data=t)
	hdu3.writeto(wd+'/data/'+obsid[i]+'/acisf'+stem+'_brand+05_renamedcols.fits',overwrite=True)

	s.call('reproject_aspect infile='+wd+'/data/'+obsid[i]+'/src_to_reproject.fits refsrcfile='+wd+'/data/'+obsid[i]+'/acisf'+stem+'_brand+05_renamedcols.fits wcsfile='+wd+'/data/'+obsid[i]+'/repro/acisf'+stem+'_repro_evt2.fits updfile='+wd+'/data/'+obsid[i]+'/repro/*_asol1.fits outfile='+wd+'/data/'+obsid[i]+'/acisf'+stem+'_new_asol.fits radius=2 verbose=1 clobber=yes',shell=True)
	'''
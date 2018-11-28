import numpy as np
import sys
import subprocess as s
import os
from astropy.io import fits

dirs='CDWFS'
obsid='19684'
wd='/Users/alberto/Desktop/XBOOTES/'+dirs+'/'+obsid+'/'

#chandra repro
print('Reprocessing '+obsid+'')
s.call('chandra_repro indir='+wd+' outdir='+wd+'repro/ verbose=0 clobber=yes', shell=True)

if len(obsid) == 4:
    stem='0'+obsid
elif len(obsid) == 3:
    stem='00'+obsid
elif len(obsid) == 5:
    stem=obsid

filename='acisf'+stem+'_repro_'

#creates full band evt file binned to native pixel scale
s.call('dmcopy \"'+wd+'repro/'+filename+'evt2.fits[events][ccd_id=0:3][energy=500:7000]\" \"'+wd+'/repro/'+filename+'05to7keV.fits\" clobber=yes',shell=True)
#creates full band *image* binned 4x4 to native pixel scale
s.call('dmcopy \"'+wd+'repro/'+filename+'05to7keV.fits[events][bin x=::4,y=::4]\" '+wd+'/repro/'+filename+'05to7keV.img opt=image clobber=yes',shell=True)

#creates 9-12 keV image to get background
s.call('dmcopy \"'+wd+'repro/'+filename+'evt2.fits[events][ccd_id=0:3][energy=9000:12000]\" \"'+wd+'/repro/'+filename+'9to12keV.fits\" clobber=yes',shell=True)


#creates exposure map not vignetting-corrected
filename='acisf'+stem+'_repro_9to12keV.fits'

s.call('fluximage \"'+wd+'/repro/'+filename+'[events][ccd_id<4]\" '+wd+'/repro/out2/ unit=time binsize=1 clobber=yes',shell=True)

R=0.1375 #radius of FoV in degrees
dmax=0.0416667 #degrees; 2.5 arcmin
radius=dmax*3600 #150"

ra=s.check_output('dmkeypar '+wd+'/repro/'+filename+' RA_PNT echo=yes',shell=True)
dec=s.check_output('dmkeypar '+wd+'/repro/'+filename+' DEC_PNT echo=yes',shell=True)
exp=s.check_output('dmkeypar '+wd+'/repro/'+filename+' LIVETIME echo=yes',shell=True)
#date=s.check_output('dmkeypar '+wd+'/repro/'+filename+' DATE-OBS echo=yes',shell=True)
#mjd=s.check_output('dmkeypar '+wd+'/repro/'+filename+' MJD_OBS echo=yes',shell=True)
bkg_sur_bri=[]
for i in range(100):
    randra=np.random.uniform(float(ra)+(0.9*R-dmax),float(ra)-(0.9*R-dmax))
    randdec=np.random.uniform(float(dec)+(0.9*R-dmax),float(dec)-(0.9*R-dmax))

    s.call('dmextract \''+wd+'/repro/'+filename+'[bin sky=circle('+str(randra)+'d,'+str(randdec)+'d,'+str(radius)+'\")]\' mode=h outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
    out=s.check_output('dmlist \'counts.fits[col SUR_BRI]\' data,clean | grep -v SUR_',shell=True)
    bkg_sur_bri.append(float(out)/float(exp))
    s.call('rm -f counts.fits',shell=True)
                 
bkg912=np.median(bkg_sur_bri)
e_bkg912=np.std(bkg_sur_bri)
print('='*15)
print('The median 9-12 keV bkg surface brightness for '+obsid+' is')
print(bkg912,e_bkg912)
print('='*15)

bkg057=(0.12+0.27+1.1)*bkg912

dat=fits.open(wd+'/repro/out2/broad_thresh.expmap')
expo=dat[0].data
header=dat[0].header
dat.close()
#rescale expomap
new=expo*bkg057
#write bkgmap
hdu = fits.PrimaryHDU(new,header=header)
hdu.writeto(wd+'/repro/out2/acisf'+stem+'_bkgmap.fits',overwrite=True)

#produce vignetting-corrected expomaps with units=s
#filename='acisf'+stem+'_repro_evt2.fits'
                
s.call('fluximage \"'+wd+'/repro/'+filename+'[events][ccd_id<4]\" '+wd+'/repro/eff_area unit=area binsize=1 clobber=yes',shell=True)
s.call('fluximage \"'+wd+'/repro/'+filename+'[events][ccd_id<4]\" '+wd+'/repro/expo binsize=1 clobber=yes',shell=True)

dat=fits.open(wd+'/repro/eff_area_broad_thresh.expmap')
max_effarea=np.max(dat[0].data)
dat.close()
dat2=fits.open(wd+'/repro/expo_broad_thresh.expmap')
oldexpo=dat2[0].data
oldhead=dat2[0].header
dat2.close()
expo=oldexpo/max_effarea
hdu=fits.PrimaryHDU(expo,header=oldhead)
hdu.writeto(wd+'/repro/out2/acisf'+stem+'_expomap.fits',overwrite=True)

#produce fov file for ACIS-I
filename='acisf'+stem+'_repro_fov1.fits'

fovfits=fits.open(wd+'/repro/'+filename)

fov=fovfits[1].data
fovfits[1].data=fov[fov["CCD_ID"]<4]
fovfits.writeto(wd+'/repro/fov_acisI.fits',overwrite=True)
fovfits.close()

#change header of bkgmap
data_head=fits.getheader(wd+'/repro/acisf'+stem+'_repro_05to7keV.fits', 1)
data_expo=data_head['EXPOSURE']
#take evt file header and wcs
head=fits.getheader(wd+'/repro/out2/acisf'+stem+'_expomap.fits', 0)
head['EXPOSURE']=data_expo
print(head['EXPOSURE'])

### remember to stop the program here if the expomap has acis-s chips, this may bring a shift in the bkg map and in the simulation as well.
#print('Stopping the script before changing header to bkgmap.')
#sys.exit()

#replace old bkgmap header with this one
data, oldheader = fits.getdata(wd+'/repro/out2/acisf'+stem+'_bkgmap.fits', header=True)
hdu = fits.PrimaryHDU(data,header=head)
hdu.writeto(wd+'/repro/out2/acisf'+stem+'_bkgmap.fits',overwrite=True)

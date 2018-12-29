import numpy as np
import sys
import subprocess as s
from astropy.io import fits
from ciao_contrib.region.check_fov import FOVFiles

wd="/Users/alberto/Box Sync/XBOOTES/"

obs=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=1,dtype='str') # List of observations
cts_f=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=2,dtype='float') # Total 0.5-7 keV counts

cutf,cuts,cuth=1.4e-2,1e-2,3.5e-3 # These are the probability cuts in F,S,H bands at 97% rel -> 9240 srcs
# Take the catalog of detected sources (cut at 97% reliability)
cat=fits.open('/Users/alberto/Desktop/prova_cdwfs_merged_cat0.fits')
raf=cat[1].data['RA_F']
decf=cat[1].data['DEC_F']
r90f=cat[1].data['R90_F']
probf=cat[1].data['PROB_F']

ras=cat[1].data['RA_S']
probs=cat[1].data['PROB_S']
rah=cat[1].data['RA_H']
probh=cat[1].data['PROB_H']

probf[raf==0.0]=9999
probs[ras==0.0]=9999
probh[rah==0.0]=9999

ra_src,dec_src,r90_src=[],[],[]
for j in range(len(probf)):
    #if (probf[j] <= cutf or probs[j] <= cuts or probh[j] <= cuth):
    if probf[j] <= cutf: # Select only srcs detected above threshold in F band
        ra_src.append(raf[j])
        dec_src.append(decf[j])
        r90_src.append(r90f[j])


for i in range(len(obs)):
    if len(obs[i]) == 5:
        stem=obs[i]
    elif len(obs[i]) == 4:
        stem='0'+obs[i]
    elif len(obs[i]) == 3:
        stem='00'+obs[i]

    sources_ra,sources_dec=[],[]
    my_obs = FOVFiles(wd+'data/'+obs[i]+'/repro_new_asol/fov_acisI.fits')
    for k in range(len(ra_src)):
        myobs = my_obs.inside(ra_src[k], dec_src[k])
        if myobs !=[]:
            sources_ra.append(ra_src[k])
            sources_dec.append(dec_src[k])
            sources_r90.append(r90_src[k])
    print('In obsid '+obs[i]+' there are '+str(len(sources_ra))+' sources')

    # Create region file with list of sources
    w=open('/Users/alberto/Desktop/try_exclude/acisf'+stem+'_broad_src.reg','w')
    for k in range(len(sources_ra)):
        radius=1.1*sources_r90[k]
        w.write('circle('+str(sources_ra[k])+'d,'+str(sources_dec[k])+'d,'+str(radius)+'\")\n')
    w.close()
      
    s.call('punlearn dmextract')
    s.call('dmextract infile="'+wd+'data/'+obs[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img[bin sky=region(/Users/alberto/Desktop/try_exclude/acisf'+stem+'_broad_src.reg)]" outfile=out2.fits opt=generic mode=h clobber=yes',shell=True)
    
    src_cts=s.check_output('dmlist "out2.fits[cols COUNTS]" data,clean | grep -v COUNTS',shell=True)
    src_area=s.check_output('dmlist "out2.fits[cols AREA]" data,clean | grep -v AREA',shell=True)
    
    x=(cts_f[i]-float(src_cts))*(2048**2/(2048**2-scr_area))
    print(cts_f[i],src_cts,src_area,x)


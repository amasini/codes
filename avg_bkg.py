import numpy as np
import sys
import subprocess as s
import os
import matplotlib.pyplot as plt

def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*xx)**2 +((pointa[1]-pointb[1]))**2)

wd='/Users/alberto/Desktop/XBOOTES/'

R=0.1375 #radius of FoV in degrees
dmax=0.0416667 #degrees; 2.5 arcmin
radius=dmax*3600 #150"
y,yerr=[],[]
#w=open(wd+'/avg_bkg.dat','w')
#w.write('field\tobsid\tsur_bri[cts/pixel2/sec]\terr_sur_bri\n')
k=0
exposure1,exposure2,exposure3,exposure4,exposure5,exposure6,exposure7,exposure8=[],[],[],[],[],[],[],[]
for dirs in os.listdir(wd):
    if os.path.isdir(dirs) == True:
        print(dirs)
        for obsid in os.listdir(wd+dirs+'/'):
            if os.path.isdir(wd+dirs+'/'+obsid+'/') == True:
                if obsid != 'AAA':
                    k=k+1
                    print(k,dirs,obsid)
                    if len(obsid) == 4:
                        stem='0'+obsid
                    elif len(obsid) == 3:
                        stem='00'+obsid
                    elif len(obsid) == 5:
                        stem=obsid
                    else:
                        print('Something\'s wrong with '+dirs+'/'+obsid+'/')
               	        sys.exit()
                    filename='acisf'+stem+'_repro_9to12keV.fits'
                
                    #ra=s.check_output('dmkeypar '+wd+dirs+'/'+obsid+'/repro/'+filename+' RA_PNT echo=yes',shell=True)
                    #dec=s.check_output('dmkeypar '+wd+dirs+'/'+obsid+'/repro/'+filename+' DEC_PNT echo=yes',shell=True)
		    #exp=s.check_output('dmkeypar '+wd+dirs+'/'+obsid+'/repro/'+filename+' LIVETIME echo=yes',shell=True)
                    #date=s.check_output('dmkeypar '+wd+dirs+'/'+obsid+'/repro/'+filename+' DATE-OBS echo=yes',shell=True)
                    mjd=s.check_output('dmkeypar '+wd+dirs+'/'+obsid+'/repro/'+filename+' MJD_OBS echo=yes',shell=True)
                    var=mjd
                    if dirs=='bootes2':
                        exposure1.append(float(var))
                    elif dirs=='bootes_cycle17':
                        exposure2.append(float(var))
                    elif dirs=='bootes3':
                        exposure3.append(float(var))
                    elif dirs=='CDWFS':
                        exposure4.append(float(var))
                    elif dirs=='bootes_mod_depth_cycle16':
                        exposure5.append(float(var))
                    elif dirs=='miscellaneous':
                        exposure6.append(float(var))
                    elif dirs=='murray05':
                        exposure7.append(float(var))
                    elif dirs=='xbootes_IRAGN':
                        exposure8.append(float(var))
                    
                    #print(date)
                    '''
                    bkg_sur_bri=[]
                    for i in range(100):
                        randra=np.random.uniform(float(ra)+(0.9*R-dmax),float(ra)-(0.9*R-dmax))
                        randdec=np.random.uniform(float(dec)+(0.9*R-dmax),float(dec)-(0.9*R-dmax))

                        s.call('dmextract \''+wd+dirs+'/'+obsid+'/repro/'+filename+'[bin sky=circle('+str(randra)+'d,'+str(randdec)+'d,'+str(radius)+'\")]\' mode=h outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
                        out=s.check_output('dmlist \'counts.fits[col SUR_BRI]\' data,clean | grep -v SUR_',shell=True)
                        bkg_sur_bri.append(float(out)/float(exp))
                        #s.call('rm -f counts.fits',shell=True)
                 
                    y.append(np.median(bkg_sur_bri))
                    yerr.append(np.std(bkg_sur_bri))
                    print(y[-1],yerr[-1])
                    #w.write(str(dirs)+'\t'+str(obsid)+'\t'+str(np.median(bkg_sur_bri))+'\t'+str(np.std(bkg_sur_bri))+'\n')
                    
                    plt.figure()
                    plt.hist(bkg_sur_bri,histtype='step',bins=10)
                    plt.show()
                    sys.exit()
                    '''
#sys.exit()  
#w.close()
(y,yerr)=np.genfromtxt(wd+'/avg_bkg.dat',unpack=True,usecols=[2,3],skip_header=1)
name=np.genfromtxt(wd+'/avg_bkg.dat',unpack=True,usecols=[0],skip_header=1,dtype='str')
bootes2=y[name=='bootes2']
ebootes2=yerr[name=='bootes2']
bootescy17=y[name=='bootes_cycle17']
ebootescy17=yerr[name=='bootes_cycle17']
bootes3=y[name=='bootes3']
ebootes3=yerr[name=='bootes3']
cdwfs=y[name=='CDWFS']
ecdwfs=yerr[name=='CDWFS']
bootes_moddepth=y[name=='bootes_mod_depth_cycle16']
ebootes_moddepth=yerr[name=='bootes_mod_depth_cycle16']
misc=y[name=='miscellaneous']
emisc=yerr[name=='miscellaneous']
murray05=y[name=='murray05']
emurray05=yerr[name=='murray05']
iragn=y[name=='xbootes_IRAGN']
eiragn=yerr[name=='xbootes_IRAGN']

x1=np.arange(1,len(bootes2)+1,1)
x2=np.arange(x1[-1]+1,x1[-1]+1+len(bootescy17),1)
x3=np.arange(x2[-1]+1,x2[-1]+len(bootes3)+1,1)
x4=np.arange(x3[-1]+1,x3[-1]+len(cdwfs)+1,1)
x5=np.arange(x4[-1]+1,x4[-1]+len(bootes_moddepth)+1,1)
x6=np.arange(x5[-1]+1,x5[-1]+len(misc)+1,1)
x7=np.arange(x6[-1]+1,x6[-1]+len(murray05)+1,1)
x8=np.arange(x7[-1]+1,x7[-1]+len(iragn)+1,1)

plt.figure()
plt.errorbar(exposure1,bootes2,yerr=ebootes2,fmt='o',color='red',label='Bootes2')
plt.errorbar(exposure2,bootescy17,yerr=ebootescy17,fmt='o',color='g',label='Bootes_cy17')
plt.errorbar(exposure3,bootes3,yerr=ebootes3,fmt='o',color='b',label='Bootes3')
plt.errorbar(exposure4,cdwfs,yerr=ecdwfs,fmt='o',color='k',label='CDWFS')
plt.errorbar(exposure5,bootes_moddepth,yerr=ebootes_moddepth,fmt='o',color='pink',label='Bootes_mod_depth')
plt.errorbar(exposure6,misc,yerr=emisc,fmt='o',color='y',label='Misc')
plt.errorbar(exposure7,murray05,yerr=emurray05,fmt='o',color='brown',label='Murray05')
plt.errorbar(exposure8,iragn,yerr=eiragn,fmt='o',color='orange',label='XB_IRAGN')
plt.xlabel('MJD')
plt.ylabel(r'Avg Bkg surf brightness [cts pixel$^{-2}$ s$^{-1}$]')
#plt.xscale('log')
plt.yscale('log')
#plt.axis([4000,200000,1e-7,6e-7])
plt.legend(ncol=2)
plt.tight_layout()
plt.savefig(wd+'avg_bkg_MJD.png',format='png')

sys.exit()

plt.figure()
plt.errorbar(x1,bootes2,yerr=ebootes2,fmt='o',color='red',label='Bootes2')
plt.errorbar(x2,bootescy17,yerr=ebootescy17,fmt='o',color='g',label='Bootes_cy17')
plt.errorbar(x3,bootes3,yerr=ebootes3,fmt='o',color='b',label='Bootes3')
plt.errorbar(x4,cdwfs,yerr=ecdwfs,fmt='o',color='k',label='CDWFS')
plt.errorbar(x5,bootes_moddepth,yerr=ebootes_moddepth,fmt='o',color='pink',label='Bootes_mod_depth')
plt.errorbar(x6,misc,yerr=emisc,fmt='o',color='y',label='Misc')
plt.errorbar(x7,murray05,yerr=emurray05,fmt='o',color='brown',label='Murray05')
plt.errorbar(x8,iragn,yerr=eiragn,fmt='o',color='orange',label='XB_IRAGN')
plt.xlabel('Bootes OBSID')
plt.ylabel(r'Avg Bkg surf brightness [cts pixel$^{-2}$ s$^{-1}$]')
plt.yscale('log')
plt.legend(ncol=2)
plt.tight_layout()
plt.savefig(wd+'avg_bkg_OBSID.png',format='png')

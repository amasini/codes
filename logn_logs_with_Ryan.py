import numpy as np
import matplotlib.pyplot as plt
import sys
import subprocess as s
import os
import scipy.special
from scipy.integrate import quad
def integrand(x,a):
    return x**(-a+1)

#from ciao_contrib.region.check_fov import FOVFiles

#Luo+17 CDFS dn/ds for AGN
def dnds(s,band):
    n=[]
    if band=='s':
        k,fb,b1,b2=161.96,7.1e-15,1.52,2.45
    elif band=='h':
        k,fb,b1,b2=453.70,8.9e-15,1.46,2.72
    else:
        print('Band not recognised. Options are s and h.')
    for i in range(len(s)):
        if s[i] <= fb:
            n.append(k*1e14*(s[i]/1.e-14)**(-b1))
        else:
            n.append(k*1e14*(fb/1.e-14)**(b2-b1)*(s[i]/1.e-14)**(-b2))
    return n

wd='/Users/alberto/Desktop/XBOOTES/'

flux=np.logspace(np.log10(5e-17),np.log10(1e-12),101) #this is the full band flux
centers0=list((flux[i+1]+flux[i])/2 for i in range(0,len(flux)-1))

#photon index as a function of flux, from Hickox & Markevitch 2006
conv_factor_soft,conv_factor_hard=[],[]
gamma_s=0.59+1.7*scipy.special.erf(np.log10((np.array(flux)/1e-15)+4.6/2.0))
gamma_h=-16.+18.*scipy.special.erf(np.log10((np.array(flux)/1e-15)+12./0.41))
for k in range(len(gamma_s)):
    I = quad(integrand, 0.5, 2, args=(gamma_s[k])) #flux in soft band
    I1 = quad(integrand, 0.5, 7, args=(gamma_s[k])) #flux in full band
    conv_factor_soft.append(I[0]/I1[0])

    I = quad(integrand, 2, 7, args=(gamma_h[k])) #flux in soft band
    I1 = quad(integrand, 0.5, 7, args=(gamma_h[k])) #flux in full band
    conv_factor_hard.append(I[0]/I1[0])

conv_factor_soft=np.array(conv_factor_soft)
conv_factor_hard=np.array(conv_factor_hard)
#f_soft=0.4516*flux #Gamma=1.8
#f_soft=0.3283*flux #Gamma=1.4
f_soft=conv_factor_soft*flux
centers1=list((f_soft[i+1]+f_soft[i])/2 for i in range(0,len(f_soft)-1))
#f_hard=0.5484*flux #Gamma=1.8
#f_hard=0.6716*flux #Gamma=1.4
f_hard=conv_factor_hard*flux
centers2=list((f_hard[i+1]+f_hard[i])/2 for i in range(0,len(f_hard)-1))

#AGNs per unit flux per squared degree, in soft and hard bands
dnds_s=dnds(centers1,'s')
dnds_h=dnds(centers2,'h')
	
#Write file with sources
w=open(wd+'poiss_rand_new.reg','w')
w2=open(wd+'poiss_rand_new.dat','w')
w2.write('Full flux\tRA\tDEC\n')
#choose rectangular area of 4x3.5 deg2 centered on the center of the field
(minra,maxra)=(215.82,220.1)
(minde,maxde)=(32.2,36.2)
area=((maxra-minra)/57.29*(np.sin(maxde/57.29)-np.sin(minde/57.29)))*57.29**2
#print(area)
dn_s,dn_h=[],[]
try0,try1=[],[]
for i in range(len(dnds_s)):
	
    #sources per square degree in each flux bin
    dn_s.append(dnds_s[i]*(f_soft[i+1]-f_soft[i]))
    dn_h.append(dnds_h[i]*(f_hard[i+1]-f_hard[i]))

    #sources in total in each flux bin
    n_s=dnds_s[i]*(f_soft[i+1]-f_soft[i])*area
    n_h=dnds_h[i]*(f_hard[i+1]-f_hard[i])*area
    try0.append(n_s)

    #Poissonian realization of the total number of sources
    N_s=np.random.poisson(n_s)
    N_h=np.random.poisson(n_h)
    try1.append(N_s)

    j=0
    while j < N_s:
        randdec=np.random.uniform(minde,maxde)
        prob=np.random.uniform(0,1)
        if np.cos(randdec*np.pi/180.) > prob:
            j=j+1
            randra=np.random.uniform(minra,maxra)
            w.write('circle('+str(randra)+'d,'+str(randdec)+'d,1\")\n')
            w2.write(str(flux[i])+'\t'+str(randra)+'\t'+str(randdec)+'\n')
    '''
    k=0
    while k < N_h:
        randdec=np.random.uniform(minde,maxde)
        prob=np.random.uniform(0,1)
        if np.cos(randdec*np.pi/180.) > prob:
            k=k+1
            randra=np.random.uniform(minra,maxra)
            w.write('circle('+str(randra)+'d,'+str(randdec)+'d,1\")\n')
            w2.write(str(centers0[i])+'\t'+str(randra)+'\t'+str(randdec)+'\n')
	'''
w.close()
w2.close()



##########################################################################
#recover logn-logs of input sources written in poiss_rand_new.dat
full_flux=np.genfromtxt(wd+'poiss_rand_new.dat',skip_header=1,usecols=0)
bins00=np.logspace(np.log10(5e-17),np.log10(1e-12),101)
centers00=list((bins00[i+1]+bins00[i])/2 for i in range(0,len(bins00)-1))

conv_factor_soft2,conv_factor_hard2=[],[]
gamma_s2=0.59+1.7*scipy.special.erf(np.log10((np.array(full_flux)/1e-15)+4.6/2.0))
gamma_h2=-16.+18.*scipy.special.erf(np.log10((np.array(full_flux)/1e-15)+12./0.41))
for k in range(len(gamma_s2)):
    I = quad(integrand, 0.5, 2, args=(gamma_s2[k])) #flux in soft band
    I1 = quad(integrand, 0.5, 7, args=(gamma_s2[k])) #flux in full band
    conv_factor_soft2.append(I[0]/I1[0])
    
    I = quad(integrand, 2, 7, args=(gamma_h2[k])) #flux in soft band
    I1 = quad(integrand, 0.5, 7, args=(gamma_h2[k])) #flux in full band
    conv_factor_hard2.append(I[0]/I1[0])

conv_factor_soft2=np.array(conv_factor_soft2)
conv_factor_hard2=np.array(conv_factor_hard2)

#f_soft_in=0.3283*full_flux #Gamma=1.4
#bins11=0.3283*bins00
#centers11=0.3283*np.array(centers00)

f_soft_in=conv_factor_soft2*full_flux
bins11=conv_factor_soft*bins00
centers11=list((bins11[i+1]+bins11[i])/2 for i in range(0,len(bins11)-1))
#centers11=conv_factor_soft*np.array(centers00)

#f_hard_in=0.6716*full_flux #Gamma=1.4
#bins22=0.6716*bins00
#centers22=0.6716*np.array(centers00)

f_hard_in=conv_factor_hard2*full_flux
bins22=conv_factor_hard*bins00
#centers22=list((f_hard_in[i+1]+f_hard_in[i])/2 for i in range(0,len(f_hard_in)-1))
#centers22=(1-conv_factor_soft)*np.array(centers00)
centers22=list((bins22[i+1]+bins22[i])/2 for i in range(0,len(bins22)-1))

quante_soft,bincenters_s=np.histogram(f_soft_in,bins=bins11)
quante_hard,bincenters_h=np.histogram(f_hard_in,bins=bins22)

quante_soft_perarea=quante_soft/area
quante_hard_perarea=quante_hard/area

ncum_soft_in=list(reversed(np.cumsum(list(reversed(quante_soft_perarea)))))
ncum_hard_in=list(reversed(np.cumsum(list(reversed(quante_hard_perarea)))))
############################################################################

#check the cumulatives
ncum_soft=list(reversed(np.cumsum(list(reversed(dn_s)))))
ncum_hard=list(reversed(np.cumsum(list(reversed(dn_h)))))

ratio=np.array(quante_soft_perarea)/np.array(dn_s)
plt.figure()
plt.plot(centers1,ratio,'k-')
plt.xscale('log')
plt.show()

ncum0=list(reversed(np.cumsum(list(reversed(try0)))))
ncum1=list(reversed(np.cumsum(list(reversed(try1)))))


#make the plots
f,(ax1,ax2)=plt.subplots(1,2)
ax1.plot(centers1,ncum_soft,'r-',linewidth=2,label='Luo+17 0.5-2 keV')
#ax1.plot(centers1,ncum0,'r-',linewidth=2,label='Luo+17 0.5-2 keV')
#ax1.plot(centers1,ncum1,'r--',linewidth=2,label='Luo+17 0.5-2 keV')
ax1.plot(centers11,ncum_soft_in,'r--',linewidth=2, label='In input file')
ax1.set_xlabel('S [cgs]')
ax1.set_ylabel(r'N(>S) [deg$^-2$]')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.axis([2e-18,1e-13,5,2.1e5])
ax1.legend()

ax2.plot(centers2,ncum_hard,'b-',linewidth=2,label='Luo+17 2-7 keV')
ax2.plot(centers22,ncum_hard_in,'b--',linewidth=2, label='In input file')
ax2.set_xlabel('S [cgs]')
ax2.set_ylabel(r'N(>S) [deg$^-2$]')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.axis([8e-18,1e-13,5,8e4])
ax2.legend()
plt.tight_layout()
plt.show()
sys.exit()

''''
print(np.cumsum(N))
sys.exit()

#PIMMS Gamma=1.8, NH=1.2e20
f_soft=0.4516*flux
f_hard=0.5484*flux 
#PIMMS Gamma=1.4, NH=1.2e20
f_soft2=0.3283*flux
f_hard2=0.6716*flux 

centers=list((f_soft[i+1]+f_soft[i])/2 for i in range(0,len(f_soft)-1))
centers2=list((f_hard[i+1]+f_hard[i])/2 for i in range(0,len(f_hard)-1))

centers3=list((f_soft2[i+1]+f_soft2[i])/2 for i in range(0,len(f_soft2)-1))
centers4=list((f_hard2[i+1]+f_hard2[i])/2 for i in range(0,len(f_hard2)-1))

y=dnds(centers,'s') #sources per flux unit per square degree
y2=dnds(centers2,'h')

y3=dnds(centers3,'s')
y4=dnds(centers4,'h')

dn_s,dn_h=[],[]
dn_s2,dn_h2=[],[]
for i in range(len(y)):
    dn_s.append((f_soft[i+1]-f_soft[i])*y[i]) #sources per square degree
    dn_h.append((f_hard[i+1]-f_hard[i])*y2[i])
    
    dn_s2.append((f_soft2[i+1]-f_soft2[i])*y3[i])
    dn_h2.append((f_hard2[i+1]-f_hard2[i])*y4[i])

'''

'''  
def dnds(s,band):
    if band=='s':
        k,fb,b1,b2=161.96,7.1e-15,1.52,2.45
    elif band=='h':
        k,fb,b1,b2=453.70,8.9e-15,1.46,2.72
    else:
        print('Band not recognised. Options are s and h.')
    if s <= fb:
        n=k*1e14*(s/1.e-14)**(-b1)
    else:
        n=k*1e14*(fb/1.e-14)**(b2-b1)*(s/1.e-14)**(-b2)
    return n
'''
'''

#choose rectangular area of 3x3.6 deg2 centered on the center of the field
#(minra,maxra)=(216.14,219.86)
#(minde,maxde)=(32.44,36.09)
#choose rectangular area of 4x3.5 deg2 centered on the center of the field
(minra,maxra)=(215.82,220.1)
(minde,maxde)=(32.2,36.2)
area=((maxra-minra)/57.29*(np.sin(maxde/57.29)-np.sin(minde/57.29)))*57.29**2
dnds_soft=area*np.array(dn_s) #effective number of sources in each flux bin in the area
dnds_hard=area*np.array(dn_h)
#poisson realization of lognlogs
poiss_soft=np.random.poisson(dnds_soft)
poiss_hard=np.random.poisson(dnds_hard)
w=open(wd+'poiss_rand_new.reg','w')
w2=open(wd+'poiss_rand_new.dat','w')
w2.write('Full flux\tRA\tDEC\n')

#CHECK THIS, MAY BE WRONG -> FULL BAND FLUX[i] (centers0[i]) == to soft_flux[i] and hard_flux[i]?
#probably is better to interpolate, starting from the two separate logn logs creating one "full" logn logs
for i in range(len(centers0)):
    j=0
    while j < int(poiss_soft[i]):
        randdec=np.random.uniform(minde,maxde)
        prob=np.random.uniform(0,1)
        if np.cos(randdec*np.pi/180.) > prob:
            j=j+1
            randra=np.random.uniform(minra,maxra)
            w.write('circle('+str(randra)+'d,'+str(randdec)+'d,1\")\n')
            w2.write(str(centers0[i])+'\t'+str(randra)+'\t'+str(randdec)+'\n')
    k=0
    while k < int(poiss_hard[i]):
        randdec=np.random.uniform(minde,maxde)
        prob=np.random.uniform(0,1)
        if np.cos(randdec*np.pi/180.) > prob:
            k=k+1
            randra=np.random.uniform(minra,maxra)
            w.write('circle('+str(randra)+'d,'+str(randdec)+'d,1\")\n')
            w2.write(str(centers0[i])+'\t'+str(randra)+'\t'+str(randdec)+'\n')
w.close()
w2.close()

sys.exit()

#multiply for total area in deg^2
a=(16.9/60)**2
dnds_soft=a*np.array(dn_s)
dnds_hard=a*np.array(dn_h)

#aa,bb=[],[]
# LOOP on all the OBSIDs
for dirs in os.listdir(wd):
    if os.path.isdir(dirs) == True:
        print(dirs)
        for obsid in os.listdir(wd+dirs+'/'):
            if os.path.isdir(wd+dirs+'/'+obsid+'/') == True:
                if dirs == 'murray05':
                    print(dirs,obsid)
                    if len(obsid) == 4:
                        stem='0'+obsid
                    elif len(obsid) == 3:
                        stem='00'+obsid
                    elif len(obsid) == 5:
                        stem=obsid
                    else:
                        print('Something\'s wrong with '+dirs+'/'+obsid+'/')
               	        sys.exit()

                    filename='acisf'+stem+'_repro_fov1.fits'
                    poiss_soft=np.random.poisson(dnds_soft)
                    poiss_hard=np.random.poisson(dnds_hard)

                    print('In obsID '+obsid+', there are '+str(np.sum(poiss_soft)+np.sum(poiss_hard))+' sources; '+str(np.sum(poiss_soft))+' in the S band, '+str(np.sum(poiss_hard))+' in the H band.')

                    ra=s.check_output('dmkeypar '+wd+dirs+'/'+obsid+'/repro/'+filename+' RA_PNT echo=yes',shell=True)
                    dec=s.check_output('dmkeypar '+wd+dirs+'/'+obsid+'/repro/'+filename+' DEC_PNT echo=yes',shell=True)
                    #aa.append(float(ra))
                    #bb.append(float(dec))
                    
                    R=16.9/60.
                    w=open(wd+dirs+'/'+obsid+'/repro/poiss_rand.reg','w')
                    w2=open(wd+dirs+'/'+obsid+'/repro/poiss_rand.dat','w')
                    w2.write('Full flux\tRA\tDEC\n')
                    for i in range(len(centers0)):
                        for j in range(int(poiss_soft[i])):
                            ii=[]
                            while ii==[]:
                                randra=np.random.uniform(float(ra)-(R/1),float(ra)+(R/1))
                                randdec=np.random.uniform(float(dec)-(R/1),float(dec)+(R/1))
                                my_obs = FOVFiles(wd+dirs+'/'+obsid+'/repro/fov_acisI.fits')
                                ii = my_obs.inside(randra, randdec)
                                if ii!=[]:
                                    w.write('circle('+str(randra)+'d,'+str(randdec)+'d,1\")\n')
                                    w2.write(str(centers0[i])+'\t'+str(randra)+'\t'+str(randdec)+'\n')

                        for k in range(int(poiss_hard[i])):
                            ii=[]
                            while ii==[]:
                                randra=np.random.uniform(float(ra)-(R/1),float(ra)+(R/1))
                                randdec=np.random.uniform(float(dec)-(R/1),float(dec)+(R/1))
                                my_obs = FOVFiles(wd+dirs+'/'+obsid+'/repro/fov_acisI.fits')
                                ii = my_obs.inside(randra, randdec)
                                if ii!=[]:
                                    w.write('circle('+str(randra)+'d,'+str(randdec)+'d,1\")\n')
                                    w2.write(str(centers0[i])+'\t'+str(randra)+'\t'+str(randdec)+'\n')
                    w.close()
                    w2.close()
                    
#print(min(aa),min(bb))
#print(max(aa),max(bb))
'''


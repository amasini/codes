# THIS SCRIPT ASSUMES A LOGN-LOGS AND CREATES A LIST OF SOURCES TO BE INPUT IN A SIMULATION
import numpy as np
import matplotlib.pyplot as plt
import sys
import subprocess as s
import os
import scipy.special
from scipy.integrate import quad
def integrand(x,a):
    return x**(-a+1)

def gauss(x,mu,sigma):
	g=np.exp(-(x-mu)**2/(2*sigma**2))
	return g

#Lehmer+12 CDFS dn/ds for AGN
def dnds(fx,band):
	if band == 'broad':
		k,b1,b2,fb,fref = 562.2e14,-1.34,-2.35,8.1e-15,1e-14
	elif band == 'soft':
		k,b1,b2,fb,fref = 169.56e14,-1.49,-2.48,6.0e-15,1e-14
	elif band == 'hard':
		k,b1,b2,fb,fref = 573.13e14,-1.32,-2.55,6.4e-15,1e-14
	else:
		print('Band not recognized. Exit')
		sys.exit()
		
	k1 = k*(fb/fref)**(b1-b2)
	
	if type(fx) == float:
		if fx <= fb:
			return k*(fx/fref)**b1
		else:
			return k1*(fx/fref)**b2
	elif (type(fx) == list)	or (type(fx) == np.ndarray):
		if type(fx) == list:
			fx = np.array(fx)
		aux = k*(fx/fref)**b1
		aux[fx > fb] = k1*(fx[fx > fb]/fref)**b2
		return aux
'''		
def dnds(s):
    n=[]
    k,fb,b1,b2=562.20,8.1e-15,1.34,2.35 #full band Lehmer
    #k,fb,b1,b2=674.64,8.1e-15,1.34,2.35 #full band Lehmer x1.2
    #k,fb,b1,b2=169.56,6.0e-15,1.49,2.48 #soft band Lehmer
    #k,fb,b1,b2=573.13,6.4e-15,1.32,2.55 #hard band Lehmer
    
    for i in range(len(s)):
        if s[i] <= fb:
            n.append(k*1e14*(s[i]/1.e-14)**(-b1))
        else:
            n.append(k*1e14*(fb/1.e-14)**(b2-b1)*(s[i]/1.e-14)**(-b2))
    return n
'''
def dnds_k(s): # logn-logs by Kenter+05 for Xbootes
    n=[]
    k,fb,b1,b2=674.64,8.1e-15,1.34,2.35 #full band Lehmer x1.2
    #k,fb,b1,b2=129,8.2e-15,1.74,2.60 #0.5-2.0 keV band Kenter
    for i in range(len(s)):
        if s[i] <= fb:
            n.append(k*1e14*(s[i]/1.e-14)**(-b1))
        else:
            n.append(k*1e14*(fb/1.e-14)**(b2-b1)*(s[i]/1.e-14)**(-b2))
    return n
    
def dnds_b(s): # logn-logs by Brandt+05 for CDFN
    n=[]
    k,b1=3970,0.67 #0.5-2.0 keV band Brandt (already cumulated)
    for i in range(len(s)):
        n.append(k*(s[i]/1.e-16)**(-b1))

    return n

def dnds_g(s): # logn-logs by Georgakakis+08 for combined fields
    n=[]
    k,fb,b1,b2=3.74,2.63e-14,1.58,2.48 #full band 0.5-10 keV
    for i in range(len(s)):
        if s[i] <= fb:
            n.append(k*1e16*(s[i]/1.e-14)**(-b1))
        else:
            n.append(k*1e16*(fb/1.e-14)**(b2-b1)*(s[i]/1.e-14)**(-b2))
    return n

wd='/Users/alberto/Desktop/XBOOTES/'
'''
n2=4
N2=[]
for i in range(1000):
	N2.append(np.random.poisson(n2))
plt.figure()
plt.hist(N2,bins=30)
plt.show()
sys.exit()
'''

# Define the Gamma PDF
xvals = np.linspace(1.3, 2.3, 101)
mu=1.8
sigma=0.2
p=[]
for ii in range(len(xvals)):
	p.append(gauss(xvals[ii],mu,sigma))
p=np.array(p)
new=p/np.sum(p)
#plt.figure()
#plt.plot(xvals,new,'k-')
#plt.axvline(mu)
#plt.axvline(mu-sigma)
#plt.axvline(mu+sigma)
#plt.show()

flux=np.logspace(np.log10(5e-17),np.log10(1e-12),101) #this is the full band flux
centers0=list((flux[i+1]+flux[i])/2 for i in range(0,len(flux)-1))

#AGNs per unit flux per squared degree, in soft and hard bands
dnds=dnds(centers0,'hard') #Lehmer
dnds_k=dnds_k(centers0) #Kenter
dnds_g=dnds_g(centers0) #Georgakakis
	
#Write file with sources
w=open(wd+'poiss_rand_hard.reg','w')
w2=open(wd+'poiss_rand_hard.dat','w')
w2.write('Hard flux \t RA \t DEC \n')
#choose rectangular area of 4x3.5 deg2 centered on the center of the field
(minra,maxra)=(215.82,220.1)
(minde,maxde)=(32.2,36.2)
area=((maxra-minra)/57.29*(np.sin(maxde/57.29)-np.sin(minde/57.29)))*57.29**2
#print(area)
dn,dn_k,dn_g=[],[],[]
for i in range(len(dnds)):
	
	#sources per square degree in each flux bin
	dn.append(dnds[i]*(flux[i+1]-flux[i]))
	dn_k.append(dnds_k[i]*(flux[i+1]-flux[i]))
	dn_g.append(dnds_g[i]*(flux[i+1]-flux[i]))

	#sources in total in each flux bin
	n=dnds[i]*(flux[i+1]-flux[i])*area

	#Poissonian realization of the total number of sources
	N=np.random.poisson(n)

	j=0
	while j < N:
		randdec=np.random.uniform(minde,maxde)
		prob=np.random.uniform(0,1)
		if np.cos(randdec*np.pi/180.) > prob:
			j=j+1
			randra=np.random.uniform(minra,maxra)
			w.write('circle('+str(randra)+'d,'+str(randdec)+'d,1\")\n')
			w2.write(str(flux[i])+'\t'+str(randra)+'\t'+str(randdec)+'\n')
			#random_gamma=np.random.choice(xvals, p=new)
			#w2.write(str(flux[i])+' \t '+str(randra)+' \t '+str(randdec)+' \t '+str(round(random_gamma,2))+'\n')
w.close()
w2.close()

sys.exit()

##########################################################################
#recover logn-logs of input sources written in poiss_rand.dat
full_flux=np.genfromtxt(wd+'poiss_rand_lehmer_broad.dat',skip_header=1,usecols=0)
bins00=np.logspace(np.log10(5e-17),np.log10(1e-12),101)
centers00=list((bins00[i+1]+bins00[i])/2 for i in range(0,len(bins00)-1))

quante,bincenters_s=np.histogram(full_flux,bins=bins00)

quante_perarea=quante/area # Here, area is larger than 9.3 deg2!

ncum_in=list(reversed(np.cumsum(list(reversed(quante_perarea)))))

#recover logn-logs of input sources written in poiss_rand_lehmer_filtered.dat
full_flux=np.genfromtxt(wd+'poiss_rand_lehmer_filtered.dat',skip_header=1,usecols=0)
bins00=np.logspace(np.log10(5e-17),np.log10(1e-12),101)
centers00=list((bins00[i+1]+bins00[i])/2 for i in range(0,len(bins00)-1))

quante,bincenters_s=np.histogram(full_flux,bins=bins00)
area=9.3
quante_perarea=quante/area

ncum_in_2=list(reversed(np.cumsum(list(reversed(quante_perarea)))))
############################################################################

#check the cumulatives
ncum=list(reversed(np.cumsum(list(reversed(dn)))))
ncum_k=list(reversed(np.cumsum(list(reversed(dn_k)))))
ncum_b=dnds_b(centers0)
ncum_g=list(reversed(np.cumsum(list(reversed(dn_g)))))

#take civano logn-logs from file
#(civano_f,civano_ns)=np.genfromtxt(wd+'lognlogs_civano_cosmoslegacy.dat',unpack=True,skip_header=1)
#cf=1.019 #from 2-10 to 0.5-7, Gamma=1.4
#cf=1.351 #from 2-10 to 0.5-7, Gamma=1.8
#civano_f=civano_f*cf 

#take soft lehmer logn-logs from file
#(lehmer_f,lehmer_ns)=np.genfromtxt(wd+'logn-logs_lehmer_full.dat',unpack=True)

#make the plots
f,ax1=plt.subplots(1,1)
ax1.plot(centers0,ncum,'r-',linewidth=2,label='Lehmer+12 0.5-7 keV')
ax1.plot(centers0,ncum_k,'g-',linewidth=2,label='Lehmer+12 0.5-7 X1.2 keV')
#ax1.plot(lehmer_f,lehmer_ns,'c-',linewidth=2,label='Lehmer+12 0.5-2 keV')
#ax1.plot(centers0,ncum_g,'c-',linewidth=2,label='Georgakakis+08 0.5-10 keV')
#ax1.plot(civano_f,civano_ns,'b*',ms=15,label=r'Civano+06 2-10 keV, $\Gamma=1.8$')
#ax1.plot(centers0,ncum_b,'b-',linewidth=2,label='Brandt+01 0.5-2 keV')
ax1.plot(centers00,ncum_in,'k--',linewidth=2, label='In input file')
ax1.plot(centers00,ncum_in_2,'b--',linewidth=2, label='In input file, filtered')
ax1.set_xlabel('S [cgs]')
ax1.set_ylabel(r'N(>S) [deg$^-2$]')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.axis([5e-17,1e-13,5,1e5])
ax1.legend()
plt.tight_layout()
plt.show()


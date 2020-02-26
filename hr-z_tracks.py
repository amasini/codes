# Try to compute HR-z tracks with more physically motivated models (PL, MYTORUS, CT Refl Dom)
import numpy as np
import sys
import matplotlib.pyplot as plt
import subprocess as s


def write_commands(model,nh,redshift,gamma):
	w=open(wd+'hr-track.xcm','w')
	if model == 'pl':
	
		w.write('mo pha*zpha*(zpow+const*zpow)\n')
		w.write('0.0104 -1\n')	# 1 Nh gal
		w.write(str(nh/1.e22)+' -1\n') # 2 NH
		w.write(str(redshift)+' -1\n') # 3 z
		w.write(str(gamma)+' -1\n') # 4 Gamma
		w.write('=p3\n') # 5 = 3
		w.write('1 -1\n') # 6 Norm
		w.write('0.02 -1\n') #7 scattered fraction - 2%
		w.write('=p4\n') # 8 Gamma
		w.write('=p5\n') # 9 Redshift
		w.write('=p6\n') # 10 Norm
	
	elif model == 'myt':
	
		w.write('mo pha*(zpow*etable{../../../mytorus/mytorus_Ezero_v00.fits}+const*atable{../../../mytorus/mytorus_scatteredH500_v00.fits} + const*atable{../../../mytorus/mytl_V000010nEp000H500_v00.fits}+const*zpow)\n')
		w.write('0.0104 -1\n')	# 1 Nh gal
		w.write(str(gamma)+' -1\n') # 2 Gamma
		w.write(str(redshift)+' -1\n') # 3 z
		w.write('1 -1\n') # 4 Norm
		w.write(str(nh/1.e24)+' -1\n') # 5 NH
		w.write('90 -1\n') # 6 IncAng
		w.write('=p3\n') # 7 = 3
		w.write('1. -1\n') # 8 A_s
		w.write('=p5\n') # 9 = 5
		w.write('=p6\n') # 10 = 6
		w.write('=p2\n') # 11 = 2
		w.write('=p3\n') # 12 = 3
		w.write('=p4\n') # 13 = 4
		w.write('=p8\n') # 14 = 8 (A_l)
		w.write('=p5\n') # 15 = 5
		w.write('=p5\n') # 16 = 6
		w.write('=p2\n') # 17 = 2
		w.write('=p3\n') # 18 = 3
		w.write('=p4\n') # 19 = 4
		w.write('0.02 -1\n') # 20 fs
		w.write('=p2\n') # 21 = 2
		w.write('=p3\n') # 22 = 3
		w.write('=p4\n') # 23 = 4	
	
	w.write('set id [open '+wd+'grid_hr-z_gamma'+str(gamma)+'/hr-track_'+str(round(np.log10(nh),2))+'.dat a]\n')
	w.write('flux 0.5 2\n')
	w.write('puts $id "[tcloutr flux 1]"\n') 

	w.write('flux 2 7\n')
	w.write('puts $id "[tcloutr flux 1]"\n')

	w.write('close $id\n')
	w.write('exit\n')
	w.close()
	
wd ='/Users/alberto/Desktop/XBOOTES/HR-z_tracks/'

# INPUT PARAMETERS
z=np.arange(0,5.1,0.1)
NH = np.logspace(20,25,251)
Gamma = 1.8
'''
for j in range(len(NH)):
	
	if NH[j] < 1e22:
		model_to_use = 'pl'
	elif NH[j] <= 1e25:
		model_to_use = 'myt'
	
	for i in range(len(z)):
		
		write_commands(model_to_use,NH[j],z[i],Gamma)
		
		s.call('xspec - '+wd+'hr-track.xcm',shell=True)

sys.exit()
'''

# CT TRACK USING NGC5765B TEMPLATE
'''
for i in range(len(z)):
	w=open('hr-z-CT.xcm','w')
	w.write('@CT_NGC5765B-template.xcm\n')
	w.write('new 8 '+str(z[i])+'\n')

	w.write('set id [open hr-z-CT.dat a]\n')
	w.write('flux 0.5 2\n')
	w.write('puts $id "[tcloutr flux 1]"\n') 

	w.write('flux 2 7\n')
	w.write('puts $id "[tcloutr flux 1]"\n')

	w.write('close $id\n')
	w.write('exit\n')
	w.close()

	s.call('xspec - hr-z-CT.xcm',shell=True)
'''

# TRACK USING MYTORUS
'''
for i in range(len(z)):
	w=open('hr-z-myt.xcm','w')
	w.write('mo pha*(zpow*etable{../../mytorus/mytorus_Ezero_v00.fits}+const*atable{../../mytorus/mytorus_scatteredH500_v00.fits} + const*atable{../../mytorus/mytl_V000010nEp000H500_v00.fits}+const*zpow)\n')
	w.write('0.013 -1\n')	# 1 Nh gal
	w.write(str(Gamma)+' -1\n') # 2 Gamma
	w.write(str(z[i])+' -1\n') # 3 z
	w.write('1e-3 -1\n') # 4 Norm
	w.write(str(NH/1.e24)+' -1\n') # 5 NH
	w.write('90 -1\n') # 6 IncAng
	w.write('=p3\n') # 7 = 3
	w.write('1. -1\n') # 8 A_s
	w.write('=p5\n') # 9 = 5
	w.write('=p6\n') # 10 = 6
	w.write('=p2\n') # 11 = 2
	w.write('=p3\n') # 12 = 3
	w.write('=p4\n') # 13 = 4
	w.write('=p8\n') # 14 = 8 (A_l)
	w.write('=p5\n') # 15 = 5
	w.write('=p5\n') # 16 = 6
	w.write('=p2\n') # 17 = 2
	w.write('=p3\n') # 18 = 3
	w.write('=p4\n') # 19 = 4

	w.write('0.02 -1\n') # 20 fs
	w.write('=p2\n') # 21 = 2
	w.write('=p3\n') # 22 = 3
	w.write('=p4\n') # 23 = 4

	w.write('set id [open hr-z-myt_'+str(Gamma)+'-'+str(int(np.log10(NH)))+'.dat a]\n')
	w.write('flux 0.5 2\n')
	w.write('puts $id "[tcloutr flux 1]"\n') 

	w.write('flux 2 7\n')
	w.write('puts $id "[tcloutr flux 1]"\n')

	w.write('close $id\n')
	w.write('exit\n')
	w.close()
	s.call('xspec - hr-z-myt.xcm',shell=True)
'''

# TRACK USING SIMPLE PL
'''
for i in range(len(z)):
    w=open('hr-z-simple-pl.xcm','w')
    w.write('mo zpow\n')
    
    w.write('1.4 -1\n') # 4 Gamma
    w.write(str(z[i])+' -1\n') # 3 z
    w.write('1e-3 -1\n') # 4 Norm
    
    w.write('set id [open hr-z-simple-pl.dat a]\n')
    w.write('flux 0.5 2\n')
    w.write('puts $id "[tcloutr flux 1]"\n')
    
    w.write('flux 2 7\n')
    w.write('puts $id "[tcloutr flux 1]"\n')
    
    w.write('close $id\n')
    w.write('exit\n')
    w.close()
    s.call('xspec - hr-z-simple-pl.xcm',shell=True)
'''

# TRACK USING OBSCURED PL
'''
for j in range(len(NH)):
	for i in range(len(z)):
		w=open('hr-z-pl.xcm','w')
		w.write('mo pha*zpha*zpow\n')
		w.write('0.0104 -1\n')	# 1 Nh gal
		w.write(str(NH[j]/1.e22)+' -1\n') # 2 NH
		w.write(str(z[i])+' -1\n') # 3 z
		w.write(str(Gamma)+' -1\n') # 4 Gamma
		w.write('=p3\n') # 5 = 3
		w.write('1 -1\n') # 4 Norm

		w.write('set id [open grid_hr-z_gamma1.8/hr-z-pl_'+str(round(np.log10(NH[j]),2))+'.dat a]\n')
		w.write('flux 0.5 2\n')
		w.write('puts $id "[tcloutr flux 1]"\n') 

		w.write('flux 2 7\n')
		w.write('puts $id "[tcloutr flux 1]"\n')

		w.write('close $id\n')
		w.write('exit\n')
		w.close()
		s.call('xspec - hr-z-pl.xcm',shell=True)
'''


cycle=3
plt.figure()
colors = plt.cm.jet(np.linspace(0,1,len(NH)))
soft_cf = 1./4.527000176345934e-12 
hard_cf = 1./1.999000038082175e-11
for j in range(len(NH)):
	
	#ecfs,ecfh = np.genfromtxt('/Users/alberto/Desktop/XBOOTES/cdwfs_ecf_flux-to-cr_CY'+str(cycle)+'.dat',skip_header=1,unpack=True,usecols=[2,3])
	#ecfs=ecfs*1e10
	#ecfh=ecfh*1e10
	#soft_cf = ecfs[int((Gamma-0.9)*10)]
	#hard_cf = ecfh[int((Gamma-0.9)*10)]

	#flux=np.genfromtxt(wd+'grid_hr-z_gamma1.8/hr-z-pl_'+str(round(np.log10(NH[j]),2))+'.dat',unpack=True,usecols=0)
	flux=np.genfromtxt(wd+'grid_hr-z_gamma1.8/hr-track_'+str(round(np.log10(NH[j]),2))+'.dat',unpack=True,usecols=0)
	soft=flux[::2]*soft_cf
	hard=flux[1::2]*hard_cf

	hr=(hard-soft)/(hard+soft)

	plt.plot(z,hr,linestyle='-',color=colors[j])

#flux=np.genfromtxt(wd+'grid_hr-z_gamma1.8/hr-z-pl_22.68.dat',unpack=True,usecols=0)
flux=np.genfromtxt(wd+'grid_hr-z_gamma1.8/hr-track_22.68.dat',unpack=True,usecols=0)
soft=flux[::2]*soft_cf
hard=flux[1::2]*hard_cf
	
hr=(hard-soft)/(hard+soft)
plt.scatter(1.440000057220459, -0.176166)
plt.plot(z,hr,linestyle='-',linewidth=3)
plt.xlabel('Redshift')
plt.ylabel('HR')
plt.tight_layout()
plt.show()

flux_unob=np.genfromtxt(wd+'grid_hr-z_gamma1.8/hr-z-pl_unobsc.dat',unpack=True,usecols=0)
#flux_unob=np.genfromtxt(wd+'grid_hr-z_gamma1.8/hr-track_unobsc.dat',unpack=True,usecols=0)
soft_unob=flux_unob[::2]*soft_cf
hard_unob=flux_unob[1::2]*hard_cf

ks = soft/soft_unob
kh = hard/hard_unob

# More precisely interpolate
zinterp=np.linspace(0.1,5.1,1001)

ksinterp = np.interp(zinterp,z,ks)
khinterp = np.interp(zinterp,z,kh)

print(ksinterp[abs(zinterp-1.440000057220459)==np.min(abs(zinterp-1.440000057220459))])
print(khinterp[abs(zinterp-1.440000057220459)==np.min(abs(zinterp-1.440000057220459))])

plt.plot(z,ks,linestyle='-')
plt.plot(z,kh,linestyle='--')
plt.axvline(x=1.440000057220459,color='k')
#plt.annotate(r'Obscured PL, logNH=22, $\Gamma=1.4-2.0$',xy=(1,0.5))
plt.xlabel('Redshift')
plt.ylabel('K')
plt.tight_layout()
plt.show()

# this code should take the list of detected sources (for simulation, should also match to 
# input sources in order to have input flux), and for each source: 
# define in which obsids it is contained; for each obsid, define theta to estimate r90;
# extract total and bkg counts from circular aperture of radius r90; sum up the total
# net counts and exposure time; compute flux.
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from ciao_contrib.region.check_fov import FOVFiles
import subprocess as s
import sys
import time

def distance(pointa, pointb):
    xx = np.cos(pointa[1]/180.*np.pi)
    return np.sqrt(((pointa[0]-pointb[0])*xx*3600.)**2 +((pointa[1]-pointb[1])*3600.)**2)

wd='/Users/alberto/Desktop/XBOOTES/'

t0 = time.time()
t1 = t0 + 60*58
print('The loop should end around')
print time.strftime("%I %M %p",time.localtime(t1))

t_in=time.time()
'''
obs=np.genfromtxt(wd+'data_counts.dat',usecols=1,dtype='str')
w=open(wd+'fov.lis','w')
for j in range(len(obs)):
	w.write(wd+'data/'+obs[j]+'/repro_new_asol/fov_acisI.fits\n')
w.close()

my_obs = FOVFiles('@'+wd+'fov.lis')
myobs = my_obs.inside(217.9, 33.9)
for j in range(len(myobs)):
	obs=myobs[j][36:-30]
	if len(obs) == 4:
		stem='0'+obs
	elif len(obs) == 3:
		stem='00'+obs
	elif len(obs) == 5:
		stem=obs
	#get theta to compute r90
	path=myobs[j][:-14]+'out/acisf'+stem+'_broad_expomap.fits'
	print(path)
sys.exit()
'''
#simulations: take simulated sources positions and fluxes to match with detected sim srcs
'''
maxdist=3.0

#take input simulated sources
(sim_in_ra,sim_in_dec,sim_in_flux)=np.genfromtxt(wd+'poiss_rand_lehmerx20.dat',unpack=True,usecols=[1,2,0],skip_header=1)

#take detected simulated sources
dat=fits.open(wd+'simul_broad_cat0.fits')
sim_ra=dat[1].data['RA']
sim_dec=dat[1].data['DEC']
sim_sign=dat[1].data['DETML']
sim_tot=dat[1].data['TOT']
sim_bkg=dat[1].data['BKG']
sim_cts=dat[1].data['NET']
sim_exp=dat[1].data['EXP']
sim_cr=dat[1].data['CR']
sim_flux=dat[1].data['FLUX']

out_ra,out_dec,out_inpflux,unmatched_cts,out_multi,out_dist=[],[],[],[],[],[]
multi=0
t_in=time.time()
#w=open(wd+'matched_sim_full.dat','w')
##w.write('RA \t DEC \t sign \t F_in\n')
for i in range(len(sim_ra)):
	multiple=0
	src_det=[sim_ra[i],sim_dec[i]]
	for j in range(len(sim_in_ra)):
		src_inp=[sim_in_ra[j],sim_in_dec[j]]
		dist=distance(src_det,src_inp)
		if (dist < maxdist and multiple==0): # found the first match!
			counter_ra,counter_dec,counter_inpflux,counter_d=[],[],[],[]
			# detected source data
			sim_out_ra=sim_ra[i]
			sim_out_dec=sim_dec[i]
			src_sign=sim_sign[i]
			
			# counterpart (input of simulation) data
			counter_ra.append(sim_in_ra[j])
			counter_dec.append(sim_in_dec[j])
			counter_inpflux.append(sim_in_flux[j])
			counter_d.append(dist)
			
			multiple=multiple+1
			
		elif (dist < maxdist and multiple!=0):
			multiple=multiple+1
			print('Found the '+str(multiple)+'nd/rd/th source...')
			#########################################
			counter_ra.append(sim_in_ra[j])
			counter_dec.append(sim_in_dec[j])
			counter_inpflux.append(sim_in_flux[j])
			counter_d.append(dist)
			
	if multiple !=0:
		multi=multi+1
		# choose the brightest among the counterparts
		if multiple==1: #if there's just one counterpart, pick that
			out_multi.append(1)
			out_ra.append(counter_ra[0])
			out_dec.append(counter_dec[0])
			out_inpflux.append(counter_inpflux[0])
			out_dist.append(counter_d[0])
		else: #we have 2 or more counterparts, pick the brightest
			list_to_sort=[]
			for k in range(len(counter_ra)):
				list_to_sort.append((counter_ra[k],counter_dec[k],counter_inpflux[k],counter_d[k]))
			#sort with criterion of maximum flux; brightest count is the first element
			newlist=sorted(list_to_sort, key=lambda cr: cr[2],reverse=True)
			out_multi.append(1)
			out_ra.append(newlist[0][0])
			out_dec.append(newlist[0][1])
			out_inpflux.append(newlist[0][2])
			out_dist.append(newlist[0][3])
	else:
		out_multi.append(0)
		out_ra.append(-99)
		out_dec.append(-99)
		out_inpflux.append(-99)
		out_dist.append(-99)
		unmatched_cts.append(sim_cts[i])
#w.close()

out_ra=np.array(out_ra)
out_dec=np.array(out_dec)
out_inpflux=np.array(out_inpflux)
out_dist=np.array(out_dist)

print(str(multi)+' out of '+str(len(sim_ra))+' sources were matched.')

out_multi=np.array(out_multi)
sim_ra_y=sim_ra[out_multi==1]
sim_dec_y=sim_dec[out_multi==1]
sim_sign_y=sim_sign[out_multi==1]
sim_tot_y=sim_tot[out_multi==1]
sim_bkg_y=sim_bkg[out_multi==1]
sim_cts_y=sim_cts[out_multi==1]
sim_exp_y=sim_exp[out_multi==1]
sim_cr_y=sim_cr[out_multi==1]
sim_flux_y=sim_flux[out_multi==1]
out_ra_y=out_ra[out_multi==1]
out_dec_y=out_dec[out_multi==1]
out_inpflux_y=out_inpflux[out_multi==1]
out_dist_y=out_dist[out_multi==1]
#write catalog
cat=Table([sim_ra_y,sim_dec_y,sim_sign_y,sim_tot_y,sim_bkg_y,sim_cts_y,sim_exp_y,sim_cr_y,sim_flux_y,out_ra_y,out_dec_y,out_inpflux_y,out_dist_y],names=('RA','DEC','DETML','TOT','BKG','NET','EXP','CR','FLUX','CTRP_RA','CTRP_DEC','CTRP_FLUX','CTRP_DIST'))
cat.write(wd+'simul_broad_matched_cat0.fits',format='fits',overwrite=True) 


print(len(unmatched_cts),'were not matched.')
sim_ra_n=sim_ra[out_multi==0]
sim_dec_n=sim_dec[out_multi==0]
sim_sign_n=sim_sign[out_multi==0]
sim_tot_n=sim_tot[out_multi==0]
sim_bkg_n=sim_bkg[out_multi==0]
sim_cts_n=sim_cts[out_multi==0]
sim_exp_n=sim_exp[out_multi==0]
sim_cr_n=sim_cr[out_multi==0]
sim_flux_n=sim_flux[out_multi==0]
out_ra_n=out_ra[out_multi==0]
out_dec_n=out_dec[out_multi==0]
out_inpflux_n=out_inpflux[out_multi==0]
out_dist_n=out_dist[out_multi==0]
#write catalog of unmatched sources
cat=Table([sim_ra_n,sim_dec_n,sim_sign_n,sim_tot_n,sim_bkg_n,sim_cts_n,sim_exp_n,sim_cr_n,sim_flux_n,out_ra_n,out_dec_n,out_inpflux_n,out_dist_n],names=('RA','DEC','DETML','TOT','BKG','NET','EXP','CR','FLUX','CTRP_RA','CTRP_DEC','CTRP_FLUX','CTRP_DIST'))
cat.write(wd+'simul_broad_unmatched_cat0.fits',format='fits',overwrite=True) 

plt.figure()
plt.hist(unmatched_cts,bins=int(np.max(unmatched_cts)))
plt.show()

dat.close()
sys.exit()
#############################
'''
sim_out_ra,sim_out_dec,sim_out_sig,sim_out_inpflux=np.genfromtxt(wd+'matched_sim_full.dat', unpack=True, skip_header=1)


inp_cts,out_cts,out_flux=[],[],[]
#now, for each detected source compute flux from data and bkg maps
for i in range(len(sim_out_ra)):
	print(i+1,len(sim_out_ra))
	#see in which obsids this source is contained 
	my_obs = FOVFiles('@'+wd+'fov.lis')
	myobs = my_obs.inside(sim_out_ra[i], sim_out_dec[i])
	cts,bkg,exp=[],[],[]
	for j in range(len(myobs)):
		#get theta to compute r90
		obs=myobs[j][36:-30]
		#print(obs)
		if len(obs) == 4:
			stem='0'+obs
		elif len(obs) == 3:
			stem='00'+obs
		elif len(obs) == 5:
			stem=obs
		path=myobs[j][:-14]+'out/acisf'+stem+'_broad_expomap.fits'
		s.call('punlearn dmcoords',shell=True)
		s.call('dmcoords '+path+' asol=none opt=cel celfmt=deg ra='+str(sim_out_ra[i])+' dec='+str(sim_out_dec[i])+'',shell=True)
		res=s.check_output('pget dmcoords theta',shell=True)
		theta=res.splitlines()
		theta=float(theta[0])
		r90=1+10*(theta/10.)**2
		#print(r90)
		
		#extract counts and background
		imagemap=wd+'sim_full/acisf'+stem+'_sim_poiss_bitpix-32.fits'
		backmap=wd+'data/'+obs+'/repro_new_asol/out/acisf'+stem+'_broad_bkgmap.fits'
		s.call('pset dmextract infile="'+imagemap+'[bin pos=circle('+str(sim_out_ra[i])+'d,'+str(sim_out_dec[i])+'d,'+str(r90*0.000277778)+'d)]" mode=h bkg="'+backmap+'[bin pos=circle('+str(sim_out_ra[i])+'d,'+str(sim_out_dec[i])+'d,'+str(r90*0.000277778)+'d)]" outfile=counts.fits opt=generic mode=h clobber=yes',shell=True)
		s.call('dmextract',shell=True)
		s.call('dmlist "counts.fits[cols COUNTS, BG_COUNTS]" data,clean | grep -v COUNTS > counts.dat', shell=True)
		(cts_i,bkg_i)=np.genfromtxt('counts.dat',unpack=True)
		cts.append(cts_i)
		bkg.append(bkg_i)
		
		#extract exposure from vignetting-corrected expomap
		expomap=wd+'data/'+obs+'/repro_new_asol/out/acisf'+stem+'_broad_expomap.fits'
		s.call('dmextract infile="'+expomap+'[bin pos=circle('+str(sim_out_ra[i])+'d,'+str(sim_out_dec[i])+'d,'+str(r90*0.000277778)+'d)]" mode=h outfile=expo.fits opt=generic mode=h clobber=yes',shell=True)
		s.call('dmlist "expo.fits[cols COUNTS, AREA]" data,clean | grep -v COUNTS > expo.dat',shell=True)
		(totexpo,npix)=np.genfromtxt('expo.dat',unpack=True)
		exp_i=totexpo/npix
		exp.append(exp_i)
	
	net=np.sum(np.array(cts))-np.sum(np.array(bkg))
	count_rate=net*1.1/np.sum(np.array(exp))
	flux=count_rate*1.825E-11
	out_flux.append(flux)
	out_cts.append(net*1.1)
	inp_cts.append(sim_out_inpflux[i]*np.sum(np.array(exp))*5.438E+10)
	#print(cts)
	#print(bkg)
	#print(exp)
	#print(sim_out_inpflux[i],flux,count_rate)
	#sys.exit()
t_out=time.time()
print((t_out-t_in)/60.,'Minutes for the loop')
'''
out_flux=np.array(out_flux)

x=np.logspace(np.log10(1e-17),np.log10(1e-11),30)
plt.figure()
plt.plot(x,x,'k--')
plt.plot(sim_out_inpflux,out_flux,'g.')
#plt.plot(sim_out_inpflux0,out_flux0,'r.',legend='<0.5"')
#plt.plot(sim_out_inpflux1,out_flux1,'g.',legend='<1.0"')
#plt.plot(sim_out_inpflux2,out_flux2,'c.',legend='<1.5"')
#plt.plot(sim_out_inpflux3,out_flux3,'b.',legend='<2.0"')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('F_in (cgs)')
plt.ylabel('F_out (cgs)')
plt.tight_layout()
plt.show()
'''
x=np.linspace(0,300,300)
plt.figure()
plt.plot(x,x,'k--')
plt.plot(inp_cts,out_cts,'g.')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Cts_in')
plt.ylabel('Cts_out')
plt.tight_layout()
plt.show()
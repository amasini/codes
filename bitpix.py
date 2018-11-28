# DMEXTRACT GIVES PROBLEMS WHEN TRYING TO EXTRACT COUNTS FROM SIMULATIONS, BECAUSE FOR SOME
# REASON THEY HAVE A BITPIX KEYWORD IN THEIR HEADER +64 INSTEAD OF -32.
# THIS HAPPENS AT LINE #162 (poiss_sim=np.random.poisson(simulation)) IN SIMUL.PY AND 
# THERE'S NO WAY TO PREVENT THIS - THAT'S WHY I MADE THIS SCRIPT, THAT INSTEAD OF CALLING
# np.random.poisson FOR THE WHOLE FILE calls it pixel by pixel. THIS CREATES ANOTHER 
# PROBLEM WITH XIMAGE, THOUGH. WEIRD MOSAICS RESULT IN THIS WAY.
import numpy as np
import sys
from astropy.io import fits

wd='/Users/alberto/Desktop/XBOOTES/'
#obsid=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=1,dtype='str')
obsid=['6998']
for k in range(len(obsid)):
	if len(obsid[k]) == 4:
		stem='0'+obsid[k]
	elif len(obsid[k]) == 3:
		stem='00'+obsid[k]
	elif len(obsid[k]) == 5:
		stem=obsid[k]

	sim=fits.open(wd+'sim_full/acisf'+stem+'_sim.fits')   
	simulation=sim[0].data
	backheader=sim[0].header
	
	print(obsid[k])
	
	print(simulation.dtype.name)
	
	#poiss_sim=np.zeros_like(simulation, dtype=float)
	#for i in range(len(simulation)):
	#	for j in range(len(simulation[0])):
	#		poiss_sim[i][j]=np.random.poisson(simulation[i][j])
	
	#old method
	poiss_sim2=np.random.poisson(simulation)
	############# this solves the DMEXTRACT issue, getting from int64 to float66
	poiss_sim3=poiss_sim2.astype(float)
	#####################################
	print(poiss_sim2.dtype.name)
	print(poiss_sim3.dtype.name)

	hdu1 = fits.PrimaryHDU(poiss_sim3,header=backheader)
	hdu1.writeto(wd+'PROVA_acisf'+stem+'_sim_poiss.fits',overwrite=True)
	
	sim.close()
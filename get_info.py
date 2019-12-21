import numpy as np
import subprocess as s
import sys

wd = '/Users/alberto/Desktop/XBOOTES/'

(obs, cy) = np.genfromtxt(wd+'data_counts.dat', unpack = True, usecols = [1, 6], dtype = 'str')

w =open(wd+'data_info.txt', 'w')

for i in range(len(obs)):
	if len(obs[i]) == 3:
		stem = '00'+obs[i]
	elif len(obs[i]) == 4:
		stem = '0'+obs[i]
	else:
		stem = obs[i]
	
	path = wd+'data/'+obs[i]+'/repro_new_asol/'
	filename = 'acisf'+stem+'_repro_evt2.fits'
	
	what = ['RA_PNT', 'DEC_PNT', 'ROLL_PNT', 'EXPOSURE']
	precision = [3, 3, 1, 1]
	vals = []
	for j in range(len(what)):
		
		command = 'dmkeypar '+path+filename+' '+what[j]+' echo+'
	
		out = s.check_output(command, shell= True)
		if what[j] == 'EXPOSURE':
			out = float(out)/1e3
		
		vals.append(round(float(out), precision[j]))
	
	w.write(obs[i]+', '+str(vals)+', '+cy[i]+'\n')

w.close()

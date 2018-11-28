import numpy as np
import sys
import subprocess as s

o=open('/home/alberto/Downloads/browse_download_script.txt','r')
for _ in range(6):
	next(o)
for l in o:
	l=l.strip()
	c=l.split()
	if len(l) != 0:
		if c[0] != '#':
			#print(l)
			s.call(l,shell=True)
o.close()
s.call('rm -f /home/alberto/Downloads/browse_download_script.txt', shell=True)

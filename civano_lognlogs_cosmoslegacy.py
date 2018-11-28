import numpy as np
import sys
import matplotlib.pyplot as plt

wd='/Users/alberto/Desktop/'

(f,ns0)=np.genfromtxt(wd+'civano_2-10keV_lognlogs.dat',unpack=True,skip_header=1)

ns=ns0/(f/1e-14)**1.5

w=open(wd+'XBOOTES/lognlogs_civano_cosmoslegacy.dat','w')
w.write('#2-10keVflux \t N(>S)\n')
for i in range(len(f)):
    w.write(str(f[i])+'\t'+str(ns[i])+'\n')
w.close()

plt.figure()
plt.plot(f,ns,'r--')
plt.xscale('log')
plt.yscale('log')
plt.show()

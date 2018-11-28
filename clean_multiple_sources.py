# THIS SCRIPT REMOVES DUPLICATES IN A LIST WITH A USER DEFINED MAX DISTANCE AND RETAINS THE
# BEST POINT BASED ON A CRITERION - IN THIS CASE THE LOWEST PROB (of being spurious)
import numpy as np
import datetime
from itertools import combinations
import sys
from astropy.io import fits
from astropy.table import Table
#import subprocess
##from sherpa.astro.ui import *
#from pycrates import *

def distance(pointa, pointb):
    """
    Function that calculates the distance
    It is here because it needs to be called many many times
    """
    xx = np.cos(pointa[0][1]/180*3.141592)
    return np.sqrt(((pointa[0][0]-pointb[0][0])*3600*xx)**2 +((pointa[0][1]-pointb[0][1])*3600)**2)

def choose_best(pointa, pointb):
    """
    Function that chooses the best between two points
    This version is based on the lower prob of being spurious
    """
    if pointa[1] <= pointb[1]:
        return pointa, pointb
    else:
        return pointb, pointa 
    
def unify_list(mylist):
    """
    function that given a list of points, removes the points that at a distance 
    less than maxdist
 
    WARNING: this function need to be tested
    """
    #create list of synonims
    nearby_points = {}
    #for each couple of elements I check if their distance is less than my limit
    for couple in combinations(mylist, 2):
        if (distance(couple[0], couple[1]) <= 1.1*couple[0][0][3] or distance(couple[0], couple[1]) <= 1.1*couple[1][0][3]):
            #choose which point is the main and which is the secondary
            main_point, secondary_point = choose_best(couple[0], couple[1])
            #if the distance is less than the limit, I assign the second as synonim of the first
            nearby_points.setdefault(main_point, []).append(secondary_point)
    #at the end of the previous for, what I have is the following:
    #if I have A, B and C and they are all under the limit
    #then I will have that nearby_points is {A:[B, C]}
    
    #create inverted dictionary of synonims: I need it to have a fast lookup
    inverted_syn = {}
    for key, syn_list in nearby_points.iteritems():
        for synonim in syn_list:
            inverted_syn[synonim] = key
    #what I will have after this section is that for each "secondary" 
    #point I assign to it the main one       
    #basically, if I have nearby_points == {A:[B, C]}
    #I will have inverted_syn = {B:A, C:A}
    
    #reconstruct the list
    ret_dict = {}
    for elem in mylist:
        if elem not in inverted_syn:
            ret_dict[elem] = True
        else:
            ret_dict[inverted_syn[elem]] = True
    #what the previous code does is:
    #it takes every element of the original list and it checks if 
    #it is a synonim of a more "important" one (meaning if it is in inverted_syn
    #if it is not in inverted_syn, then it means that it has no synonims and it can be used as it is
    #otherwise the synonim is taken instead of the original
    #in other words if mylist is [A, B, C, D]
    #D is not in inverted_syn and is taken as it is
    #A is not in inverted_syn (it is the main one), and is taken as it is
    #for B and C I take A instead of them
    #I use ret_dict instead of a list because it unifies the results
            
    ret_list = ret_dict.keys()
    ret_list.sort()
    if ret_list == mylist:
        return ret_list
    else:
        return unify_list(ret_list)
    #the previous code takes all the keys of ret_dict and put them in a list
    #then it checks if the new list is the same of the one in input
    #if it is than we are done
    #otherwise we call the same function recursively to run the reduction on the new set


def main():
    """main function"""

wd='/Users/alberto/Desktop/XBOOTES/'

tstart= datetime.datetime.now()
#read the file and put it in memory
file_in_memory = []

hdul=fits.open(wd+'mosaic_broad_sim_poiss_cat0_3.fits') #file to clean, with duplicates (blendings)
ra_k=hdul[1].data['RA']
dec_k=hdul[1].data['DEC']
detml=hdul[1].data['DETML']
r90=hdul[1].data['AV_R90']
cts_full=hdul[1].data['TOT']
flux_k=hdul[1].data['FLUX']
	
file=hdul[1].data
hdul.close()
#cut detected sample at a given probability threshold (1e-3)
cut=1e-3
prob=np.e**(-detml)

ra_k=ra_k[prob<=cut]
dec_k=dec_k[prob<=cut]
cts_full=cts_full[prob<=cut]
r90=r90[prob<=cut]
detml=detml[prob<=cut]
flux_k=flux_k[prob<=cut]
file=file[prob<=cut]
prob=prob[prob<=cut]

for line in range(len(file)):
	file_in_memory.append((file[line],prob[line]))
	
#file_in_memory[0][0] is the point (line from catalog), file_in_memory[0][1] is prob

print('Starting with ',len(file_in_memory))

print('calling function to remove duplicates...')
#call the function to remove the points too close within (1.1*R90)"
final_list = unify_list(file_in_memory)

print('Ended with ',len(final_list))

file_in_memory2=[]
for line in range(len(final_list)):
	line_in_memory=[]
	#file_in_memory.append((file[line],theta[line]))
	for j in range(len(final_list[line][0])):
		line_in_memory.append(final_list[line][0][j])
	line_in_memory.append(final_list[line][1])
	file_in_memory2.append(line_in_memory)

ra,dec,prob,av_r90,tot,bkg,net,exptime,cr,flux=[],[],[],[],[],[],[],[],[],[]
for i in range(len(file_in_memory2)):
	ra.append(file_in_memory2[i][0])
	dec.append(file_in_memory2[i][1])
	
	av_r90.append(file_in_memory2[i][3])
	tot.append(file_in_memory2[i][4])
	bkg.append(file_in_memory2[i][5])
	net.append(file_in_memory2[i][6])
	exptime.append(file_in_memory2[i][7])
	cr.append(file_in_memory2[i][8])
	flux.append(file_in_memory2[i][9])
	
	prob.append(file_in_memory2[i][10])

#write catalog
cat=Table([ra,dec,prob,av_r90,tot,bkg,net,exptime,cr,flux],names=('RA','DEC','PROB','AV_R90','TOT','BKG','NET','EXP','CR','FLUX'))
cat.write(wd+'mosaic_broad_sim_poiss_cat1_3.fits',format='fits',overwrite=True)

tstop= datetime.datetime.now()
time_elapse = tstop - tstart
print 'timer ',divmod(time_elapse.days * 86400 + time_elapse.seconds, 60)		    

if __name__ == '__main__':
	main()

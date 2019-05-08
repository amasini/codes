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
    xx = np.cos(pointa[1]/180*3.141592)
    return np.sqrt(((pointa[0]-pointb[0])*3600*xx)**2 +((pointa[1]-pointb[1])*3600)**2)

def choose_best(pointa, pointb):
	"""
	Function that chooses the best between two points
	This version is based on the lower prob of being spurious
	"""
	#if pointa[1] <= pointb[1]:
	if pointa[2] <= pointb[2]:
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
		if (distance(couple[0], couple[1]) <= 1.1*couple[0][3] or distance(couple[0], couple[1]) <= 1.1*couple[1][3]):
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

band='hard'

tstart= datetime.datetime.now()
#read the file and put it in memory
file_in_memory = []

hdul=fits.open(wd+'cdwfs_'+band+'_cat0.fits') #file to clean, with duplicates (blendings)

prob=hdul[1].data['PROB']
file=hdul[1].data

hdul.close()

'''
#cut detected sample at a given probability threshold (1e-3)
cut=1e-3
#prob=np.e**(-detml)

ra_k=ra_k[prob<=cut]
dec_k=dec_k[prob<=cut]
cts_full=cts_full[prob<=cut]
r90=r90[prob<=cut]
#detml=detml[prob<=cut]
flux_k=flux_k[prob<=cut]
file=file[prob<=cut]
prob=prob[prob<=cut]
'''
for line in range(len(file)):
	#file_in_memory.append((file[line],prob[line]))
	file_in_memory.append(file[line])

#file_in_memory[0][0] is the point (line from catalog), file_in_memory[0][1] is prob

print('Starting with ',len(file_in_memory))

print('calling function to remove duplicates...')
#call the function to remove the points too close within (1.1*R90)"
final_list = unify_list(file_in_memory)

print('Ended with ',len(final_list))

ra,dec,prob,av_r90,tot,bkg,net,e_net_up,e_net_lo,exptime,cr,e_cr_up,e_cr_lo,flux,e_flux_up,e_flux_lo=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
for i in range(len(final_list)):
	ra.append(final_list[i][0])
	dec.append(final_list[i][1])
	
	prob.append(final_list[i][2])
	av_r90.append(final_list[i][3])
	tot.append(final_list[i][4])
	bkg.append(final_list[i][5])
	net.append(final_list[i][6])
	e_net_up.append(final_list[i][7])
	e_net_lo.append(final_list[i][8])
	exptime.append(final_list[i][9])
	cr.append(final_list[i][10])
	e_cr_up.append(final_list[i][11])
	e_cr_lo.append(final_list[i][12])
	flux.append(final_list[i][13])
	e_flux_up.append(final_list[i][14])
	e_flux_lo.append(final_list[i][15])	

#write catalog
cat=Table([ra,dec,prob,av_r90,tot,bkg,net,e_net_up,e_net_lo,exptime,cr,e_cr_up,e_cr_lo,flux,e_flux_up,e_flux_lo],names=('RA','DEC','PROB','AV_R90','TOT','BKG','NET','E_NET_+','E_NET_-','EXP','CR','E_CR_+','e_CR_-','FLUX','E_FLUX_+','E_FLUX_-'))
cat.write(wd+'cdwfs_'+band+'_cat1.fits',format='fits',overwrite=True)

tstop= datetime.datetime.now()
time_elapse = tstop - tstart
print 'timer ',divmod(time_elapse.days * 86400 + time_elapse.seconds, 60)		    

if __name__ == '__main__':
	main()

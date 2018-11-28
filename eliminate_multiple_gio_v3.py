# THIS SCRIPT REMOVES DUPLICATES IN A LIST WITH A USER DEFINED MAX DISTANCE AND RETAINS THE
# BEST POINT BASED ON A CRITERION - IN THIS CASE THE LOWEST THETA (SEP. FROM AIMPOINT)
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
    This version is based on the lower theta
    """
    if pointa[1] <= pointb[1]:
        return pointa, pointb
    else:
        return pointb, pointa 
    
def unify_list(mylist, maxdist):
    """
    function that given a list of points, removes the points that at a distance 
    less than maxdist
 
    WARNING: this function need to be tested
    """
    #create list of synonims
    nearby_points = {}
    #for each couple of elements I check if their distance is less than my limit
    for couple in combinations(mylist, 2):
        if 0 <= distance(couple[0], couple[1]) <= maxdist:
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
        return unify_list(ret_list, maxdist)
    #the previous code takes all the keys of ret_dict and put them in a list
    #then it checks if the new list is the same of the one in input
    #if it is than we are done
    #otherwise we call the same function recursively to run the reduction on the new set
    


def main():
    """main function"""

wd='/Users/alberto/Desktop/XBOOTES/'
obsid=np.genfromtxt(wd+'data_counts.dat',unpack=True,usecols=1,dtype='str')

tstart= datetime.datetime.now()
#read the file and put it in memory
file_in_memory = []
for i in range(len(obsid)):
	if len(obsid[i]) == 4:
		stem='0'+obsid[i]
	elif len(obsid[i]) == 3:
		stem='00'+obsid[i]
	elif len(obsid[i]) == 5:
		stem=obsid[i]
	
	hdul=fits.open(wd+'/wav_xbootes_singleacis_5e-5/acisf'+stem+'_src_new_3.fits')
	expofile=fits.open(wd+'data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img')
	exp=expofile[0].data
	expheader=expofile[0].header
	crpix1=expheader["CRPIX1"] #aimpoint pixel in x 1024.5 -> 4096.5
	crpix1_phys=expheader["CRPIX1"]+3072 #aimpoint pixel in x 1024.5 -> 4096.5
	crpix2=expheader["CRPIX2"] #aimpoint pixel in y 1024.5 -> 4096.5
	crpix2_phys=expheader["CRPIX2"]+3072 #aimpoint pixel in y 1024.5 -> 4096.5
		
	file=hdul[1].data
	expofile.close()
	hdul.close()
	
	x=file["X"] #these are phys quantities
	y=file["Y"]
	
	deltax_phys=crpix1_phys-x
	deltay_phys=crpix2_phys-y
	deltax_ima=deltax_phys/4.
	deltay_ima=deltay_phys/4.
	x_ima=crpix1-deltax_ima
	y_ima=crpix2-deltay_ima
		
	theta=np.sqrt(((x_ima-crpix1)**2+(y_ima-crpix2)**2)*(1.968/60.)**2)

	for line in range(len(file)):
		file_in_memory.append((file[line],theta[line]))
	
	'''
	with fits.open(wd+'/wav_xbootes_singleacis_5e-5/acisf'+stem+'_src_new_3.fits') as hdul:
		expofile=fits.open(wd+'data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img')
		exp=expofile[0].data
		expheader=expofile[0].header
		crpix1=expheader["CRPIX1"] #aimpoint pixel in x 1024.5 -> 4096.5
		crpix1_phys=expheader["CRPIX1"]+3072 #aimpoint pixel in x 1024.5 -> 4096.5
		crpix2=expheader["CRPIX2"] #aimpoint pixel in y 1024.5 -> 4096.5
		crpix2_phys=expheader["CRPIX2"]+3072 #aimpoint pixel in y 1024.5 -> 4096.5
		
		file=hdul[1].data
		
		#ra=file["RA"]
		#dec=file["DEC"]
		#ra_err=file["RA_ERR"]
		#dec_err=file["DEC_ERR"]
		x=file["X"] #these are phys quantities
		y=file["Y"]
		#x_err=file["X_ERR"]
		#y_err=file["Y_ERR"]
		#npixsou=file["NPIXSOU"]
		#net_counts=file["NET_COUNTS"]
		#net_counts_err=file["NET_COUNTS_ERR"]
		#bkg_counts=file["BKG_COUNTS"]
		#bkg_counts_err=file["BKG_COUNTS_ERR"]
		#net_rate=file["NET_RATE"]
		#net_rate_err=file["NET_RATE_ERR"]
		#bkg_rate=file["BKG_RATE"]
		#bkg_rate_err=file["BKG_RATE_ERR"]
		#exptime=file["EXPTIME"]
		#exptime_err=file["EXPTIME_ERR"]
		#src_significance=file["SRC_SIGNIFICANCE"]
		#psf_size=file["PSF_SIZE"]
		#multi_correl_max=file["MULTI_CORREL_MAX"]
		#shape=file["SHAPE"]
		#r=file["R"]
		#rotang=file["ROTANG"]
		#psfratio=file["PSFRATIO"]
		#component=file["COMPONENT"]
		
		deltax_phys=crpix1_phys-x
		deltay_phys=crpix2_phys-y
		deltax_ima=deltax_phys/4.
		deltay_ima=deltay_phys/4.
		x_ima=crpix1-deltax_ima
		y_ima=crpix2-deltay_ima
		
		theta=np.sqrt(((x_ima-crpix1)**2+(y_ima-crpix2)**2)*(1.968/60.)**2)

	for line in range(len(file)):
		file_in_memory.append((file[line],theta[line]))
	'''
#file_in_memory[0][0] is the point (line from catalog), file_in_memory[0][1] is theta

'''
with open(''+d+'/Output_'+fpm+'_'+energy+'_prob_convALL_resorted.list', 'r') as m:
	for linem in m:
		linem = linem.strip()
		columnsm = linem.split()
		if len(columnsm)>0:
		#transform the data to be usable
			id_sex=int(columnsm[0])
			probm=float(columnsm[1])
			pix=float(columnsm[2])
			ram=float(columnsm[3])
			decm=float(columnsm[4])
			ra2=float(columnsm[5])
			dec2=float(columnsm[6])
			prob_pois=columnsm[7]
			detm=float(columnsm[8])
			counts=columnsm[9]
			bg_counts=columnsm[10]
			tot=columnsm[11]
			bg=columnsm[12]
			net=columnsm[13]
			#create a new version of the row and put it in the list in memory
			new_row = (id_sex, probm, pix, ram, decm, ra2, dec2, prob_pois, detm, counts, bg_counts,tot,bg,net)
			file_in_memory.append(new_row)
'''
print('Starting with ',len(file_in_memory))

print('calling function to remove duplicates...')
#call the function to remove the points too close within 2"
final_list = unify_list(file_in_memory, 2.0)

print('Ended with ',len(final_list))

file_in_memory2=[]
for line in range(len(final_list)):
	line_in_memory=[]
	#file_in_memory.append((file[line],theta[line]))
	for j in range(len(final_list[line][0])):
		line_in_memory.append(final_list[line][0][j])
	line_in_memory.append(final_list[line][1])
	file_in_memory2.append(line_in_memory)

ra,dec,ra_err,dec_err,x,y,x_err,y_err,npixsou,net_counts,net_counts_err,bkg_counts,bkg_counts_err,net_rate,net_rate_err,bkg_rate,bkg_rate_err,exptime,exptime_err,src_significance,psf_size,multi_correl_max,shape,r,rotang,psfratio,component,theta=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
for i in range(len(file_in_memory2)):
	ra.append(file_in_memory2[i][0])
	dec.append(file_in_memory2[i][1])
	ra_err.append(file_in_memory2[i][2])
	dec_err.append(file_in_memory2[i][3])
	x.append(file_in_memory2[i][4])
	y.append(file_in_memory2[i][5])
	x_err.append(file_in_memory2[i][6])
	y_err.append(file_in_memory2[i][7])
	npixsou.append(file_in_memory2[i][8])
	net_counts.append(file_in_memory2[i][9])
	net_counts_err.append(file_in_memory2[i][10])
	bkg_counts.append(file_in_memory2[i][11])
	bkg_counts_err.append(file_in_memory2[i][12])
	net_rate.append(file_in_memory2[i][13])
	net_rate_err.append(file_in_memory2[i][14])
	bkg_rate.append(file_in_memory2[i][15])
	bkg_rate_err.append(file_in_memory2[i][16])
	exptime.append(file_in_memory2[i][17])
	exptime_err.append(file_in_memory2[i][18])
	src_significance.append(file_in_memory2[i][19])
	psf_size.append(file_in_memory2[i][20])
	multi_correl_max.append(file_in_memory2[i][21])
	shape.append(file_in_memory2[i][22])
	r.append(file_in_memory2[i][23])
	rotang.append(file_in_memory2[i][24])
	psfratio.append(file_in_memory2[i][25])
	component.append(file_in_memory2[i][26])
	theta.append(file_in_memory2[i][27])

#write catalog
cat=Table([ra,dec,ra_err,dec_err,x,y,x_err,y_err,npixsou,net_counts,net_counts_err,bkg_counts,bkg_counts_err,net_rate,net_rate_err,bkg_rate,bkg_rate_err,exptime,exptime_err,src_significance,psf_size,multi_correl_max,shape,r,rotang,psfratio,component,theta],names=('RA','DEC','RA_ERR','DEC_ERR','X','Y','X_ERR','Y_ERR','NPIXSOU','NET_COUNTS','NET_COUNTS_ERR', 'BKG_COUNTS','BKG_COUNTS_ERR','NET_RATE','NET_RATE_ERR','BKG_RATE','BKG_RATE_ERR','EXPTIME','EXPTIME_ERR','SRC_SIGNIFICANCE','PSF_SIZE','MULTI_CORREL_MAX','SHAPE','R','ROTANG','PSFRATIO','COMPONENT','THETA'))
cat.write(wd+'murray_sens/cdwfs_unique_3.fits',format='fits',overwrite=True)

'''
with open(''+d+'/Output_'+fpm+'_'+energy+'_prob_nodouble30.list', 'a') as f:  #remember to remove this after each time you run the script
	for elem in final_list:
		line = '%s %s %s %s %s %s %s %s %s %s %s  %s %s %s\n' % (elem[0], elem[1], elem[2], elem[3], elem[4], elem[5], elem[6], elem[7], elem[8],  elem[9],  elem[10],  elem[11],  elem[12],  elem[13])
		f.write(line)
'''
tstop= datetime.datetime.now()
time_elapse = tstop - tstart
print 'timer ',divmod(time_elapse.days * 86400 + time_elapse.seconds, 60)		    

if __name__ == '__main__':
	main()

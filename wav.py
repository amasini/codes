import numpy as np
import subprocess as s
import sys

#wd='/Users/alberto/Desktop/XBOOTES/chunks_of_mosaics_fullres'
wd='/Users/alberto/Desktop/XBOOTES/murray_sens/'
#wd='/Users/alberto/Desktop/XBOOTES/test_merge/'
#wd='/Users/alberto/Desktop/XBOOTES/data/'

'''
# WAVDETECT ON CHUNKS OF SIMULATIONS
for chunk_number in range(0,6):
	s.call('punlearn wavdetect', shell=True)
	infile=wd+'/mosaic_broad_chunk_'+str(chunk_number)+'_new.fits'
	bkginput=wd+'/mosaic_broad_chunk_'+str(chunk_number)+'_bkg.fits'
	
	psffile=''
	#psffile=wd+'/mosaic_broad_chunk_'+str(chunk_number)+'_averager90.fits'
	
	outfile=wd+'/mosaic_'+str(chunk_number)+'_src_new_bothscales.fits'
	scellfile=wd+'/mosaic_'+str(chunk_number)+'_scell_new_bothscales.fits'
	imagefile=wd+'/mosaic_'+str(chunk_number)+'_ima_new_bothscales.fits'
	defnbkgfile=wd+'/mosaic_'+str(chunk_number)+'_defnbkg_new_bothscales.fits'
	regfile=wd+'/mosaic_'+str(chunk_number)+'_src_new_bothscales.reg'
	#scales='\"2.0 4.0 8.0 16.0\"' #large scales
	#scales='\"1.0 1.4142 2.0 4.0\"' #small scales
	scales='\"1.0 2.0 4.0 8.0 16.0\"'
	sigthresh='1e-6'
	clobber='yes'
	
	s.call('wavdetect infile='+infile+' bkginput='+bkginput+' psffile='+psffile+' outfile='+outfile+' scellfile='+scellfile+' imagefile='+imagefile+' defnbkgfile='+defnbkgfile+' regfile='+regfile+' scales='+scales+' sigthresh='+sigthresh+' clobber='+clobber+'',shell=True)
'''
# WAVDETECT ON TOTAL MOSAIC 4X4 REBINNED
#wavdetect infile=mosaic_broad_4rebinned.fits bkginput=mosaic_broad_bkgmap_4rebinned_murray.fits psffile= outfile=mosaic_broad_src_3.fits scellfile=out_scell.fits imagefile=out_ima.fits defnbkgfile=out_defnbkg.fits regfile=mosaic_broad_src_3.reg scales="1 2 4" sigthresh=5e-5 clobber=yes

# WAVDETECT ON ALL 126 XBOOTES ACIS-I 4X4 BINNED IMAGES
di='/Users/alberto/Desktop/XBOOTES/'
########## #choose which scales use
scales='\"1 2 4\"'
obsid=np.genfromtxt(di+'data_counts.dat',unpack=True,usecols=1,dtype='str')
#obsid=['4254']
for band in ['broad']:
	for i in range(len(obsid)):
		if len(obsid[i]) == 4:
			stem='0'+obsid[i]
		elif len(obsid[i]) == 3:
			stem='00'+obsid[i]
		elif len(obsid[i]) == 5:
			stem=obsid[i]
		
		s.call('punlearn wavdetect', shell=True)
		infile=di+'/data/'+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_05to7keV_4rebinned.img'
		bkginput=di+'/data/'+obsid[i]+'/repro_new_asol/out/acisf'+stem+'_broad_bkgmap_4rebinned.img'
		#psffile=di+'murray_psfmaps/'+stem+'_r90_4rebinned.fits'
		psffile=''

		# in the output directory, _src.fits is the first WITH NO PSF AND NO BKG
		# _src_new.fits HAS ONLY PSF
		# _src_new_2.fits HAS BOTH PSF AND BKG
		# _src_new_3.fits ONLY BKG
		outfile=di+'/wav_xbootes_singleacis_5e-5/acisf'+stem+'_src_new_3.fits'
		scellfile=di+'/wav_xbootes_singleacis_5e-5/acisf'+stem+'_scell_new_3.fits'
		imagefile=di+'/wav_xbootes_singleacis_5e-5/acisf'+stem+'_ima_new_3.fits'
		defnbkgfile=di+'/wav_xbootes_singleacis_5e-5/acisf'+stem+'_defnbkg_new_3.fits'
		regfile=di+'/wav_xbootes_singleacis_5e-5/acisf'+stem+'_src_new_3.reg'
		sigthresh='5e-5'
		clobber='yes'
	
		print('wavdetecting '+obsid[i]+'...')
		s.call('wavdetect infile='+infile+' bkginput='+bkginput+' psffile='+psffile+' outfile='+outfile+' scellfile='+scellfile+' imagefile='+imagefile+' defnbkgfile='+defnbkgfile+' regfile='+regfile+' scales='+scales+' sigthresh='+sigthresh+' clobber='+clobber+'',shell=True)

sys.exit()

# WAVDETECT ON CHUNKS OF REAL DATA
########## #choose which scales use
sc=''
###########
if sc=='large':
	scales='\"2.0 4.0 8.0 16.0\"'
	name=''
elif sc=='small':
	scales='\"1.0 1.4142 2.0 4.0\"'
	name='_smallscales'
elif sc=='both':
	scales='\"1.0 2.0 4.0 8.0 16.0 32.0\"'
	name='_bothscales'
name=''	
scales='\"1 2 4 8\"'
for band in ['broad']:
	for chunk_number in range(5,6):
		s.call('punlearn wavdetect', shell=True)
		infile=wd+'/mosaic_'+band+'_chunk_'+str(chunk_number)+'_dat_4rebinned.fits'
		#bkginput=wd+'/mosaic_'+band+'_chunk_'+str(chunk_number)+'_bkg_4rebinned.fits'
		bkginput=''
		#infile=wd+'/mosaic_'+band+'_4rebinned_not-expo-corr.fits'
		#bkginput=wd+'/mosaic_'+band+'_bkgmap_4rebinned.fits'
		
		psffile=''
		#psffile=wd+'/mosaic_broad_chunk_'+str(chunk_number)+'_averager90.fits'
	
		outfile=wd+'/mosaic_'+band+'_'+str(chunk_number)+'_src'+name+'_4reb.fits'
		#outfile=wd+'/mosaic_'+band+'_4rebinned_src.fits'
		scellfile=wd+'/mosaic_'+band+'_'+str(chunk_number)+'_scell'+name+'_4reb.fits'
		#scellfile=wd+'/deleteme2.fits'
		imagefile=wd+'/mosaic_'+band+'_'+str(chunk_number)+'_ima'+name+'_4reb.fits'
		#imagefile=wd+'/deleteme3.fits'
		defnbkgfile=wd+'/mosaic_'+band+'_'+str(chunk_number)+'_defnbkg'+name+'_4reb.fits'
		#defnbkgfile=wd+'/deleteme4.fits'
		regfile=wd+'/mosaic_'+band+'_'+str(chunk_number)+'_src'+name+'_4reb.reg'
		#regfile=wd+'/mosaic_'+band+'_4rebinned_src.reg'
		sigthresh='5e-5'
		clobber='yes'
	
		s.call('wavdetect infile='+infile+' bkginput='+bkginput+' psffile='+psffile+' outfile='+outfile+' scellfile='+scellfile+' imagefile='+imagefile+' defnbkgfile='+defnbkgfile+' regfile='+regfile+' scales='+scales+' sigthresh='+sigthresh+' clobber='+clobber+'',shell=True)

sys.exit()
##############################
##############################
for chunk_number in range(0,6):
	s.call('punlearn wavdetect', shell=True)
	infile=wd+str(chunk_number)+'_merged_evt.fits'
	bkginput=''
	psffile=''
	outfile=wd+'/mosaic_'+str(chunk_number)+'_src.fits'
	scellfile=wd+'/mosaic_'+str(chunk_number)+'_scell.fits'
	imagefile=wd+'/mosaic_'+str(chunk_number)+'_ima.fits'
	defnbkgfile=wd+'/mosaic_'+str(chunk_number)+'_defnbkg.fits'
	regfile=wd+'/mosaic_'+str(chunk_number)+'_src.reg'
	scales='\"1.0 1.4142 2.0 4.0\"'
	sigthresh='1e-6'
	clobber='yes'
	
	s.call('wavdetect infile='+infile+' bkginput='+bkginput+' psffile='+psffile+' outfile='+outfile+' scellfile='+scellfile+' imagefile='+imagefile+' defnbkgfile='+defnbkgfile+' regfile='+regfile+' scales='+scales+' sigthresh='+sigthresh+' clobber='+clobber+'',shell=True)

sys.exit()

# WAVDETECT ON NEW OBSID
#(obsid)=np.genfromtxt('/Users/alberto/Desktop/XBOOTES/data_counts.dat',unpack=True, usecols=[1],dtype='str')
obsid=['19676','19777','19778','19779']

for i in range(len(obsid)):

	if len(obsid[i]) == 4:
		stem='0'+obsid[i]
	elif len(obsid[i]) == 3:
		stem='00'+obsid[i]
	elif len(obsid[i]) == 5:
		stem=obsid[i]
	'''
	s.call('punlearn wavdetect', shell=True)
	infile=wd+obsid[i]+'/repro/acisf'+stem+'_repro_evt2.fits'
	psffile=''
	outfile=wd+'/junk_output/'+stem+'_src.fits'
	scellfile=wd+'/junk_output/'+stem+'_scell.fits'
	imagefile=wd+'/junk_output/'+stem+'_ima.fits'
	defnbkgfile=wd+'/junk_output/'+stem+'_defnbkg.fits'
	regfile=wd+'/junk_output/'+stem+'_src.reg'
	scales='\"1.0 2.0 4.0\"'
	#sigthresh='1e-6'
	falsesrc='3.0' #number of false source detections allowed per scale per field
	clobber='yes'
	'''
	s.call('punlearn wavdetect', shell=True)
	infile=wd+obsid[i]+'/repro_new_asol/acisf'+stem+'_repro_evt2.fits'
	#bkginput=wd+obsid[i]+'/repro_new_asol/out/acisf'+stem+'_bkgmap_ximage.fits',shell=True)
	psffile=''
	outfile=wd+'/junk_output/'+stem+'_src_new_asol.fits'
	scellfile=wd+'/junk_output/'+stem+'_scell_new_asol.fits'
	imagefile=wd+'/junk_output/'+stem+'_ima_new_asol.fits'
	defnbkgfile=wd+'/junk_output/'+stem+'_defnbkg_new_asol.fits'
	regfile=wd+'/junk_output/'+stem+'_src_new_asol.reg'
	scales='\"1.0 2.0 4.0\"'
	#sigthresh='5e-5'
	falsesrc='3.0' #number of false source detections allowed per scale per field
	clobber='yes'
	
	s.call('wavdetect infile='+infile+' psffile='+psffile+' outfile='+outfile+' scellfile='+scellfile+' imagefile='+imagefile+' defnbkgfile='+defnbkgfile+' regfile='+regfile+' scales='+scales+' clobber='+clobber+'',shell=True)
#INPUT
#s.call('pset wavdetect infile='+wd+'/sim_full/acisf'+stem+'_sim_poiss_ximage.fits',shell=True)
#s.call('pset wavdetect infile='+wd+dirs+'/'+obsid+'/repro/acisf'+stem+'_repro_05to7keV_bin1_ximage.img',shell=True)
#s.call('pset wavdetect bkginput='+wd+dirs+'/'+obsid+'/repro/out2/acisf'+stem+'_bkgmap_ximage.fits',shell=True)
#s.call('pset wavdetect expfile='+wd+dirs+'/'+obsid+'/repro/out2/acisf'+stem+'_expomap_ximage.fits',shell=True)
#s.call('pset wavdetect psffile='+wd+dirs+'/'+obsid+'/repro/out2/acisf'+stem+'_psfmap_ximage.fits',shell=True)

#OUTPUT
#s.call('pset wavdetect outfile='+wd+'/sim_full/acisf'+stem+'_sim_src.fits',shell=True)
#s.call('pset wavdetect scellfile='+wd+'/sim_full/acisf'+stem+'_sim_scell.fits',shell=True)
#s.call('pset wavdetect imagefile='+wd+'/sim_full/acisf'+stem+'_sim_ima.fits',shell=True)
#s.call('pset wavdetect defnbkgfile='+wd+'/sim_full/acisf'+stem+'_sim_defnbkg.fits',shell=True)
#s.call('pset wavdetect regfile='+wd+'/sim_full/acisf'+stem+'_sim_src.reg',shell=True)

#s.call('pset wavdetect log=yes',shell=True)
#s.call('pset wavdetect verbose=5',shell=True)
#s.call('pset wavdetect maxiter=4',shell=True)

PFILES = $(PWD)/pfiles;../
HP_DATA=hpdata
export HP_DATA

process: 
#test
	python ../../run/process.py 2095 hpdata
# hpdata/hp_5/hp_123/hp_209*.fits

test:
	mkdir -pv pfiles

	plist varmosaic


	../varmosaic filelist=image_list.txt outimage=mosaic.fits outregion=mosaic.reg coordtype=-TAN pixdivide=5



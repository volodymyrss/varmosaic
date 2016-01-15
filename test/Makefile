PFILES = $(PWD)/pfiles;../
HP_DATA=hpdata
export HP_DATA

test:
	mkdir -pv pfiles

	plist varmosaic


	../varmosaic filelist=image_list.txt outimage=mosaic.fits outregion=mosaic.reg coordtype=-TAN pixdivide=3

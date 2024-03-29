#include <stdio.h>
#include <string.h>
#include <math.h>
#include "pil.h"
#include "headas.h"        
#include "chealpix.h"
#include <sys/stat.h>

                                                                                
int printsub(const char *);
                                                                                
#define TOOLSUB varmosaic
/* headas_main() requires that TOOLSUB be defined first */
#include "headas_main.c"

#define VARMOSAIC_VERSION "varmosaic  ver 3.0 (2005-03-07)"

#define  MAXFILE 1000      /* Maximum number of input files */
#define  MAXNAMELENGTH 256  /* Maximum character length of file names */
#define  BufLen_2      MAXNAMELENGTH /* Required for pfile.h */
#define  MAX_E_RANGES 20    /* Maximum number of energy ranges for the input images */
    
#define offset_first_output_image 2
#define output_images_per_band 4

/* Carry out quick mosaic of integral images taking account of variance 
   2004-11-10 version 1.0 Ken Ebisawa 
   2004-11-11 version 1.1 Ken Ebisawa 
   In ISDC Barn for INTEGRAL Quick Look Analysis, somehow fits_open_image does not work.
   If there is a blank line in the input file, it was considered
   to be a file with null name.

   2005-03-07 version 2.0. First release an ftool from GSFC
*/


/* Function prototypes */
//int varmosaic(void);
int write_wcs(fitsfile *outfptr,double xrval,double yrval,double xrpix,double yrpix,double xinc,double yinc,double rot, char ctype[], char ctype2[], int *status);
void writebintable ( char fn[], int n_scw, double * ijds, long int * i_hp1024, double *flux, double *variance, double *exposure, int *nsubpx );

int varmosaic(void){
  int status = 0, i, j, nfile, hdutype,jj, l, m;
  fitsfile *infptr,*outfptr;
  FILE *fp,*region_file;
  float total_exposure;
  char filelist[MAXNAMELENGTH];
  char inimages[MAXFILE][MAXNAMELENGTH], tmp[MAXNAMELENGTH], coordtype[10];
  int inimages_mode[MAXFILE];
  double xrefval[MAXFILE], yrefval[MAXFILE], xrefpix, yrefpix, xinc, yinc, rot;
  double pointx[MAXFILE],pointy[MAXFILE],pointz[MAXFILE],ptstart[MAXFILE];
  double avex=0.0,avey=0.0,avez=0.0,deg2rad=1.745329252E-2;
  double ngood=0;
  double tstart;
  char instrume[20], comment[MAXNAMELENGTH],outfile[MAXNAMELENGTH], outregion[MAXNAMELENGTH], out_coordtype[10], ctype1[20],  ctype2[20];
  char text[9999]; /* Output text messages */

  double maxoffset=60;
  
  double x_ref=0,y_ref=0,z_ref=0;
  double aveRA, aveDEC, xrefpixout, yrefpixout;   /* WCS information of the output file. Note, xinc and yinc


						     are the same as those of input files  */
  int imagebin;    /* Output image bins (X and Y same)  */
  int num_E_range=-1; /* Number of energy ranges in the input images/output image */

  
  int clobber; /* clobber the output image or not */
  int use_ref=0; /* clobber the output image or not */
  int exist;  /*output file alreayd exists or not*/
  int pixdivide; /* input image pixel is divided into pixdivide**2 subpixels */

  long inimsize,inimsize1,inimsize2; /* Input imagesize. Readfrom NAXIS1 */
    double offset;

  double dist=0.0, maxdist=-999.0,imagesize=0;
    
    long int i_hp1;
    long int i_hp4;
    long int i_hp16;
    long int i_hp1024;
    long int i_hp1024_in_hp16;
    
    double theta,phi;
    double ra,dec;


  struct ScWimage{ /* Defined for each energy range               */
    int    intHDU; /* HDU of the intensity map                    */
    int    varHDU; /* HDU of the variance map                     */
    int    expoHDU; /* HDU of the variance map                     */
    double E_min;  /* E_min of the energy range                   */
    double E_max;  /* E_max of the energy range                   */
  } ScW[MAX_E_RANGES];

  static char taskname[80] = "varmosaic";
  static char version[8] = "3.0";
  /* Register taskname and version. */
  set_toolname(taskname);
  set_toolversion(version);

  /* Read parameters from the paramete file */
  status=PILGetString("filelist", filelist);

  /* Read ascii image list*/
  if((fp = fopen(filelist,"r"))==NULL){
    sprintf(text,"##### Input file \"%s\" NOT opened.\n",filelist);
    HD_printf(text);
    return status;
  }else{
    sprintf(text,"##### Input file \"%s\" opened\n",filelist);
    HD_printf(text);
    i = 0;
    inimages_mode[i]=0;
    while(fgets(tmp, MAXNAMELENGTH,fp)!=NULL&&sscanf(tmp,"%s",inimages[i])>0){
      i++;
      inimages_mode[i]=0;
      sprintf(text,"file in the input list: %s %d mode %i\n", inimages[i],nfile,inimages_mode[i]);  
      HD_printf(text);
    }
    fclose(fp);
  }
  nfile = i;
  sprintf(text,"##### Total number of input files: %d\n", nfile);  
  HD_printf(text);

  /* Read other parameters */
  status=PILGetString("outregion", outregion);
  status=PILGetString("outimage", outfile);
  status=PILGetString("coordtype", out_coordtype);

  sprintf(ctype1,"RA--%s",out_coordtype);
  sprintf(ctype2,"DEC-%s",out_coordtype);

  sprintf(text,"output mosaic in %s and region in %s\n", outfile, outregion);  
  HD_printf(text);
  status=PILGetInt("pixdivide", &pixdivide);
  sprintf(text,"##### Input pixel is divided into %d x %d sub-pixels\n", pixdivide, pixdivide);  
  HD_printf(text);

  status=PILGetBool("clobber",&clobber);

  {FILE *openfile; /* Check if the output file exists or not */
  openfile = fopen(outfile,"r");
  if(openfile != NULL) { 
    fclose(openfile);
    exist = 1;
  }else    {
    exist = 0;
  }}

  if(exist){
    if(clobber){
      sprintf(text,"##### Output file exists and will be overwritten since clobber=yes");
      headas_clobberfile(outfile);
    }
    else{
      sprintf(text,"##### Output file exists and will NOT be overwritten since clobber=no");
      fprintf(stderr,text);
      return status;
    }
  }
      
  status=PILGetReal("maxoffset",&maxoffset);
  sprintf(text,"##### Maximum allowed offset %.5lg",maxoffset);

  status=PILGetInt("outimagesize",&imagebin);
  if(imagebin > 0){
    /* Output image size and reference points read from the parameter file */
      sprintf(text,"##### Output image size is given in the parameter file as %d\n", imagebin);
      HD_printf(text);
      status=PILGetReal("ra_ref",&aveRA);
      HD_printf(text);
      status=PILGetReal("dec_ref",&aveDEC);
      sprintf(text,"##### Reference point RA, DEC = %10.3f %10.3f\n", aveRA, aveDEC);
      HD_printf(text);
      x_ref=cos(aveRA*deg2rad)*cos(aveDEC*deg2rad);
      y_ref=sin(aveRA*deg2rad)*cos(aveDEC*deg2rad);
      z_ref=sin(aveDEC*deg2rad);
      use_ref=1;
  }

  
  HD_printf("");

  /* Read the input image FITS file one by one, determine the average of the pointings, which 
     will be the center of the mosaic image.  Also, we deterime the mosaic image extent,
     to determine the size of the mosaic image.  Note, we use the fixed  imagebin size, 
     same as the ScW level images. Fom the first input file, we get the number of energy bands. 
     Assumption is that all the input files have the same number of energy bands */

  region_file=fopen(outregion,"w");
  
  fprintf(region_file,"# comments\n");
  fprintf(region_file,"global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n");
  fprintf(region_file,"fk5\n");
  
  for(i=0;i<nfile;i++){ /* loop DDD for each input file*/
    /*workaround in Barn. Somehow fits_open_image detect the empty primary
      fits_open_image(&infptr, inimages[i], 0, &status);
      {int hdunum;
      fits_get_hdu_num(infptr, &hdunum);
      }
    */
    status=0;

    fits_open_data(&infptr, inimages[i], 0, &status);
    sprintf(text,"Input File[%d]: %s\n", i+1,inimages[i]);
    HD_printf(text);

    if (status!=0) {
        sprintf(text,"failed to open, skipping\n");
        inimages_mode[i]=1;
        HD_printf(text);
        status=0;
        continue;
    };


    fits_movabs_hdu(infptr, 3, &hdutype, &status);

    fits_read_key_dbl(infptr, "TSTART", &tstart,comment,&status);
    fits_read_key_str(infptr, "INSTRUME", instrume,comment,&status);
    fits_read_img_coord(infptr, &xrefval[i], &yrefval[i], &xrefpix, &yrefpix, &xinc,  &yinc,
			&rot, coordtype, &status);
    pointx[i]=cos(xrefval[i]*deg2rad)*cos(yrefval[i]*deg2rad);
    pointy[i]=sin(xrefval[i]*deg2rad)*cos(yrefval[i]*deg2rad);
    pointz[i]=sin(yrefval[i]*deg2rad);
    ptstart[i]=tstart;

    if (use_ref) {
        offset=acos((pointx[i]*x_ref+pointy[i]*y_ref+pointz[i]*z_ref))/deg2rad;
    } else {
        if (ngood>0) {
            offset=acos((pointx[i]*avex+pointy[i]*avey+pointz[i]*avez)/pow(avex*avex+avey*avey+avez*avez,0.5))/deg2rad;
        } else offset=0;
    };
    
    sprintf(text,"Pointing %.15lg direction: (%7.3f,%7.3f) rotation %.5g, offset from average  (%.5g,%.5g) %.5g\n",tstart, xrefval[i], yrefval[i],rot,avex,avey,avez,offset);
    HD_printf(text);

    if (offset>180.) {
        inimages_mode[i]=2;
        sprintf(text,"offset too large (%lg), skipping\n", offset);
        HD_printf(text);
        fits_close_file(infptr, &status);
        continue;
    };

    avex= avex + pointx[i];
    avey= avey + pointy[i];
    avez= avez + pointz[i];
    ngood+=1;

    
    
    if(num_E_range<=0){/* block CCC */
      /* Look into the first image file structure to know the energy bands*/
      {int hdunum,ii=0;
      char imatype[20];
      double E_min, E_max;
      fits_get_num_hdus(infptr, &hdunum, &status);
      sprintf(text,"Looking into the structure of the first image. Number of HDU = %d\n", hdunum);
      HD_printf(text);
      for(j=1;j<=hdunum;j++){
	fits_movabs_hdu(infptr, j, &hdutype, &status);
	fits_read_key_lng(infptr, "NAXIS1",   &inimsize, comment, &status);
	fits_read_key_str(infptr, "IMATYPE", imatype,comment, &status);
	fits_read_key_dbl(infptr, "E_MIN",   &E_min, comment, &status);
	fits_read_key_dbl(infptr, "E_MAX",   &E_max, comment, &status);
	if(strstr(imatype,"INTENSITY") != NULL){
	  ScW[ii].intHDU = j;
	  ScW[ii].E_min = E_min;
	  ScW[ii].E_max = E_max;
	}else if(strstr(imatype,"VARIANCE") != NULL) {
	  ScW[ii].varHDU = j;
	  //ii++;
	}else if(strstr(imatype,"EXPOSURE") != NULL && strstr(imatype,"EXPOSUREFLAT") == NULL) {
	  ScW[ii].expoHDU = j;
	  ii++;
	}
	if(status !=0){
	  status = 0;
	}else{
	  sprintf(text,"HDU = %2d IMATYPE = %12s E_min = %7.2f keV E_max = %7.2f keV\n", j, imatype, E_min, E_max);
	  HD_printf(text);
}
      }
      num_E_range = ii;
      sprintf(text,"Number of energy ranges = %d\n", num_E_range);
      HD_printf(text);
      for(ii=0;ii<num_E_range;ii++){
	sprintf(text,"Energy range:%d Intensity HDU:%d Variance HDU:%d Exposure HDU:%d E_min %4.2f keV E_max %4.2f keV\n", ii+1, ScW[ii].intHDU,ScW[ii].varHDU,ScW[ii].expoHDU,ScW[ii].E_min,ScW[ii].E_max);
	HD_printf(text);
      }
      HD_printf("");
      }
    }/*End CCC for the first input file */
      fits_close_file(infptr, &status);

     fprintf(region_file,"vector(%.5g,%.5g,0.5,%.5g)\n",xrefval[i],yrefval[i],rot);
  } /* End DDD for each input file */
  
  fclose(region_file);

  if(imagebin <=0){ /* block LLL when imagebin in the parameter file <=0,
		       Determine the average poiting position and output image size from the input files*/
    { 
      double length=0.0;
      length = sqrt(avex*avex+avey*avey+avez*avez);
      avex=avex/length;
      avey=avey/length;
      avez=avez/length;
    }
    
    aveRA = atan2(avey,avex)/deg2rad;
    if(aveRA<0){aveRA = aveRA + 360.0;}
    aveDEC= asin(avez)/deg2rad;
    
    sprintf(text,"Average Pointing direction is calculated as: (%7.3f,%7.3f)\n", aveRA, aveDEC);
    HD_printf(text);
    
    /* Determine the maximum distance between input image centers and the pointing average. 
       Create the output image mosaic file having the same number of energy bands as
       input imagess,  write necessary keywords to the output file, in particular WCS keywords.*/
    
    if (nfile>1) {
        for(i=0;i<nfile;i++){
            if (inimages_mode[i]!=0) {
                continue;
            };
          dist = sqrt((avex-pointx[i])*(avex-pointx[i])+(avey-pointy[i])*(avey-pointy[i])+(avez-pointz[i])*(avez-pointz[i]));
          if(dist > maxdist){maxdist = dist;}
        }
        maxdist = 2.0*asin(maxdist/2.0)/deg2rad;
    } else maxdist=0;
    sprintf(text,"\nMaximum distance of the mosaic center and each ScW pointing :%8.3f degree\n", maxdist);
    HD_printf(text);
    /* ISGRI FOV = 30.6 deg x 30.6 deg which is covered with a circle with 21.6 deg radius.
       JEMX  FOV = 6.6 deg radius */
    if(strchr(instrume,'J')!=NULL){imagesize = 2.0*(maxdist +  6.6);}
    if(strchr(instrume,'I')!=NULL){imagesize = 2.0*(maxdist + 21.6);}
    imagebin =  floor(imagesize/yinc)+1;
    sprintf(text,"We define the output image size as %9.3f deg.\n", imagesize);
    HD_printf(text);
    sprintf(text,"Imagebinsize = %9.4f deg (same as input), thus output image size = %d pixels\n", yinc,imagebin);
    HD_printf(text);
  } /* end of block LLL */
  
  {/* start block FFF */
    long naxes[2];
    char extname[20];
    
    fits_create_file(&outfptr, outfile, &status);
    naxes[0] = imagebin;naxes[1]=imagebin;
    xrefpixout = (1.0+ imagebin)/2.0;
    yrefpixout = (1.0+ imagebin)/2.0;
    fits_create_img(outfptr, FLOAT_IMG, 2, naxes, &status);
    fits_create_img(outfptr, FLOAT_IMG, 2, naxes, &status); // offset here
    for(jj=1;jj<=num_E_range*4+1;jj++){ /* start loop EEE for each output HDU*/
      fits_create_img(outfptr, FLOAT_IMG, 2, naxes, &status);
      fits_write_key_str(outfptr, "INSTRUME", instrume, "Instrument name", &status);
      fits_write_key_str(outfptr, "CREATOR", VARMOSAIC_VERSION, "Program to create this file", &status);
      if(jj>=num_E_range*4+1){
        if(jj==num_E_range*4+1){
            fits_write_key_str(outfptr, "IMATYPE", "EXPOSUREFLAT", "Type of image", &status);
            fits_write_key_str(outfptr, "EXTNAME", "EXPOSUREFLAT", "Extension name", &status);
            fits_write_key_str(outfptr, "BUNIT", "sec", "Unit for pixel values", &status);
        };
      }else{
	fits_write_key_dbl(outfptr,"E_MIN",ScW[(jj-1)/4].E_min,-4,"Lower energy boundary [keV]",&status);
	fits_write_key_dbl(outfptr,"E_MAX",ScW[(jj-1)/4].E_max,-4,"Upper energy boundary [keV]",&status);
	if(jj%4==1){
	  fits_write_key_str(outfptr, "IMATYPE", "INTENSITY", "Type of image", &status);
	  sprintf(extname,"INTENSITY");
	  //sprintf(extname,"INTENSITY%-d",(jj-1)/4+1);
	  fits_write_key_str(outfptr, "EXTNAME", extname, "Extension name", &status);
	  fits_write_key_str(outfptr, "BUNIT", "count/s", "Unit for pixel values", &status);
	}else if(jj%4==2){
	  fits_write_key_str(outfptr, "IMATYPE", "VARIANCE", "Type of image", &status);
	  sprintf(extname,"VARIANCE");
	  fits_write_key_str(outfptr, "EXTNAME", extname, "Extension name", &status);
	  fits_write_key_str(outfptr, "BUNIT", "(count/s)**2", "Unit for pixel values", &status);
	}else if(jj%4==3){
	  fits_write_key_str(outfptr, "IMATYPE", "SIGNIFICANCE", "Type of image", &status);
	  sprintf(extname,"SIGNIFICANCE");
      fits_write_key_str(outfptr, "EXTNAME", extname, "Extension name", &status);
	  fits_write_key_str(outfptr, "BUNIT", "no unit", "Unit for pixel values", &status);
	}else if(jj%4==0){
	  fits_write_key_str(outfptr, "IMATYPE", "EXPOSURE", "Type of image", &status);
	  sprintf(extname,"EXPOSURE");
	  fits_write_key_str(outfptr, "EXTNAME", extname, "Extension name", &status);
	  fits_write_key_str(outfptr, "BUNIT", "seconds", "effective seconds", &status);
	}
      }
      write_wcs(outfptr,aveRA,aveDEC,xrefpixout,yrefpixout,xinc,yinc,0.0E0,ctype1,ctype2,&status);
      sprintf(comment, "Number of subpixel divisions: %2dx%2d",pixdivide,pixdivide);
      fits_write_comment(outfptr, comment, &status);
      sprintf(comment, "Total number of input files: %5d",nfile);
      fits_write_comment(outfptr, comment, &status);
      for(i=0;i<nfile;i++){
	sprintf(comment, "%4d: %s", i+1,inimages[i]);
	fits_write_comment(outfptr, comment, &status);
      }
    }/* end of EEE for each output HDU*/
    fits_close_file(outfptr, &status);    
  }/* end of block FFF. End creating the empty output FITS file. Necessary extensions and keywords are
      written to the output file.*/
  
  /* Allocate input and output image buffer.  We need output images, the number of which is
     3 times (flux, variance, significance) the number of energy ranges plus 1 (exposure) */
  {/* start block ZZZ*/
    int anynul,ii,kk,ienergy;
    float *in_flx_image[MAX_E_RANGES], *in_var_image[MAX_E_RANGES], *in_expo_image[MAX_E_RANGES];
    float *flx_image[MAX_E_RANGES], *var_image[MAX_E_RANGES], *sig_image[MAX_E_RANGES], *expo_image[MAX_E_RANGES], *exp_image;
    double ScWxrefval,ScWyrefval,ScWxrefpix,ScWyrefpix,ScWrot,xpos,ypos;
    double XpixOut,YpixOut;
    float flux, variance, exposure, exposure_flat;
    float sub_exposure, sub_variance, sub_exposure_flat;    
    int intXpixOut,intYpixOut;
    
    long int * hp1024_in_hp16=malloc(nside2npix(1024)*sizeof(long int));
    long int * inv_hp1024_in_hp16=malloc(nside2npix(1024)*sizeof(long int));
    long int aux_hp1024_in_hp16[16*16*12];

    double **flux_hp16[MAX_E_RANGES];
    double **var_hp16[MAX_E_RANGES];
    double **expo_hp16[MAX_E_RANGES];
    int **nsubpx_hp16[MAX_E_RANGES];
    double **ijd_hp16[MAX_E_RANGES];
    int *nscw_hp16[MAX_E_RANGES];

    for (i_hp16=0;i_hp16<nside2npix(16);i_hp16++) {
        aux_hp1024_in_hp16[i_hp16]=0;
    };

    for (i_hp1024=0;i_hp1024<nside2npix(1024);i_hp1024++) {
        pix2ang_ring(1024,i_hp1024,&theta,&phi);
        ang2pix_ring(16,theta,phi,&i_hp16);

        hp1024_in_hp16[i_hp1024]=aux_hp1024_in_hp16[i_hp16];
        aux_hp1024_in_hp16[i_hp16]+=1;
        inv_hp1024_in_hp16[i_hp16*4096+hp1024_in_hp16[i_hp1024]]=i_hp1024;
    };
    

    for(ienergy=0;ienergy<num_E_range;ienergy++){
      expo_image[ienergy]=malloc(imagebin*imagebin*sizeof(float));
      /*in_flx_image[ienergy]=malloc(inimsize*inimsize*sizeof(float));
      in_var_image[ienergy]=malloc(inimsize*inimsize*sizeof(float));    
      in_expo_image[ienergy]=malloc(inimsize*inimsize*sizeof(float));    */
      flx_image[ienergy]=malloc(imagebin*imagebin*sizeof(float));
      sig_image[ienergy]=malloc(imagebin*imagebin*sizeof(float));
      var_image[ienergy]=malloc(imagebin*imagebin*sizeof(float));
      expo_image[ienergy]=malloc(imagebin*imagebin*sizeof(float));

    // allocate healpix
      if((ijd_hp16[ienergy]=malloc(nside2npix(16)*sizeof(double*)))==NULL) {
            printf("unable to allocate!");
      };

      if((nscw_hp16[ienergy]=malloc(nside2npix(16)*sizeof(double)))==NULL) {
            printf("unable to allocate!");
      };
      if((flux_hp16[ienergy]=malloc(nside2npix(16)*sizeof(double*)))==NULL) {
            printf("unable to allocate!");
      };
      if((var_hp16[ienergy]=malloc(nside2npix(16)*sizeof(double*)))==NULL) {
            printf("unable to allocate!");
      };
      if((expo_hp16[ienergy]=malloc(nside2npix(16)*sizeof(double*)))==NULL) {
            printf("unable to allocate!");
      };
      if((nsubpx_hp16[ienergy]=malloc(nside2npix(16)*sizeof(int*)))==NULL) {
            printf("unable to allocate!");
      };

      for (i_hp16=0;i_hp16<nside2npix(16);i_hp16++) {
        nscw_hp16[ienergy][i_hp16]=0;
        ijd_hp16[ienergy][i_hp16]=malloc(MAXFILE*sizeof(double));
        flux_hp16[ienergy][i_hp16]=malloc(4096*sizeof(double));
        var_hp16[ienergy][i_hp16]=malloc(4096*sizeof(double));
        expo_hp16[ienergy][i_hp16]=malloc(4096*sizeof(double));
        nsubpx_hp16[ienergy][i_hp16]=malloc(4096*sizeof(int));
      };

    }
    exp_image=malloc(imagebin*imagebin*sizeof(float));
    
    /* initialize the output images */
    for(i=0;i<imagebin;i++){
      for(j=0;j<imagebin;j++){
	for(ienergy=0;ienergy<num_E_range;ienergy++){
	  flx_image[ienergy][i+imagebin*j]=0.0;
	  sig_image[ienergy][i+imagebin*j]=0.0;
	  var_image[ienergy][i+imagebin*j]=0.0;
	  expo_image[ienergy][i+imagebin*j]=0.0;
	} 
      }
    }
    
    /* Open the input images one by one, read soft flux, soft variance, hard flux and
       hard variance images.  Add them properly taking account of variance weight,
       and write into the output file.
       Note, we do not divide pixels to save time; this causes quantization error
       below pixel level.*/

    total_exposure=0;
    for(i=0;i<nfile;i++){/* start loop AA for each input file*/
        sprintf(text,"Image %i bands %i\n", i, num_E_range);
        HD_printf(text);
        if (inimages_mode[i]!=0) {
            sprintf(text,"image disabled, code %i\n", inimages_mode[0]);
            HD_printf(text);
            continue;
        };

      fits_open_image(&infptr, inimages[i], 0, &status);


      for(ienergy=1;ienergy<=num_E_range;ienergy++){/* Read all the flux and variance maps in the input file */

          fits_movabs_hdu(infptr,ScW[ienergy-1].intHDU, &hdutype, &status); /* flux HDU */
          // allocate for each
          fits_read_key_lng(infptr, "NAXIS1",   &inimsize1, comment, &status);
          fits_read_key_lng(infptr, "NAXIS2",   &inimsize2, comment, &status);

          if (status!=0) {
              sprintf(text,"failed to read image dimentions!\n");
              HD_printf(text);
              break;
          }; 

          printf("image %i energy band %i size %li x %li, allocating\n",i,ienergy,inimsize1,inimsize2);
          in_flx_image[ienergy-1]=malloc(inimsize1*inimsize2*sizeof(float));
          in_var_image[ienergy-1]=malloc(inimsize1*inimsize2*sizeof(float));    
          in_expo_image[ienergy-1]=malloc(inimsize1*inimsize2*sizeof(float));    

          fits_read_img_coord(infptr, &ScWxrefval, &ScWyrefval, &ScWxrefpix, &ScWyrefpix,&xinc,&yinc,
                  &ScWrot, coordtype, &status);

          fits_read_key_flt(infptr, "EXPOSURE",   &exposure_flat, comment, &status);
          if(strchr(instrume,'I')&&exposure_flat<=0.0){
              /* Occationally, ISGRI has a bug that EXPOSURE = 0.0 */

              fits_read_key_flt(infptr, "TELAPSE",   &exposure_flat, comment, &status);
          };

          if (ienergy==1)
            total_exposure+=exposure_flat;

          printf("Exposure: %.5g coordtype %s; status %i\n",exposure_flat,coordtype,status);

          fits_read_2d_flt(infptr, 0L, 0, inimsize1, inimsize1, inimsize2, in_flx_image[ienergy-1], &anynul, &status);
          fits_movabs_hdu(infptr,ScW[ienergy-1].varHDU, &hdutype, &status); /* variance HDU */
          fits_read_2d_flt(infptr, 0L, 0, inimsize1, inimsize1, inimsize2, in_var_image[ienergy-1], &anynul, &status);
          fits_movabs_hdu(infptr,ScW[ienergy-1].expoHDU, &hdutype, &status); /* variance HDU */
          fits_read_2d_flt(infptr, 0L, 0, inimsize1, inimsize1, inimsize2, in_expo_image[ienergy-1], &anynul, &status);
      }
      if (status!=0) {
          sprintf(text,"failed to read image dimentions!\n");
          HD_printf(text);
          status=0;
          break;
      }; 

      fits_close_file(infptr, &status);
      sprintf(text,"Reading image:%5d %s, total exposure %.5g ks, status %i\n", i, inimages[i],total_exposure/1e3,status);
      HD_printf(text);
      if (status!=0) {
          sprintf(text,"failed to open!\n");
          HD_printf(text);
          status=0;
          continue;
      }; 

      for(ii=1;ii<=inimsize1;ii++){for(jj=1;jj<=inimsize2;jj++){/*start loop YYY for input image pixels */
	  /* Now we are dividing this input image pixel into sub-pixels.*/

        long int subpixels_filled=0;

	for(l=1; l<=pixdivide; l++){for(m=1; m<=pixdivide; m++){ /* Start loop CCC repeat for sub-pixels */
	  for(ienergy=1;ienergy<=num_E_range;ienergy++){/*start loop XXX for energy bands */
	    /* Read all the flux and variance maps in the input file 
	       Here, ii and jj are the X and Y pixel number of the input images,
	       starting from 1, not zero.
	    */
	    flux     = in_flx_image[ienergy-1][(ii-1)+(jj-1)*inimsize1]; // guess
	    variance = in_var_image[ienergy-1][(ii-1)+(jj-1)*inimsize1];
	    exposure = in_expo_image[ienergy-1][(ii-1)+(jj-1)*inimsize1];
	    if(flux==flux&&variance==variance&&variance>0.0){ /* not NULL */
	      {
		double subpixelX, subpixelY, subpixelsize;
		subpixelsize = 1.0/ (double) pixdivide;
		subpixelX = (double) ii - 0.5 + (l-0.5) * subpixelsize;
		subpixelY = (double) jj - 0.5 + (m-0.5) * subpixelsize;
        if (status!=0) {
            printf("failed before fits_pix_to_world\n");
            exit(1);
        };
		fits_pix_to_world(subpixelX,subpixelY,ScWxrefval,ScWyrefval,ScWxrefpix,ScWyrefpix,xinc,yinc,ScWrot,coordtype,&xpos,&ypos,
				  &status);
        if (status!=0) {
            printf("failed to fits_pix_to_world status %i\n",status);
            exit(1);
        };
	      }
		
        /* Each sub-pixel carries the associated values,
		   sub_exposure=exposure/sub_pixel_numbers, sub_pixel_flux=flux, and sub_variance=sub_pixel_numbers*variance.*/
		sub_exposure = exposure; /// (float)(pow(pixdivide,2));
		sub_exposure_flat = exposure_flat; // (float)(pow(pixdivide,2));
		sub_variance = variance*(float) pow(pixdivide,2);

        //  <healpix
        

        double theta=(ypos+90.)/180.*M_PI;
        double phi=xpos/180.*M_PI;
        long ipix_1;
        long ipix_4;
        long ipix_16;
        long ipix_1024;
        
        ang2pix_ring(1,theta,phi,&ipix_1);
        ang2pix_ring(4,theta,phi,&ipix_4);
        ang2pix_ring(16,theta,phi,&ipix_16);
        ang2pix_ring(1024,theta,phi,&ipix_1024);

        /*if (ipix_1024==8895424) {
            printf("hp1024 pixel %li\n",ipix_1024);
        };*/
            

        int scw_in_pix=nscw_hp16[ienergy-1][ipix_16];
        int i_s;


        if (scw_in_pix==0 || ijd_hp16[ienergy-1][ipix_16][scw_in_pix-1]!=ptstart[i]) {
            scw_in_pix=++nscw_hp16[ienergy-1][ipix_16];

            ijd_hp16[ienergy-1][ipix_16][scw_in_pix-1]=ptstart[i];

           // flux_hp16[ienergy-1][ipix_16][scw_in_pix-1]=malloc((nside2npix(1024)/(nside2npix(16))*sizeof(double)));
           // var_hp16[ienergy-1][ipix_16][scw_in_pix-1]=malloc((nside2npix(1024)/(nside2npix(16))*sizeof(double)));
                
            flux_hp16[ienergy-1][ipix_16]=realloc(flux_hp16[ienergy-1][ipix_16],(nside2npix(1024)/(nside2npix(16))*sizeof(double)*scw_in_pix));
            var_hp16[ienergy-1][ipix_16]=realloc(var_hp16[ienergy-1][ipix_16],(nside2npix(1024)/(nside2npix(16))*sizeof(double)*scw_in_pix));
            nsubpx_hp16[ienergy-1][ipix_16]=realloc(nsubpx_hp16[ienergy-1][ipix_16],(nside2npix(1024)/(nside2npix(16))*sizeof(double)*scw_in_pix));
            expo_hp16[ienergy-1][ipix_16]=realloc(expo_hp16[ienergy-1][ipix_16],(nside2npix(1024)/(nside2npix(16))*sizeof(double)*scw_in_pix));

           // printf("new hp16 pixel ijd %.15g %.5g %.5g ipix 1,4,16,1024 %li %li %li %li row %i\n",ptstart[i],xpos,ypos,ipix_1,ipix_4,ipix_16,ipix_1024,nscw_hp16[ienergy-1][ipix_16]);
            for (i_s=0;i_s<4096;i_s++) {
                flux_hp16[ienergy-1][ipix_16][(scw_in_pix-1)*4096+i_s]=0;
                var_hp16[ienergy-1][ipix_16][(scw_in_pix-1)*4096+i_s]=0;
                nsubpx_hp16[ienergy-1][ipix_16][(scw_in_pix-1)*4096+i_s]=0;
                expo_hp16[ienergy-1][ipix_16][(scw_in_pix-1)*4096+i_s]=0;
            };
        };


        //printf("%.li %.li %.li\n",ipix_16,ipix_1024,hp1024_in_hp16[ipix_1024]);
        //flux_hp16[ienergy-1][ipix_16][scw_in_pix-1][0]=flux;
        int ns=nsubpx_hp16[ienergy-1][ipix_16][(scw_in_pix-1)*4096+hp1024_in_hp16[ipix_1024]];

        flux_hp16[ienergy-1][ipix_16][(scw_in_pix-1)*4096+hp1024_in_hp16[ipix_1024]]*=ns;
        flux_hp16[ienergy-1][ipix_16][(scw_in_pix-1)*4096+hp1024_in_hp16[ipix_1024]]+=flux;
        flux_hp16[ienergy-1][ipix_16][(scw_in_pix-1)*4096+hp1024_in_hp16[ipix_1024]]/=(ns+1);
        

/*
        if (var_hp16[ienergy-1][ipix_16][(scw_in_pix-1)*4096+hp1024_in_hp16[ipix_1024]]>0)
            var_hp16[ienergy-1][ipix_16][(scw_in_pix-1)*4096+hp1024_in_hp16[ipix_1024]]=1./var_hp16[ienergy-1][ipix_16][(scw_in_pix-1)*4096+hp1024_in_hp16[ipix_1024]];
        var_hp16[ienergy-1][ipix_16][(scw_in_pix-1)*4096+hp1024_in_hp16[ipix_1024]]+=1./sub_variance;
        var_hp16[ienergy-1][ipix_16][(scw_in_pix-1)*4096+hp1024_in_hp16[ipix_1024]]=1./var_hp16[ienergy-1][ipix_16][(scw_in_pix-1)*4096+hp1024_in_hp16[ipix_1024]];
*/
        var_hp16[ienergy-1][ipix_16][(scw_in_pix-1)*4096+hp1024_in_hp16[ipix_1024]]*=ns;
        var_hp16[ienergy-1][ipix_16][(scw_in_pix-1)*4096+hp1024_in_hp16[ipix_1024]]+=variance;
        var_hp16[ienergy-1][ipix_16][(scw_in_pix-1)*4096+hp1024_in_hp16[ipix_1024]]/=(ns+1);
        
        expo_hp16[ienergy-1][ipix_16][(scw_in_pix-1)*4096+hp1024_in_hp16[ipix_1024]]*=ns;
        expo_hp16[ienergy-1][ipix_16][(scw_in_pix-1)*4096+hp1024_in_hp16[ipix_1024]]+=exposure;
        expo_hp16[ienergy-1][ipix_16][(scw_in_pix-1)*4096+hp1024_in_hp16[ipix_1024]]/=(ns+1);

        nsubpx_hp16[ienergy-1][ipix_16][(scw_in_pix-1)*4096+hp1024_in_hp16[ipix_1024]]=ns+1;


        subpixels_filled+=1;

        if (flux/pow(sub_variance,0.5)>5.) {
            printf("%.5lg %.5lg ipix 1,4,16,1024 %li %li %li %li flux %.5lg var %.5lg sig %.5lg\n",xpos,ypos,ipix_1,ipix_4,ipix_16,ipix_1024,flux,sub_variance,flux/pow(sub_variance,0.5));
            printf("hp flux %.5lg var %.5lg ns %i\n",flux_hp16[ienergy-1][ipix_16][(scw_in_pix-1)*4096+hp1024_in_hp16[ipix_1024]],var_hp16[ienergy-1][ipix_16][(scw_in_pix-1)*4096+hp1024_in_hp16[ipix_1024]],ns);
        };

        //  healpix>
	      
	      fits_world_to_pix(xpos,ypos,aveRA,aveDEC,xrefpixout,yrefpixout,xinc,yinc,0.0E0,out_coordtype,&XpixOut,&YpixOut,&status);
        if (status!=0) {
            printf("failed to fits_world_to_pix status %i\n",status);
            exit(1);
        };
	      if(XpixOut>=0.5&&XpixOut<=imagebin+0.5&&YpixOut>=0.5&&YpixOut<=imagebin+0.5){
		intXpixOut=floor(XpixOut+0.5);
		intYpixOut=floor(YpixOut+0.5);
		/* kk is the suffix for the output image pixels */
		kk = intXpixOut-1 +(intYpixOut-1)*imagebin;
        if(ienergy==1){exp_image[kk]=exp_image[kk]+sub_exposure_flat;}
        expo_image[ienergy-1][kk]=expo_image[ienergy-1][kk]+sub_exposure;
		if(var_image[ienergy-1][kk]>0.0){
		  flx_image[ienergy-1][kk]=flx_image[ienergy-1][kk]/var_image[ienergy-1][kk]+ flux/sub_variance;
		  {float tmp;
		  tmp = 1.0/var_image[ienergy-1][kk]+1.0/sub_variance;
		  flx_image[ienergy-1][kk]=flx_image[ienergy-1][kk]/tmp;
		  var_image[ienergy-1][kk]=1.0/tmp;
		  }
		}else{
		  flx_image[ienergy-1][kk]=flux;
		  var_image[ienergy-1][kk]=sub_variance;
		}
		sig_image[ienergy-1][kk]=flx_image[ienergy-1][kk]/sqrt(var_image[ienergy-1][kk]);
	      }
	    }
	  } /* end XXX, loop for energy bin for each pixel. Note, energy loop is inside pixel loop, assuming all the enegy bands have the same WCS information */
	}}/* end CCC repeating for the number of sub-pixels */

       // printf("subpixels used: %li\n",subpixels_filled);

      }}  /* end YYY, loop for each pixel of an input image*/
    }     /* end AA, loop for input images */

    
    /* Now we are going to write the output image */
    fits_open_image(&outfptr, outfile, 1, &status);
    for(ienergy=1;ienergy<=num_E_range;ienergy++){ /* loop for enery ranges */
      fits_movabs_hdu(outfptr,4*(ienergy-1)+1+offset_first_output_image,&hdutype, &status); /*significane */
      fits_write_key_dbl(outfptr, "EXPOSURE",  (double)total_exposure, -4, "exposure", &status);
      fits_write_2d_flt(outfptr, 0L, (long) imagebin, (long) imagebin,(long)imagebin,flx_image[ienergy-1],&status);
      fits_movabs_hdu(outfptr,4*(ienergy-1)+2+offset_first_output_image,&hdutype, &status); /*flux        */
      fits_write_2d_flt(outfptr, 0L, (long) imagebin, (long) imagebin,(long)imagebin,var_image[ienergy-1],&status);
      fits_movabs_hdu(outfptr,4*(ienergy-1)+3+offset_first_output_image,&hdutype, &status); /*variance    */
      fits_write_2d_flt(outfptr, 0L, (long) imagebin, (long) imagebin,(long)imagebin,sig_image[ienergy-1],&status);
      fits_movabs_hdu(outfptr,4*(ienergy-1)+4+offset_first_output_image,&hdutype, &status); /*exposure    */
      fits_write_2d_flt(outfptr, 0L, (long) imagebin, (long) imagebin,(long)imagebin,expo_image[ienergy-1],&status);
    }
    fits_movabs_hdu(outfptr,4*num_E_range+1+offset_first_output_image,&hdutype, &status); /*Exposure    */
    fits_write_2d_flt(outfptr, 0L, (long) imagebin, (long) imagebin,(long)imagebin,exp_image,&status);
    fits_close_file(outfptr, &status);
    
    /*free memory buffer*/
    for(i=0;i<num_E_range;i++){
      free(in_flx_image[i]);
      free(in_var_image[i]);
      free(in_expo_image[i]);
      free(flx_image[i]);
      free(sig_image[i]);
      free(var_image[i]);
      free(expo_image[i]);
    }
    free(exp_image);
    
// write hp
    for (i_hp16=0;i_hp16<nside2npix(16);i_hp16++) {
        if (nscw_hp16[0][i_hp16]==0) continue;
    
        char fn[255];
        char f_fn[255];
        char dir1[255];
        char dir4[255];
        int i_scw;
        int i_s;


        pix2ang_ring(16,i_hp16,&theta,&phi);
        ang2pix_ring(4,theta,phi,&i_hp4);
        ang2pix_ring(1,theta,phi,&i_hp1);

        sprintf(dir1,"%s/hp_%li",getenv("HP_DATA"),i_hp1);
        sprintf(dir4,"%s/hp_%li",dir1,i_hp4);
        sprintf(fn,"%s/hp_%li.txt",dir4,i_hp16);
        sprintf(f_fn,"!%s/hp_%li.fits",dir4,i_hp16);
        mkdir(dir1,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        mkdir(dir4,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        
        printf("writing %s of %i\n",f_fn,nscw_hp16[0][i_hp16]);

        writebintable(f_fn,nscw_hp16[0][i_hp16],ijd_hp16[0][i_hp16],inv_hp1024_in_hp16+i_hp16*4096,flux_hp16[0][i_hp16],var_hp16[0][i_hp16],expo_hp16[0][i_hp16],nsubpx_hp16[0][i_hp16]);

        if (0) {
            printf("writing %s of %i\n",fn,nscw_hp16[0][i_hp16]);

            FILE * f = fopen(fn,"w");
            for (i_scw=0;i_scw<nscw_hp16[0][i_hp16];i_scw++) {
                for(ienergy=0;ienergy<num_E_range;ienergy++){
                    for (i_s=0;i_s<4096;i_s++) {
                        //printf("%i %i %i\n",ienergy,i_scw,i);
                        i_hp1024=inv_hp1024_in_hp16[i_hp16*4096+i_s];
                        pix2ang_ring(1024,i_hp1024,&theta,&phi);
                        long int i_tmp;
                        ang2pix_ring(16,theta,phi,&i_tmp);
                        dec=theta/M_PI*180.-90.;
                        ra=phi/M_PI*180.;
                        fprintf(f,"%i %i %.15g %i %li %.5g %.5g %.5g %.5g %.5g %li\n",ienergy,i_scw,ijd_hp16[ienergy][i_hp16][i_scw],i,i_hp1024,ra,dec,flux_hp16[ienergy][i_hp16][i_scw*496+i_s],var_hp16[ienergy][i_hp16][4096*i_scw+i_s],i_tmp);
                    };
                };
            };
            
            fclose(f);
        };
    };
// done write hp
  }/*end of block ZZZ finish writing images*/
  
  /* Close the output file. */
  if(status == 0){
    sprintf(text,"varmosaic ends status = %d\n", status);
    HD_printf(text);
  }else{
    fprintf(stderr,"varmosaic ends with error.\n");
    sprintf(text,"status = %d\n", status);
    HD_printf(text);
  }
  return status;
}

int write_wcs(fitsfile *outfptr,double xrval,double yrval,double xrpix,double yrpix,double xinc,double yinc,double rot, char ctype1[], char ctype2[], int *status){
  fits_write_key_str(outfptr,"RADECSYS","FK5","Reference frame",status);
  fits_write_key_str(outfptr,"CTYPE1",ctype1,"Projection method", status);
  fits_write_key_str(outfptr,"CTYPE2",ctype2,"Projection method", status);
  fits_write_key_dbl(outfptr,"EQUINOX",2000.0E0,-15,"Coordinate system equinox",status);
  fits_write_key_dbl(outfptr,"CRVAL1",xrval,-15,"RA at reference point", status);
  fits_write_key_dbl(outfptr,"CRVAL2",yrval,-15,"DEC at reference point",status);  
  fits_write_key_dbl(outfptr,"CRPIX1",xrpix,-15,"X at reference point",status);  
  fits_write_key_dbl(outfptr,"CRPIX2",yrpix,-15,"Y at reference point", status);  
  fits_write_key_dbl(outfptr,"CDELT1",xinc, -15,"X increment (deg)",status);  
  fits_write_key_dbl(outfptr,"CDELT2",yinc, -15,"Y increment (deg)",status);  
  fits_write_key_dbl(outfptr,"CROTA2",rot, -15,"Rotation",status);  
}

/*--------------------------------------------------------------------------*/
void printerror( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/


    if (status)
    {
       fits_report_error(stderr, status); /* print error report */

       exit( status );    /* terminate the program, returning error status */
    }
    return;
}

void writebintable ( char fn[], int n_scw, double * ijds, long int * i_hp1024, double *flux, double *variance, double *exposure, int *nsubpx ) {
    fitsfile *fptr;       
    int status, hdutype;
    long firstrow, firstelem;

    int tfields   = 1; 

    char extname[] = "IJD";           /* extension name */
    /* define the name, datatype, and physical units for the 3 columns */
    char *ttype[] = { "IJD" };
    char *tform[] = { "1E"  };
    char *tunit[] = { "days" };

    /* define the name diameter, and density of each planet */

    status=0;

    printf("f 1\n");

    /* open the FITS file containing a primary array and an ASCII table */
    if ( fits_create_file(&fptr, fn, &status) ) 
         printerror( status );
    
    printf("f 1.\n");
    

    long *paxes[] = { 10 ,10  };
    if ( fits_create_img(fptr, USHORT_IMG, 2, paxes, &status) ) 
         printerror( status );

    //if ( fits_open_file(&fptr, fn, READWRITE, &status) ) 
      //   printerror( status );
    printf("f 2\n");

    //if ( fits_movabs_hdu(fptr, 1, &hdutype, &status) ) /* move to 2nd HDU */
    //     printerror( status );
    //printf("f 3\n");

    /* append a new empty binary table onto the FITS file */
    if ( fits_create_tbl( fptr, BINARY_TBL, n_scw, tfields, ttype, tform,
                tunit, extname, &status) )
         printerror( status );

    firstrow  = 1;  /* first row in table to write   */
    firstelem = 1;  /* first element in row  (ignored in ASCII tables) */

    /* write names to the first column (character strings) */
    /* write diameters to the second column (longs) */
    /* write density to the third column (floats) */

    fits_write_col(fptr, TDOUBLE, 1, firstrow, firstelem, n_scw, ijds,
                   &status);
    
    tfields   = 3; 
    
    printf("f 3\n");
    
    long *axis[] = { 4096 };
    if ( fits_create_img(fptr, ULONG_IMG, 1, axis, &status) ) 
         printerror( status );
    if ( fits_write_img(fptr, TULONG, 1, 4096, i_hp1024, &status) )
        printerror( status );

    long *axes[] = { 4096, n_scw };
    if ( fits_create_img(fptr, DOUBLE_IMG, 2, axes, &status) ) 
         printerror( status );
    if ( fits_write_img(fptr, TDOUBLE, 1, 4096*n_scw, flux, &status) )
        printerror( status );
    
    printf("f 4\n");

    if ( fits_create_img(fptr, DOUBLE_IMG, 2, axes, &status) ) 
         printerror( status );
    if ( fits_write_img(fptr, TDOUBLE, 1, 4096*n_scw, variance, &status) )
        printerror( status );
    


    if ( fits_create_img(fptr, DOUBLE_IMG, 2, axes, &status) ) 
         printerror( status );
    if ( fits_write_img(fptr, TDOUBLE, 1, 4096*n_scw, exposure, &status) )
        printerror( status );
    
    if ( fits_create_img(fptr, SHORT_IMG, 2, axes, &status) ) 
         printerror( status );
    if ( fits_write_img(fptr, TSHORT, 1, 4096*n_scw, nsubpx, &status) )
        printerror( status );
    
    if ( fits_close_file(fptr, &status) )       /* close the FITS file */
         printerror( status );
    

    return;
}

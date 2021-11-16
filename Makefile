HD_COMPONENT_NAME	= integral

HD_COMPONENT_VERS	= 

HD_CTASK		= varmosaic

HD_CTASK_SRC_c		= ${HD_CTASK}.c

HD_CFLAGS		= ${HD_STD_CFLAGS} -g -I $(HEALPIX)/include

HD_CLIBS		= ${HD_STD_CLIBS} -g -L $(HEALPIX)/lib -lchealpix

HD_INSTALL_TASKS	= ${HD_CTASK}

HD_INSTALL_PFILES	= ${HD_CTASK}.par

HD_INSTALL_HELP		= ${HD_CTASK}.txt


include ${HD_STD_MAKEFILE}


#ifndef HEALPIX
#.ABORT
#endif

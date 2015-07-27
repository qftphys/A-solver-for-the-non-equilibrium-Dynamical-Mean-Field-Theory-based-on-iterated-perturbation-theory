##$ COMPILER: suppprted compilers are ifort, gnu >v4.7
FC=gfortran

##$ PLATFORM: supported platform are intel, gnu
PLAT=gnu


##$ LOCATION OF THE scifor and dmft_tools DIRECTORIES
##$ is placed in different locations edit INCARGS here below
LIBDIR=/opt



##$ CHOOSE THE MODEL BY SELECTING THE PROGRAM DRIVER
#EXE=neqdmft_bethe_quench
EXE=neqdmft_bethe_dos_quench
#EXE=neqdmft_2dsquare_quench
#EXE=neqdmft_2dsquare_field
#EXE=neqdmft_attractive_bethe_dos_quench

##$ SET THE LOCATION OF YOU PROGRAM DRIVER (default is ./drivers)
DIR =drivers

##$ SET THE LOCATION WHERE TO PLACE THE EXECUTABLE (default is $HOME/.bin)
DIREXE=$(HOME)/.bin


ifeq ($(PLAT),intel)
FFLAG +=
endif

ifeq ($(PLAT),gnu)
INCARGS=-I$(LIBDIR)/scifor/gnu/include -I$(LIBDIR)/dmft_tools/gnu/include
FFLAG +=-ffree-line-length-none $(INCARGS)
endif




##$ CHOOSE LINKING OPTIONS:
##$ 
##$ If you intend to use mkl:
##$ 
#MKLARGS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
#ARGS=-ldmftt -lscifor $(MKLARGS) 
##$ 
##$ ELSE:
##$ 
ARGS=  -ldmftt -lscifor -lfftpack -llapack -lblas




##$ REVISION SOFTWARE VARIABLES
##$ 
REV=$(shell git rev-parse HEAD)
BRANCH=_$(shell git rev-parse --abbrev-ref HEAD)
REV=$(shell git rev-parse HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

ifeq ($(BRANCH),_master)
BRANCH=
endif


##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90

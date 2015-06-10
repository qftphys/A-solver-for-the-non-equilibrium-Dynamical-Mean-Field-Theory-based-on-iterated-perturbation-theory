FC=gfortran
PLAT=gnu
#PRECOMPILATION FLAG (leave blank for serial code)
FPP=




#set here the absolute bath to the DMFT_TOOLS and SCIFORTRAN libraries:
LIBROOT=/scratch/fperonaci/opt

.SUFFIXES: .f90


ifeq ($(PLAT),intel)
FFLAG +=-fpp -D_$(FPP_INEQ) -D_$(FPP_ED)
endif
ifeq ($(PLAT),gnu)
INCARGS=-I$(LIBROOT)/scifor/gnu/include -I$(LIBROOT)/dmft_tools/gnu/include
FFLAG +=-ffree-line-length-none -cpp -D_$(FPP_INEQ) -D_$(FPP_ED) $(INCARGS)
endif

#CHOOSE LINKING OPTIONS:
#if you intend to use mkl:
#MKLARGS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
#ARGS=-ldmftt -lscifor $(MKLARGS) 
#ELSE:
ARGS=  -ldmftt -lscifor -lfftpack -llapack -lblas

#REVISION SOFTWARE GIT:
REV=_$(shell git rev-parse HEAD)
BRANCH=_$(shell git rev-parse --abbrev-ref HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

ifeq ($(BRANCH),_master)
BRANCH=
endif

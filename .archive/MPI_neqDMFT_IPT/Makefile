SHELL = /bin/bash
FC    = mpif90
HERE  =`pwd`
#CODE TO COMPILE
EXE   = freeNEQ
DIREXE= $(HERE)/test

#########################################################################

#COMPILER FLAGS:
OPT = -fpp -D__mix -axT -ftz -ip -march=core2 -mtune=core2 -O3 -xT  -parallel -static #>& error_output
OPT = -ftz -O3 -static -assume nobuffered_io -ip -parallel -xsse4.1 -no-prec-div -no-prec-sqrt -openmp  -std03 >& error_output

#OPT=-fpp -D__mix -O3 -assume nobuffered_io -static-intel -cxxlib -static-libgcc -static
DEB =-fpp -D__mix -p  -check all -O0 -debug -g -fpe0 -traceback -static-intel #-cxxlib -static-libgcc

#########################################################################

#LIBRARIES: (here you  add all the libraries, not only yours)
X11     = -L$(XLIB)                     -lX11 -lXt -lXext -lXau -lXdmcp # -lxcb -lX11-xcb -lxcb-xlib 
X11     = -L$(XLIB) -L/usr/lib64 -lX11 -lXt -lXext -lxcb -lX11-xcb -lxcb-xlib -lXau -lXdmcp
ifeq ($(ARCH), 32)
MKL     = -L$(MKLDIR)/lib/$(ARCH)       -lmkl_lapack95 -lmkl_blas95 -lpthread -lmkl_ia32 -lmkl_solver -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5
else
MKL     = -L$(MKLDIR)/lib/$(ARCH)       -lmkl_lapack95 -lmkl_blas95 -lpthread -lmkl_em64t -lmkl_solver -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
endif
DSL     = -L${DISLIN}                   -ldislin
FFT     = -L$(HOME)/lib/FFTW/lib    	-lfftw 
DLPLOT  = -L$(HOME)/lib/DLPLOT/lib 	-ldlplot
GRIDS   = -L$(HOME)/lib/GRIDS/lib 	-lgrids
LATTICE = -L$(HOME)/lib/LATTICE/lib 	-llattice
TOOLS   = -L$(HOME)/lib/TOOLS/lib 	-ltools
CGFIT   = -L$(HOME)/lib/CGFIT/lib 	-lcgfit
BAR     = -L$(HOME)/lib/PROGRESS/lib 	-lbar
MEM     = -L$(HOME)/lib/MEMORY/lib 	-lmem
NR77    = -L$(HOME)/lib/NR/F77/lib 	-lnumrecf77
#NR90   = -L$(HOME)/lib/NR/F90/lib 	-lnumrecf90
ifdef $(CUDA)
GPU     = -L$(HOME)/lib/LIBGPU/lib 	-lgpu
CUDA    = -L/opt/cuda_2.3/lib64 -lcufft -lcuda -lcudart -lcublas
CUDA_S  = -L/opt/cuda_sdk_2.3/lib -lcutil
#CULA = -L$(HOME)/MYLIBRARY/LIB_CULA_PREMIUM_1.2/lib64/ -lcula -lcula_fortran
else	
GPU     = 
CUDA    = 
CUDA_S  = 
#CULA   = 
endif
#Gather them in a single variable
#NB: personal libraries making use of MKL/DISLIN/CUDA/etc. *must* be above(ie left) of their respective general libs   
LIBS = 	${FFT} \
	${TOOLS} \
	${MKL} \
	${DLPLOT} \
	${DSL} \
	${CGFIT} \
	${BAR} \
	${MEM} \
	${GRIDS} \
	${LATTICE} \
	${X11} \
	${GPU} \
	${CUDA_STATIC} \
	${CUDA} 

#########################################################################
#MODULES: here are the path to the .mod files, not only yours
MODS = 	-I$(HOME)/lib/FFTW/include		\
	-I$(HOME)/lib/MEMORY/include 		\
	-I$(HOME)/lib/PROGRESS/include 		\
	-I$(HOME)/lib/TOOLS/include 		\
	-I$(HOME)/lib/LATTICE/include 		\
	-I$(HOME)/lib/GRIDS/include		\
	-I$(HOME)/lib/mkl_include		\
	-I${DISLIN}/ifc				\
	-I$(HOME)/lib/DLPLOT/include		\
	-I$(HOME)/lib/LIBGPU/include		\
#       -I$(HOME)/lib/NR/F77/include 


#########################################################################
#COMPILATION:
OBJS  = VARS_GLOBAL.o EQKELDYSH.o FUNX_NEQ.o  KADANOFBAYM.o
OBJS_DEB = VARS_GLOBAL_DEB.o EQKELDYSH_DEB.o FUNX_NEQ_DEB.o  KADANOFBAYM_DEB.o

all: 	$(OBJS)
	@echo " ........... compile: optimized ........... "
	$(FC)  $(OPT) $(OBJS) $(EXE).f90 -o $(DIREXE)/$(EXE) $(MODS) $(LIBS)
	@echo " ...................... done .............................. "
VARS_GLOBAL.o: VARS_GLOBAL.f90
	$(FC) $(OPT) -c VARS_GLOBAL.f90 $(MODS)
EQKELDYSH.o: EQKELDYSH.f90 
	$(FC) $(OPT) -c EQKELDYSH.f90 $(MODS)
FUNX_NEQ.o: FUNX_NEQ.f90
	$(FC) $(OPT) -c FUNX_NEQ.f90 $(MODS)
KADANOFBAYM.o: KADANOFBAYM.f90 
	$(FC) $(OPT) -c KADANOFBAYM.f90 $(MODS)



debug: 	$(OBJS_DEB)
	@echo " ........... compile: debug   ........... "
	$(FC) $(DEB) $(OBJS) $(EXE).f90 -o $(DIREXE)/$(EXE)_deb $(MODS) $(LIBS)
	@echo " ...................... done .............................. "
VARS_GLOBAL_DEB.o: VARS_GLOBAL.f90
	$(FC) $(DEB) -c VARS_GLOBAL.f90 $(MODS)
EQKELDYSH_DEB.o: EQKELDYSH.f90 
	$(FC) $(DEB) -c EQKELDYSH.f90 $(MODS)
FUNX_NEQ_DEB.o: FUNX_NEQ.f90
	$(FC) $(DEB) -c FUNX_NEQ.f90 $(MODS)
KADANOFBAYM_DEB.o: KADANOFBAYM.f90 
	$(FC) $(DEB) -c KADANOFBAYM.f90 $(MODS)


clean: 
	@echo "Cleaning:"
	@rm -f *.mod
	@rm -f *.o
	@rm -f *~

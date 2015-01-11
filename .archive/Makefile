
#FILE PER COMPILARE IL CODICE tDMFT
################################################

FC=ifort  
HERE=`pwd`
#EXE=neqIPT
EXE=freeNEQ
DIREXE=$(HERE)/test
NUM=`ls -l|wc -l`

OPT= -axT -ftz -ip -march=core2 -mtune=core2 -O3 -xT -parallel 
#OPT=-O3
#OPT= -axT -ftz -ip  -march=core2 -mtune=core2 -O3 -xT -parallel -ipo${NUM}
#OPT=-O0 -check -debug full -g -fpe0 -traceback #-inline-debug-info
#OPT=-check -traceback -g -O0
#OPT= -O0

LIBS=-I${DISLIN}/ifc -I$(MKLDIR)/include -L${DISLIN} -L$(MKLDIR)/lib/$(ARCH) -lmkl_lapack95 -lmkl_blas95 -lmkl -lguide -ldislin

OBJS=FFTW.o VARS_GLOBAL.o FUNX_GLOBAL.o FUNX_NEQ.o EQHMSOLUTION.o  KADANOFBAYM.o

all: $(OBJS)
	$(FC) $(OPT) $(OBJS) $(EXE).f90 -o $(DIREXE)/$(EXE) $(LIBS)

VARS_GLOBAL.o: VARS_GLOBAL.f90
	$(FC) $(OPT) -c VARS_GLOBAL.f90

FUNX_GLOBAL.o: FUNX_GLOBAL.f90
	$(FC) $(OPT) -c FUNX_GLOBAL.f90

FUNX_NEQ.o: FUNX_NEQ.f90
	$(FC) $(OPT) -c FUNX_NEQ.f90

EQHMSOLUTION.o: EQHMSOLUTION.f90
	$(FC) $(OPT) -c EQHMSOLUTION.f90

FFTW.o: FFTW.f90
	$(FC) $(OPT) -c FFTW.f90

KADANOFBAYM.o: KADANOFBAYM.f90
	$(FC) $(OPT) -c KADANOFBAYM.f90

clean: 
	@echo 'removing *.mod *.o *~'
	@rm -f *.mod
	@rm -f *.o
	@rm -f *~

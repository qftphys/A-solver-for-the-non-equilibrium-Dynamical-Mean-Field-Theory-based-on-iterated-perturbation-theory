#FC=ifort  
FC=mpif77

HERE=`pwd`
MAIN_CODE=neq_keldysh.f90
DIREXE=$(HERE)/test
OUTEXE=neqDMFT

OPT= -axT -ftz -ip -ipo -march=core2 -mtune=core2 -O3 -xT 
DEBUG=-O0 -warn -check -debug full -g -fpe0 -traceback #-inline-debug-info

LIBS=-I/opt/Intel/mkl/10.1.3.027/include -L/opt/Intel/mkl/10.1.3.027/lib/em64t -lmkl_lapack95 -lmkl_blas95 -lmkl -lguide

OBJS=VARIABILI.o FUNZIONI.o GLOC.o

all: $(OBJS)
	$(FC) $(OPT) $(OBJS) $(MAIN_CODE) -o $(DIREXE)/$(OUTEXE) $(LIBS)

VARIABILI.o: VARIABILI.f90
	$(FC) $(OPT) -c VARIABILI.f90

FUNZIONI.o: FUNZIONI.f90
	$(FC) $(OPT) -c FUNZIONI.f90

GLOC.o: GLOC.f90
	$(FC) $(OPT) -c GLOC.f90


clean: 
	@echo 'removing *.mod *.o *~'
	@rm -f *.mod
	@rm -f *.o
	@rm -f *~

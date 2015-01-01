FC=gfortran
#PRECOMPILATION FLAG (leave blank for serial code)
FPP=

#EXE=neqdmft_bethe
#EXE=neqdmft_bethe_2bands
#EXE=neqdmft_lattice_2bands
#EXE=neqdmft_bethe_dos
#EXE=neqdmft_hypercubic
EXE=neqdmft_2dsquare_quench
#EXE=neqdmft_2dsquare_field



DIR=./drivers
DIREXE= $(HOME)/.bin
.SUFFIXES: .f90

#REVISION SOFTWARE GIT:
BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc


OBJS =  CONTOUR_GF.o NEQ_VARS_GLOBAL.o ELECTRIC_FIELD.o NEQ_THERMOSTAT.o NEQ_IPT.o

#MKLARGS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm

INCARGS=-I/opt/scifor/gnu/include -I/opt/dmft_tools/gnu/include
FFLAG +=-ffree-line-length-none $(INCARGS)

#ARGS=-ldmftt -lscifor $(MKLARGS) -lminpack -larpack -lparpack 
ARGS= -ldmftt -lscifor -lfftpack -llapack -lblas -lminpack -larpack -lparpack

all:compile

compile: version $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC) $(FFLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)_$(BRANCH) $(ARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_$(BRANCH)


.f90.o:	
	$(FC) $(FFLAG) -c $<


completion:
	src_completion.sh $(DIR)/$(EXE).f90
	@echo "run: . .bash_completion.d/$(EXE) to add completion for $(EXE) in this shell"

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc

version:
	@echo $(VER)

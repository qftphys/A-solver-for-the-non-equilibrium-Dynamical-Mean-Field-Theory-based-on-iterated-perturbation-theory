#=========================================================================
include sfmake.inc
#=========================================================================
#EXE=neqdmft_bethe
#EXE=neqdmft_bethe_dos
#EXE=neqdmft_hypercubic
#EXE=neqdmft_2dsquare_quench
EXE=neqdmft_2dsquare_field
DIR=./drivers
DIREXE= $(HOME)/.bin
FC=ifort
BRANCH= $(shell git rev-parse --abbrev-ref HEAD)

OBJS =  CONTOUR_GF.o NEQ_VARS_GLOBAL.o ELECTRIC_FIELD.o NEQ_THERMOSTAT.o NEQ_IPT.o


#=================STANDARD COMPILATION====================================
all: FLAG=$(STD)
all: ARGS=$(SFLIBS)
all:compile

#================OPTIMIZED COMPILATION====================================
opt: FLAG=$(OPT)
opt: ARGS=$(SFLIBS)
opt:compile

#================DEBUGGIN COMPILATION=====================================
debug: FLAG=$(DEB)
debug: ARGS=$(SFLIBS_DEB)
debug:compile

compile: version $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC) $(FLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)_$(BRANCH) $(ARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_$(BRANCH)


.f90.o:	
	$(FC) $(FLAG) -c $< $(SFINCLUDE) 

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc

version:
	@echo $(VER)

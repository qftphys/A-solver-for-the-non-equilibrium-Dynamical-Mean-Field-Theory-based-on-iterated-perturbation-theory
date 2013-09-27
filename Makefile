#=========================================================================
include sfmake.inc
#=========================================================================
FC=ifort
EXE   = neqDMFT
DIR=./drivers
DIREXE= $(HOME)/.bin
BRANCH= $(shell git rev-parse --abbrev-ref HEAD)

OBJS =  VIDE_RK24.o CONTOUR_GF.o NEQ_VARS_GLOBAL.o ELECTRIC_FIELD.o NEQ_BATH.o NEQ_IPT.o NEQ_KADANOFF_BAYM.o NEQ_RESULTS.o


FLAG=$(STD)
#FLAG=$(OPT)
ARGS= $(SFLIBS)

#FLAG=$(DEB)
#ARGS= $(SFLIBS_DEB)


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

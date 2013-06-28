#=========================================================================
include sfmake.inc
#=========================================================================
EXE= neqDMFT
#EXE=get_data_neqDMFT
DIR=./drivers
DIREXE= $(HOME)/.bin

FLAG=$(STD)
#FLAG=${DEB}
#FLAG=$(OPT)
ARGS=$(SFLIBS)


BRANCH=  $(shell git rev-parse --abbrev-ref HEAD)

OBJS =  CONTOUR_GF.o NEQ_VARS_GLOBAL.o ELECTRIC_FIELD.o BATH.o EQ_IPT.o NEQ_IPT.o NEQ_UPDATE_WF.o KADANOFBAYM.o NEQ_KADANOFF_BAYM.o RESULTS.o

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

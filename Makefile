#=========================================================================
include sfmake.inc
#=========================================================================
FC=$(SFMPI)/mpif90
EXE   = neqDMFT
#EXE=get_data_neqDMFT
DIREXE= $(HOME)/.bin
BRANCH=  $(shell git rev-parse --abbrev-ref HEAD)
.SUFFIXES: .f90 
OBJS =  CONTOUR_GF.o VARS_GLOBAL.o ELECTRIC_FIELD.o BATH.o EQUILIBRIUM.o IPT_NEQ.o UPDATE_WF.o KADANOFBAYM.o RESULTS.o


FLAG=$(STD)
#FLAG=$(DEB)
#FLAG=$(OPT)
ARGS= $(SFLIBS)
#ARGS= $(SFLIBS_DEB)


compile: version $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC) $(FLAG) $(OBJS) $(EXE).f90 -o $(DIREXE)/$(EXE)_$(BRANCH) $(ARGS)
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

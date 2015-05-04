include make.inc

#EXE=neqdmft_bethe_quench
EXE=neqdmft_bethe_dos_quench
#EXE=neqdmft_2dsquare_quench
#EXE=neqdmft_2dsquare_field

DIR=./drivers
DIREXE=.



OBJS =  NEQ_CONTOUR.o NEQ_CONTOUR_GF.o NEQ_INPUT_VARS.o ELECTRIC_FIELD.o NEQ_THERMOSTAT.o NEQ_EQUILIBRIUM.o  NEQ_MEASURE.o NEQ_IPT.o NEQ_DMFT_IPT.o


all:compile

compile: version $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC) $(FFLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)$(BRANCH) $(ARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)$(BRANCH)


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

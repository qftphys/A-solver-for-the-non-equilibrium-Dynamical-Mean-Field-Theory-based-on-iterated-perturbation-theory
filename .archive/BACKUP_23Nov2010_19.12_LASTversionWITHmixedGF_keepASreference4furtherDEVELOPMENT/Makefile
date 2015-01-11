HERE  =`pwd`
EXE   = neqDMFT
DIREXE= $(HOME)/.bin
FC  = ifort

#=========================================================================
STD = -fpp -D_noD
OPT = -fpp -D_noD -O3 -ftz -ip -assume nobuffered_io -openmp -parallel -static
DEB = -fpp -D_noD -p  -check all -O0 -debug -g -fpe0 -traceback -static-intel #-warn unused
SHKD = -fpp -D_mix
OHKD = -fpp -D_mix -O3 -ftz -ip -assume nobuffered_io -openmp -parallel -static 
DHKD = -fpp -D_mix -p  -check all -O0 -debug -g -fpe0 -traceback -static-intel -warn unused

#=========================================================================
include $(HOME)/lib/lib.mk


#=========================================================================
OBJS     = VARS_GLOBAL.o EQKELDYSH.o FUNX_NEQ.o KADANOFBAYM.o OBSERVABLES.o
OBJS_OPT = VARS_GLOBAL_OPT.o EQKELDYSH_OPT.o FUNX_NEQ_OPT.o KADANOFBAYM_OPT.o OBSERVABLES_OPT.o
OBJS_DEB = VARS_GLOBAL_DEB.o EQKELDYSH_DEB.o FUNX_NEQ_DEB.o KADANOFBAYM_DEB.o OBSERVABLES_DEB.o
OBJS_SHKD = VARS_GLOBAL_SHKD.o EQKELDYSH_SHKD.o FUNX_NEQ_SHKD.o KADANOFBAYM_SHKD.o OBSERVABLES_SHKD.o
OBJS_OHKD = VARS_GLOBAL_OHKD.o EQKELDYSH_OHKD.o FUNX_NEQ_OHKD.o KADANOFBAYM_OHKD.o OBSERVABLES_OHKD.o
OBJS_DHKD = VARS_GLOBAL_DHKD.o EQKELDYSH_DHKD.o FUNX_NEQ_DHKD.o KADANOFBAYM_DHKD.o OBSERVABLES_DHKD.o


#=================STANDARD COMPILATION====================================
all: 	$(OBJS)
	@echo " ........... compile: normal ........... "
	$(FC)  $(STD) $(OBJS) $(EXE).f90 -o $(DIREXE)/$(EXE) $(MODS) $(LIBS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)

VARS_GLOBAL.o: VARS_GLOBAL.f90
	$(FC) $(STD) -c VARS_GLOBAL.f90 $(MODS)
EQKELDYSH.o: EQKELDYSH.f90 
	$(FC) $(STD) -c EQKELDYSH.f90 $(MODS)
OBSERVABLES.o: OBSERVABLES.f90 
	$(FC) $(STD) -c OBSERVABLES.f90 $(MODS)
FUNX_NEQ.o: FUNX_NEQ.f90
	$(FC) $(STD) -c FUNX_NEQ.f90 $(MODS)
KADANOFBAYM.o: KADANOFBAYM.f90 
	$(FC) $(STD) -c KADANOFBAYM.f90 $(MODS)


hookd: 	$(OBJS_SHKD)
	@echo " ........... compile: hooked   ........... "
	$(FC) $(SHKD) $(OBJS_SHKD) $(EXE).f90 -o $(DIREXE)/$(EXE)_hkd $(MODS) $(LIBS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_hkd
VARS_GLOBAL_SHKD.o: VARS_GLOBAL.f90
	$(FC) $(SHKD) -c VARS_GLOBAL.f90 $(MODS) -o VARS_GLOBAL_SHKD.o
EQKELDYSH_SHKD.o: EQKELDYSH.f90 
	$(FC) $(SHKD) -c EQKELDYSH.f90 $(MODS) -o EQKELDYSH_SHKD.o
OBSERVABLES_SHKD.o: OBSERVABLES.f90 
	$(FC) $(SHKD) -c OBSERVABLES.f90 $(MODS) -o OBSERVABLES_SHKD.o
FUNX_NEQ_SHKD.o: FUNX_NEQ.f90 
	$(FC) $(SHKD) -c FUNX_NEQ.f90 $(MODS) -o FUNX_NEQ_SHKD.o
KADANOFBAYM_SHKD.o: KADANOFBAYM.f90 
	$(FC) $(SHKD) -c KADANOFBAYM.f90 $(MODS) -o KADANOFBAYM_SHKD.o




#================OPTIMIZED COMPILATION====================================
opt: 	$(OBJS_OPT)
	@echo " ........... compile: optimized   ........... "
	$(FC) $(OPT) $(OBJS_OPT) $(EXE).f90 -o $(DIREXE)/$(EXE) $(MODS) $(LIBS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)

VARS_GLOBAL_OPT.o: VARS_GLOBAL.f90
	$(FC) $(OPT) -c VARS_GLOBAL.f90 $(MODS) -o VARS_GLOBAL_OPT.o
EQKELDYSH_OPT.o: EQKELDYSH.f90 
	$(FC) $(OPT) -c EQKELDYSH.f90 $(MODS) -o EQKELDYSH_OPT.o
OBSERVABLES_OPT.o: OBSERVABLES.f90 
	$(FC) $(OPT) -c OBSERVABLES.f90 $(MODS) -o OBSERVABLES_OPT.o
FUNX_NEQ_OPT.o: FUNX_NEQ.f90 
	$(FC) $(OPT) -c FUNX_NEQ.f90 $(MODS) -o FUNX_NEQ_OPT.o
KADANOFBAYM_OPT.o: KADANOFBAYM.f90 
	$(FC) $(OPT) -c KADANOFBAYM.f90 $(MODS) -o KADANOFBAYM_OPT.o

hookd_opt: 	$(OBJS_OHKD)
	@echo " ........... compile: optimized hooked   ........... "
	$(FC) $(OHKD) $(OBJS_OHKD) $(EXE).f90 -o $(DIREXE)/$(EXE)_hkd $(MODS) $(LIBS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_hkd
VARS_GLOBAL_OHKD.o: VARS_GLOBAL.f90
	$(FC) $(OHKD) -c VARS_GLOBAL.f90 $(MODS) -o VARS_GLOBAL_OHKD.o
EQKELDYSH_OHKD.o: EQKELDYSH.f90 
	$(FC) $(OHKD) -c EQKELDYSH.f90 $(MODS) -o EQKELDYSH_OHKD.o
OBSERVABLES_OHKD.o: OBSERVABLES.f90 
	$(FC) $(OHKD) -c OBSERVABLES.f90 $(MODS) -o OBSERVABLES_OHKD.o
FUNX_NEQ_OHKD.o: FUNX_NEQ.f90 
	$(FC) $(OHKD) -c FUNX_NEQ.f90 $(MODS) -o FUNX_NEQ_OHKD.o
KADANOFBAYM_OHKD.o: KADANOFBAYM.f90 
	$(FC) $(OHKD) -c KADANOFBAYM.f90 $(MODS) -o KADANOFBAYM_OHKD.o


#================DEBUGGIN COMPILATION=====================================
debug: 	$(OBJS_DEB)
	@echo " ........... compile: debug   ........... "
	$(FC) $(DEB) $(OBJS_DEB) $(EXE).f90 -o $(DIREXE)/$(EXE)_deb $(MODS) $(LIBS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_deb
VARS_GLOBAL_DEB.o: VARS_GLOBAL.f90
	$(FC) $(DEB) -c VARS_GLOBAL.f90 $(MODS) -o VARS_GLOBAL_DEB.o
EQKELDYSH_DEB.o: EQKELDYSH.f90 
	$(FC) $(DEB) -c EQKELDYSH.f90 $(MODS) -o EQKELDYSH_DEB.o
OBSERVABLES_DEB.o: OBSERVABLES.f90 
	$(FC) $(DEB) -c OBSERVABLES.f90 $(MODS) -o OBSERVABLES_DEB.o
FUNX_NEQ_DEB.o: FUNX_NEQ.f90 
	$(FC) $(DEB) -c FUNX_NEQ.f90 $(MODS) -o FUNX_NEQ_DEB.o
KADANOFBAYM_DEB.o: KADANOFBAYM.f90 
	$(FC) $(DEB) -c KADANOFBAYM.f90 $(MODS) -o KADANOFBAYM_DEB.o


hookd_deb: 	$(OBJS_DHKD)
	@echo " ........... compile: debug hooked   ........... "
	$(FC) $(DHKD) $(OBJS_DHKD) $(EXE).f90 -o $(DIREXE)/$(EXE)_deb $(MODS) $(LIBS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_deb
VARS_GLOBAL_DHKD.o: VARS_GLOBAL.f90
	$(FC) $(DHKD) -c VARS_GLOBAL.f90 $(MODS) -o VARS_GLOBAL_DHKD.o
EQKELDYSH_DHKD.o: EQKELDYSH.f90 
	$(FC) $(DHKD) -c EQKELDYSH.f90 $(MODS) -o EQKELDYSH_DHKD.o
OBSERVABLES_DHKD.o: OBSERVABLES.f90 
	$(FC) $(DHKD) -c OBSERVABLES.f90 $(MODS) -o OBSERVABLES_DHKD.o
FUNX_NEQ_DHKD.o: FUNX_NEQ.f90 
	$(FC) $(DHKD) -c FUNX_NEQ.f90 $(MODS) -o FUNX_NEQ_DHKD.o
KADANOFBAYM_DHKD.o: KADANOFBAYM.f90 
	$(FC) $(DHKD) -c KADANOFBAYM.f90 $(MODS) -o KADANOFBAYM_DHKD.o


#==============WIGNER TRANSFORMATION======================================
wigner: VARS_GLOBAL_W.o
	@echo " ........... compile: wigner ........... "
	${FC} VARS_GLOBAL_W.o ${OPT} WIGNER_postprocessing.f90 -o ${DIREXE}/wigner_postpro ${MODS} ${LIBS}
	@echo " ...................... done .............................. "

VARS_GLOBAL_W.o: VARS_GLOBAL.f90
	$(FC) $(OPT) -c VARS_GLOBAL.f90 $(MODS) -o VARS_GLOBAL_W.o

#=============CLEAN ALL===================================================
clean: 
	@echo "Cleaning:"
	@rm -f *.mod
	@rm -f *.o
	@rm -f *~
	@rm -vf	$(DIREXE)/$(EXE)*

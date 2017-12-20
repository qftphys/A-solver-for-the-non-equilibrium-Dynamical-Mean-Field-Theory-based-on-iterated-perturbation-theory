##$ COMPILER: suppprted compilers are ifort, gnu >v4.7 or use mpif90
FC=mpif90

##$ PLATFORM: supported platform are intel, gnu 
PLAT=gnu

##$ SET THE LOCATION OF YOU PROGRAM DRIVER (default is ./drivers)
DIR =./drivers

##$ SET THE TARGET DIRECTORY WHERE TO PUT THE EXECUTABLE (default if $HOME/.bin in the PATH)
DIREXE=$(HOME)/.bin

##$ CHOOSE THE DRIVER CODE:
EXE=neqdmft_2dsquare_field

##$ SET INCLUDE AND LINK OPTIONS USING pkg-config
INCARGS=$(shell pkg-config --cflags dmft_tools scifor)
LIBARGS=$(shell pkg-config --libs   dmft_tools scifor)


ifeq ($(PLAT),intel)
FFLAG=-O2 -ftz
OFLAG=-O3 -ftz
DFLAG=-p -O0 -g -fpe0 -warn -warn errors -debuEg extended -traceback -check all,noarg_temp_created
FPPFLAG =-fpp
endif

ifeq ($(PLAT),gnu)
FFLAG = -O2 -ffree-line-length-none
DFLAG = -O0 -p -g -fimplicit-none -Wsurprising  -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none
OFLAG = -O3 -ffast-math -march=native -funroll-all-loops -fno-protect-parens -flto -ffree-line-length-none
FPPFLAG =-cpp
endif



##$ REVISION SOFTWARE VARIABLES
REV=$(shell git rev-parse HEAD)
BRANCH=_$(shell git rev-parse --abbrev-ref HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

ifeq ($(BRANCH),_master)
BRANCH=
endif


##$ Extends the implicit support of the Makefile to .f90 files

.SUFFIXES: .f90

OBJS =  NEQ_INPUT_VARS.o NEQ_CONTOUR.o NEQ_CONTOUR_GF.o  NEQ_FFT_GF.o ELECTRIC_FIELD.o NEQ_THERMOSTAT.o NEQ_SETUP.o  NEQ_MEASURE.o NEQ_IPT_EQUILIBRIUM.o NEQ_IPT.o NEQ_DMFT_IPT.o



all: all version compile
debug: debug version compile

all: FPPFLAG+=-D_

debug: FFLAG=$(DFLAG)
debug: FPPFLAG+=-D_


compile: version $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC) $(FFLAG) $(INCARGS) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)$(BRANCH) $(LIBARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)$(BRANCH)


.f90.o:	
	$(FC) $(FFLAG) $(INCARGS) -c $<


completion:
	src_completion.sh $(DIR)/$(EXE).f90
	@echo "run: . .bash_completion.d/$(EXE) to add completion for $(EXE) in this shell"

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc

version:
	@echo $(VER)

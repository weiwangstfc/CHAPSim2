#################################################################################################
# Makefile for CHAPSim2, by Wei Wang, July 2021                                                 #
# Usage:                                                                                        #
#       make all        to make all files with -O2                                              #
#       make cfg=gnu    to debug for gfortran compiler                                          #
#       make cfg=intel  to debug for intel compiler                                             #
#       make cfg=cray   to make for cray                                                        #
# For debugging run:                                                                            #
# mpirun -np 4 valgrind --leak-check=full --track-origins=yes \                                 #
#                       --log-file=valgrind_output.txt ./CHAPSIM*                               #
#                          < solver_input > solver_output                                       #
# add flags for debug:
#    -DDEBUG_STEPS  : serial debug for each step                                                #
#    -DDEBUG_FFT    : fft debug                                                                 #
#################################################################################################
.SUFFIXES:
.PHONY: debug default clean

PROGRAM= CHAPSim

ifeq ($(cfg), gnu)
  FOPTS= -O -g -Wall -fbacktrace -fbounds-check -fcheck=all -ffpe-trap=invalid,zero,overflow \
   -finit-real=snan -ftrapv -ffree-line-length-512 -Wuninitialized -Wmaybe-uninitialized\
   -Wno-unused -fallow-argument-mismatch		   
  FFLGS= -DDOUBLE_PREC -fdefault-real-8 -fdefault-double-8
  FDEBG= -DDEBUG_STEPS -DDEBUG_ALGO# -DDEBUG_FFT # -DDEBUG_TEST -DDEBUG_FFT
else ifeq ($(cfg), intel)
  FOPTS= -g -assume ieee_fpe_flags -check all -check bounds -check uninit -debug all \
	-fp-stack-check fpe0 -fpe3 -fpe-all=3 -ftrapuv -ftz -warn all, nounused
  FFLGS= -DDOUBLE_PREC -DDEBUG
  FDEBG= -DDEBUG_STEPS -DDEBUG_FFT
else ifeq ($(cfg), cray)
  FOPTS= # -m 3
  FFLGS= # -s default64
  FDEBG= -DDEBUG_STEPS -DDEBUG_FFT
else ifeq ($(cfg), pg)
  FOPTS= -O3 -pg -march=native  -Wall -fimplicit-none  -ffree-line-length-512  -fwhole-file  -std=gnu \
	-ffpe-trap=invalid,zero,overflow -fall-intrinsics
  FFLGS= -DDOUBLE_PREC 
  FDEBG= -DDEBUG_STEPS -DDEBUG_FFT
else
  FOPTS= -O3  -march=native  -Wall -fimplicit-none  -ffree-line-length-512  -fwhole-file  -std=gnu \
	-ffpe-trap=invalid,zero,overflow -fall-intrinsics -fallow-argument-mismatch
  FFLGS= -DDOUBLE_PREC
  FDEBG= # -DDEBUG_STEPS # -DDEBUG_FFT -DDEBUG_VISU
endif

# this is based on the latest 2decomp&fft lib by UoE&ICL, 2022
#include ./lib/2decomp-fft-main/Makefile.settings
#INCLUDE= -I./lib/2decomp-fft-main/mod
#LIBS= -L./lib/2decomp-fft-main -ldecomp2d

# this is based on the original 2decomp&fft lib by NAG&ICL, 2012
#include ./lib/2decomp_fft/src/Makefile.inc
#INCLUDE= -I./lib/2decomp_fft/include
#LIBS= -L./lib/2decomp_fft/lib -l2decomp_fft

# this is based on the updated 2decomp&fft lib by NAG&ICL, 2013
include ./lib/2decomp_fft_updated/src/Makefile.inc
INCLUDE= -I./lib/2decomp_fft_updated/include
LIBS= -L./lib/2decomp_fft_updated/lib -l2decomp_fft

DIR_SRC= ./src
DIR_BIN= ./bin
DIR_OBJ= ./obj

OBJS1= mpi_mod.o\
      modules.o\
      tools_general.o\
      input_thermo.o\
      boundary_conditions.o\
      input_general.o\
      algorithms.o\
      operations2.o\
      tools_solver.o\
      geometry.o\
      io_tools.o\
      io_monitor.o\
      io_visulisation.o\
      statistics.o\
      domain_decomposition.o\
      poisson_interface.o\
      poisson_1stderivcomp.o\
      eq_continuity.o\
      eq_energy.o\
      eq_momentum.o\
      test_algrithms.o\
      io_restart.o\
      flow_initialization.o\
      chapsim.o
OBJS = $(OBJS1:%=$(DIR_OBJ)/%)

default :
	@cd $(DIR_BIN)
	make $(PROGRAM) -f Makefile
	@mv *.mod $(DIR_OBJ)
	@mv $(PROGRAM) $(DIR_BIN)
	@mv *.o $(DIR_OBJ)
	@echo "\n ======Successfully compiled====== \n"

$(PROGRAM): $(OBJS)
	@echo -n "\n ======Linking  ...====== \n"
	$(FC) $(FOPTS) $(FFLGS) $(FDEBG) -o $@ $(OBJS) $(LIBS)

$(DIR_OBJ)/%.o : $(DIR_SRC)/%.f90
	$(FC) $(INCLUDE) $(FOPTS) $(FFLGS) $(FDEBG) $(FCFLAGS)  $(FFLAGS) -c  -o $@ $<

all:
	@make clean
	@cd $(DIR_BIN)
	make $(PROGRAM) -f Makefile
	@mv *.mod $(DIR_OBJ)
	@mv $(PROGRAM) $(DIR_BIN)
	@echo "\n ======Successfully compiled====== \n"

clean:
	@rm -f $(DIR_OBJ)/*.o $(DIR_OBJ)/*.mod $(DIR_BIN)/$(PROGRAM)
	@rm -f *.mod *.o $(DIR_SRC)/*.mod $(DIR_SRC)/*.o

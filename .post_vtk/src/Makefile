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
#                                                                                               #
#################################################################################################

.PHONY: debug default clean
.SUFFIXES:

PROGRAM= post_vtk

ifeq ($(cfg), gnu)
	FOPTS= -g \
		   -fbacktrace \
		   -fbounds-check \
		   -fcheck=all \
	           -ffpe-trap=invalid,zero,overflow \
		   -finit-real=snan -ftrapv \
	           -Wall
	FFLGS= -DOUBLE_PREC -DDEBUG
else ifeq ($(cfg), intel)
	FOPTS= -g -assume ieee_fpe_flags -check all -check bounds -check uninit -debug all \
	-fp-stack-check fpe0 -fpe3 -fpe-all=3 -ftrapuv -ftz -warn all, nounused
	FFLGS= -DOUBLE_PREC -DDEBUG
	FOPTS= -O3  -march=native  -fimplicit-none  -Wall  -Wline-truncation  -fwhole-file  -std=gnu
else ifeq ($(cfg), cray)
	FOPTS= # -m 3
	FFLGS= # -s default64
else
	FOPTS= -O3  -march=native  -Wall -fimplicit-none  -Wline-truncation  -fwhole-file  -std=gnu \
	-ffpe-trap=invalid,zero,overflow -fall-intrinsics
	FFLGS= -DOUBLE_PREC
endif


include ./../../lib/2decomp_fft/src/Makefile.inc
INCLUDE = -I ./../../lib/2decomp_fft/include
LIBS = -L ./../../lib/2decomp_fft/lib -l2decomp_fft


OBJS= convert_decomp_data_to_vtk.o


default :
	make $(PROGRAM) -f Makefile

$(PROGRAM): $(OBJS)
	$(F90) -o $@ $(OBJS) $(FOPTS) $(FFLGS) $(LIBS)

%.o : %.f90
	$(F90) $(INCLUDE) $(FOPTS) $(FFLGS) $(F90FLAGS) -c $<

all:
	@make clean
	make $(PROGRAM) -f Makefile

clean:
	@rm -f *.o $(PROGRAM)
	@rm -f *.mod *.o



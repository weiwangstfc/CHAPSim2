.PHONY: debug default clean
.SUFFIXES:

PROGRAM= CHAPSim

include ./lib/2decomp_fft/src/Makefile.inc
INCLUDE = -I./lib/2decomp_fft/include
LIBS = -L./lib/2decomp_fft/lib -l2decomp_fft

DIR_SRC= ./src
DIR_BIN= ./bin
DIR_OBJ= ./obj
DIR_MOD= ./mod

OBJS= modules.o\
      mpi_mod.o\
      tools.o\
      input_general.o\
      input_thermo.o\
      geometry.o\
			flow_initialization.o\
			algorithms.o\
      chapsim.o

OPTIONS= -O -g -fbacktrace -fbounds-check -fcheck=all -Wall 

default :
	@cd $(DIR_BIN)
	make $(PROGRAM) -f Makefile
	@mv *.mod $(DIR_MOD)
	@mv *.o $(DIR_OBJ)
	@mv $(PROGRAM) $(DIR_BIN)

$(PROGRAM): $(OBJS)
	$(F90) -o $@ $(OBJS) $(LIBS)

%.o : $(DIR_SRC)/%.f90
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) -c $<

all:
	@make clean
	@cd $(DIR_BIN)
	make $(PROGRAM) -f Makefile
	@mv *.mod $(DIR_MOD)
	@mv *.o $(DIR_OBJ)
	@mv $(PROGRAM) $(DIR_BIN)

clean:
	@rm -f $(DIR_OBJ)/*.o $(DIR_BIN)/$(PROGRAM)
	@rm -f *.mod *.o $(DIR_SRC)/*.mod $(DIR_SRC)/*.o



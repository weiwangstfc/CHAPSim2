#include ../../src/Makefile.inc

#INCLUDE = -I../../include
FFLAGS := $(subst $(MODFLAG),$(MODFLAG)../../,$(FFLAGS))
FFLAGS := $(patsubst -I%,-I../../%,$(FFLAGS))
LIBS = -L../../ -l$(LIBDECOMP) $(LFLAGS)

NP ?= 1
MPIRUN ?= mpirun
NROW ?= 0
NCOL ?= 0

all: $(DECOMPINC) io_test io_read io_var_test io_plane_test io_bench

$(DECOMPINC):
	mkdir $@

io_test: io_test.o
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $< $(OBJ) $(LIBS)

io_read: io_read.o
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $< $(OBJ) $(LIBS)

io_var_test: io_var_test.o
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $< $(OBJ) $(LIBS)

io_plane_test: io_plane_test.o
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $< $(OBJ) $(LIBS)

io_bench: io_bench.o
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $< $(OBJ) $(LIBS)

ifeq ($(PARAMOD),gpu)
check:
	$(MPIRUN) -n $(NP) ./bind.sh ./io_test
	$(MPIRUN) -n $(NP) ./bind.sh ./io_read
	$(MPIRUN) -n $(NP) ./bind.sh ./io_var_test $(NROW) $(NCOL)
	$(MPIRUN) -n $(NP) ./bind.sh ./io_plane_test
	$(MPIRUN) -n $(NP) ./bind.sh ./io_bench
else
check:
	$(MPIRUN) -n $(NP) ./io_test
	$(MPIRUN) -n $(NP) ./io_read
	$(MPIRUN) -n $(NP) ./io_var_test $(NROW) $(NCOL)
	$(MPIRUN) -n $(NP) ./io_plane_test
	$(MPIRUN) -n $(NP) ./io_bench
endif

clean:
	rm -f *.o io_test io_read io_var_test io_plane_test io_bench

realclean: clean
	rm -f *.dat io_var_data.* *.log

%.o : %.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $< -o $@

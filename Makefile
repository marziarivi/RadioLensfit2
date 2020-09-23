###############################################################
#                                                             #
###############################################################

#-------------------------------------------------------------------------  
# OpenMP compiler switch
OMP = -fopenmp

OPT += -DUSE_MPI
OPT += -DFACET

#------------------------------------------------------------------------- 

ifeq (FACET,$(findstring FACET,$(OPT)))
 OPT += -DGRID
endif

# MERLINO
#CASACORE_INC =  -I/iranet/arcsoft/casacore/casacore-2.4.1/include -I/iranet/arcsoft/casacore/casacore-2.4.1/include/casacore
#CASACORE_LIB =  -L/iranet/arcsoft/casacore/casacore-2.4.1/lib -lcasa_casa -lcasa_measures -lcasa_tables -lcasa_ms -std=c++11

#GSL_INC = -I/iranet/arcsoft/gsl/gsl-2.5/include
#GSL_LIB = -L/iranet/arcsoft/gsl/gsl-2.5/lib


# SPLINTER
CASACORE_INC =  -I/share/splinter/cosmos/modules/sep_2019/install_dir/casacore-5.4.0/include -I/share/splinter/cosmos/modules/sep_2019/install_dir/casacore-5.4.0/include/casacore
CASACORE_LIB =  -L/share/splinter/cosmos/modules/sep_2019/install_dir/casacore-5.4.0/lib -lcasa_casa -lcasa_measures -lcasa_tables -lcasa_ms -std=c++11

GSL_INC = -I/share/splinter/cosmos/modules/sep_2019/install_dir/gsl-2.6/include
GSL_LIB = -L/share/splinter/cosmos/modules/sep_2019/install_dir/gsl-2.6/lib

SUP_INCL = -I. $(CASACORE_INC) $(GSL_INC)  


# INTEL compiler
#-------------------
OPTIMIZE = -O3 #bad results! -ip -ipo #-xAVX 
#OPTIMIZE = -O3 -ip -ipo  # using compilers/intel/15.1/133

ifeq (MPI,$(findstring MPI,$(OPT)))
#  CC  =  mpiicpc
else
#  CC  = icpc 
endif

# GNU compiler
# #-----------------
ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
  CC  =  mpiCC -g 
else
  CC  = g++
endif

OPTIONS = $(OPTIMIZE) $(OPT)

ifeq (MPI,$(findstring MPI,$(OPT)))
EXEC1 = RadioLensfit2-mpi
OBJS1  = RadioLensfit2-mpi.o galaxy_fitting.o utils.o read_catalog.o data_simulation.o distributions.o galaxy_visibilities.o  evaluate_uv_grid.o likelihood.o marginalise_r.o measurement_set.o 
else
EXEC1 = RadioLensfit2
OBJS1  = RadioLensfit2-MS.o galaxy_fitting.o utils.o read_catalog.o data_simulation.o distributions.o galaxy_visibilities.o evaluate_uv_grid.o likelihood.o marginalise_r.o measurement_set.o
endif 

#EXEC2 = Simulate
#OBJS2 = Simulate.o generate_catalog.o distributions.o utils.o generate_random_values.o galaxy_visibilities.o measurement_set.o data_simulation.o
 
EXEC2 = Simulate
OBJS2 = Simulate-from-catalog.o read_catalog.o distributions.o utils.o galaxy_visibilities.o measurement_set.o data_simulation.o

OBJS = RadioLensfit2-MS.o RadioLensfit-mpi.o galaxy_fitting-mpi.o galaxy-fitting.o Simulate.o utils.o generate_catalog.o read_catalog.o data_simulation.o distributions.o galaxy_visibilities.o evaluate_uv_grid.o generate_random_values.o galaxy_fitting.o likelihood.o  marginalise_r.o measurement_set.o


EXEC3 = RadioLensfit2-single   
 
OBJS3  = RadioLensfit2-single.o utils.o read_catalog.o measurement_set.o data_simulation.o distributions.o galaxy_visibilities.o  evaluate_uv_grid.o generate_random_values.o galaxy_fitting.o likelihood.o marginalise_r.o  
 
#OBJS3  = RadioLensfit2-single-noMS.o utils.o read_coordinates.o data_simulation.o distributions.o galaxy_visibilities.o  evaluate_uv_grid.o generate_random_values.o galaxy_fitting.o likelihood.o marginalise_r.o

INCL   = *.h Makefile
LIB_OPT =  -lgsl -lgslcblas -lm $(CASACORE_LIB) $(GSL_LIB)

CPPFLAGS = $(OPTIONS) $(SUP_INCL)  $(OMP)

LIBS   = $(LIB_OPT) $(OMP)

.SUFFIXES: .o .cc .cxx .cpp 

.cc.o:
	$(CC) -g -c $(CPPFLAGS) -o "$@" "$<"

.cxx.o:
	$(CC) -g -c $(CPPFLAGS) -o "$@" "$<"

.cpp.o:
	$(CC) -g -c $(CPPFLAGS) -o "$@" "$<"


RadioLensfit: $(OBJS1)
	$(CC)  $(OBJS1)  $(OPTIONS) $(LIBS) -o $(EXEC1)

$(OBJS1): $(INCL)

Simulate: $(OBJS2)
	$(CC)  $(OBJS2)  $(OPTIONS) $(LIBS) -o $(EXEC2)

$(OBJS2): $(INCL)

all: $(OBJS)
	$(CC)  $(OBJS1)  $(OPTIONS) $(LIBS) -o $(EXEC1)
	$(CC)  $(OBJS2)  $(OPTIONS) $(LIBS) -o $(EXEC2)

$(OBJS): $(INCL)



RadioLensfit2-single: $(OBJS3)
	$(CC)  $(OBJS3)  $(OPTIONS) $(LIBS) -o $(EXEC3)

$(OBJS3): $(INCL)

clean:
	rm -f $(OBJS)

realclean: clean
	rm -f $(EXEC)


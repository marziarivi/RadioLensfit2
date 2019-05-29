###############################################################
#                                                             #
###############################################################

#------------------------------------------------------------------------- 
 OMP  = -fopenmp  # OpenMP compiler switch

#OPT += -DUSE_MPI

 OPT += -DFACET  #fitting by faceting
ifeq (FACET,$(findstring FACET,$(OPT)))
 OPT += -DGRID
endif

# OPT += -DGRID

#------------------------------------------------------------------------- 

CASACORE_INC =  -I/iranet/arcsoft/casacore/casacore-2.4.1/include -I/iranet/arcsoft/casacore/casacore-2.4.1/include/casacore
CASACORE_LIB =  -L/iranet/arcsoft/casacore/casacore-2.4.1/lib -lcasa_casa -lcasa_measures -lcasa_tables -lcasa_ms -std=c++11

SUP_INCL = -I. $(CASACORE_INC)
OPTIMIZE = -O3 -g 

ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
  CC  =  mpiCC -g 
else
  CC  = g++
endif

OPTIONS = $(OPTIMIZE) $(OPT)
EXEC1 = RadioLensfit2-MS   
OBJS1  = RadioLensfit2-MS.o utils.o read_catalog.o data_simulation.o distributions.o galaxy_visibilities.o  evaluate_uv_grid.o galaxy_fitting.o likelihood.o marginalise_r.o measurement_set.o 
 
EXEC2 = Simulate
OBJS2 = Simulate.o generate_catalog.o distributions.o utils.o generate_random_values.o galaxy_visibilities.o measurement_set.o data_simulation.o
 
OBJS = RadioLensfit2-MS.o Simulate.o utils.o generate_catalog.o read_catalog.o data_simulation.o distributions.o galaxy_visibilities.o  evaluate_uv_grid.o generate_random_values.o galaxy_fitting.o likelihood.o  marginalise_r.o measurement_set.o


#OBJS  = RadioLensfit2.o utils.o read_coordinates.o generate_catalog.o data_simulation.o distributions.o galaxy_visibilities.o  evaluate_uv_grid.o generate_random_values.o galaxy_fitting.o likelihood.o marginalise_r.o 

 
INCL   = *.h Makefile
LIB_OPT = -lgsl -lgslcblas -lm $(CASACORE_LIB)

CPPFLAGS = $(OPTIONS) $(SUP_INCL)  $(OMP)

LIBS   = $(LIB_OPT) $(OMP)

.SUFFIXES: .o .cc .cxx .cpp .cu

.cc.o:
	$(CC) -c $(CPPFLAGS) -o "$@" "$<"

.cxx.o:
	$(CC) -c $(CPPFLAGS) -o "$@" "$<"

.cpp.o:
	$(CC) -c $(CPPFLAGS) -o "$@" "$<"

RadioLensfit: $(OBJS1)
	$(CC)  $(OBJS1)  $(OPTIONS) $(LIBS) -o $(EXEC1)

$(OBJS1): $(INCL)

Simulate: $(OBJS2)
	$(CC)  $(OBJS2)  $(OPTIONS) $(LIBS) -o $(EXEC2)

all: $(OBJS)
	$(CC)  $(OBJS1)  $(OPTIONS) $(LIBS) -o $(EXEC1)
	$(CC)  $(OBJS2)  $(OPTIONS) $(LIBS) -o $(EXEC2)

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS) 

realclean: clean
	rm -f $(EXEC1) $(EXEC2) 


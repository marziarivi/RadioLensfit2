###############################################################
#                                                             #
###############################################################

#-------------------------------------------------------------------------  
# OpenMP compiler switch
 OMP = -fopenmp

#OPT += -DUSE_MPI

 OPT += -DFACET  #fitting by faceting
# OPT += -DGRID

#------------------------------------------------------------------------- 

CASACORE_INC =  -I/iranet/arcsoft/casacore/casacore-2.4.1/include -I/iranet/arcsoft/casacore/casacore-2.4.1/include/casacore
CASACORE_LIB =  -L/iranet/arcsoft/casacore/casacore-2.4.1/lib -lcasa_casa -lcasa_measures -lcasa_tables -lcasa_ms -std=c++11

ifeq (FACET,$(findstring FACET,$(OPT)))
 OPT += -DGRID
endif

SUP_INCL = -I. $(CASACORE_INC)
OPTIMIZE = -O3 -g 

ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
  CC  =  mpiCC -g 
else
  CC  = g++
endif

OPTIONS = $(OPTIMIZE) $(OPT)
EXEC = RadioLensfit2-MS.x   

OBJS  = RadioLensfit2-MS.o utils.o generate_catalog.o data_simulation.o distributions.o galaxy_visibilities.o  evaluate_uv_grid.o generate_random_values.o galaxy_fitting.o likelihood.o random_gaussian.o marginalise_r.o ms_reader.o ms_utils.o
 
#OBJS  = RadioLensfit2.o utils.o read_coordinates.o generate_catalog.o data_simulation.o distributions.o galaxy_visibilities.o  evaluate_uv_grid.o generate_random_values.o galaxy_fitting.o likelihood.o random_gaussian.o marginalise_r.o 

 
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


$(EXEC): $(OBJS)
	$(CC)  $(OBJS)  $(OPTIONS) $(LIBS) -o $(EXEC)

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS) 

realclean: clean
	rm -f $(EXEC) 


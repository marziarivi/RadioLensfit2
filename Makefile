###############################################################
#                                                             #
###############################################################

#-------------------------------------------------------------------------  
# OpenMP compiler switch
OMP = -fopenmp

OPT += -DUSE_MPI
OPT += -DFACET

#default: Sersic profile, read scalelength
#OPT += GAUSSIAN            # Gaussian profile, read FWHM
#OPT += MATCH_EXP           # Sersic profile, Galsim Exponential matched scalelength, read FWHM 

#OPT += -DSCALELENGTH_ON    # use size provided by the catalog

#------------------------------------------------------------------------- 

# HOTCAT
CASACORE_INC =  -I/u/mrivi/casacore-3.3.0/include -I/u/mrivi/casacore-3.3.0/include/casacore
CASACORE_LIB =  -L/u/mrivi/casacore-3.3.0/lib -lcasa_casa -lcasa_measures -lcasa_tables -lcasa_ms -std=c++11

GSL_INC = -I/opt/cluster/gsl/2.5/gnu/9.3.0/include
GSL_LIB = -L/opt/cluster/gsl/2.5/gnu/9.3.0/lib

SUP_INCL = -I. $(CASACORE_INC) $(GSL_INC)  


# INTEL compiler
#-------------------
#OPTIMIZE = -O3 #bad results! -ip -ipo #-xAVX 
#OPTIMIZE = -O3 -ip -ipo  # using compilers/intel/15.1/133

ifeq (MPI,$(findstring MPI,$(OPT)))
#  CC  =  mpiicpc
else
#  CC  = icpc 
endif

# GNU compiler
# #-----------------
OPTIMIZE = -O3 -msse2 -mfpmath=sse -ffast-math -ftree-vectorizer-verbose=2

ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
  CC  =  mpiCC -g 
else
  CC  = g++
endif

OPTIONS = -std=gnu++11 -D_GLIBCXX_USE_CXX11_ABI=0 $(OPTIMIZE) $(OPT)

ifeq (MPI,$(findstring MPI,$(OPT)))
EXEC1 = RadioLensfit2-mpi
OBJS1  = RadioLensfit2-mpi.o data_processing.o source_extraction.o galaxy_fitting.o utils.o read_catalog.o data_simulation.o distributions.o galaxy_visibilities.o  evaluate_uv_grid.o likelihood.o marginalise_r.o measurement_set.o tukey_tapering.o
else
EXEC1 = RadioLensfit2
OBJS1  = RadioLensfit2.o data_processing.o source_extraction.o galaxy_fitting.o utils.o read_catalog.o data_simulation.o distributions.o galaxy_visibilities.o evaluate_uv_grid.o likelihood.o marginalise_r.o measurement_set.o tukey_tapering.o
endif 

EXEC2 = Catalog
OBJS2 = Generate_catalog.o distributions.o utils.o generate_random_values.o

EXEC3 = Simulate
OBJS3 = Simulate-from-catalog.o read_catalog.o distributions.o utils.o galaxy_visibilities.o measurement_set.o data_simulation.o

OBJS = RadioLensfit2.o RadioLensfit-mpi.o Simulate-from-catalog.o data_processing.o source_extraction.o galaxy_fitting-mpi.o galaxy-fitting.o Simulate.o utils.o generate_catalog.o read_catalog.o 
data_simulation.o distributions.o galaxy_visibilities.o evaluate_uv_grid.o generate_random_values.o galaxy_fitting.o likelihood.o  marginalise_r.o measurement_set.o tukey_tapering.o
EXEC = RadioLensfit2-mpi RadioLensfit2 Catalog Simulate

#Generate galaxy catalog, simulate visibilities one at a time, and measure its ellipticity
EXEC4 = RadioLensfit2-single    
OBJS4  = RadioLensfit2-single.o utils.o read_catalog.o measurement_set.o data_simulation.o distributions.o galaxy_visibilities.o  evaluate_uv_grid.o generate_random_values.o galaxy_fitting.o likelihood.o marginalise_r.o tukey_tapering.o  
 
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

Catalog: $(OBJS2)
	$(CC)  $(OBJS2)  $(OPTIONS) $(LIBS) -o $(EXEC2)

$(OBJS2): $(INCL)

Simulate: $(OBJS3)
	$(CC)  $(OBJS3)  $(OPTIONS) $(LIBS) -o $(EXEC3)

$(OBJS3): $(INCL)

all: $(OBJS)
	$(CC)  $(OBJS1)  $(OPTIONS) $(LIBS) -o $(EXEC1)
	$(CC)  $(OBJS2)  $(OPTIONS) $(LIBS) -o $(EXEC2)
        $(CC)  $(OBJS3)  $(OPTIONS) $(LIBS) -o $(EXEC3)

$(OBJS): $(INCL)



RadioLensfit2-single: $(OBJS4)
	$(CC)  $(OBJS4)  $(OPTIONS) $(LIBS) -o $(EXEC4)

$(OBJS4): $(INCL)

clean:
	rm -f $(OBJS)

realclean: clean
	rm -f $(EXEC)


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

ifeq (FACET,$(findstring FACET,$(OPT)))
 OPT += -DGRID
endif

SUP_INCL = -I. 
LIB_OPT = 
OPTIMIZE = -O3 -g 

ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
  CC  =  mpiCC -g 
else
  CC  = g++-4.9
endif

OPTIONS = $(OPTIMIZE) $(OPT)
EXEC = RadioLensfit2-single2.x   
#EXEC = GalaxyCatalog.x
 
OBJS  = RadioLensfit2-single2.o utils.o read_coordinates.o generate_catalog.o data_simulation.o distributions.o galaxy_visibilities.o  evaluate_uv_grid.o generate_random_values.o galaxy_fitting.o likelihood.o random_gaussian.o marginalise_r.o evaluate_image_lm_grid.o

#OBJS = GalaxyCatalog.o generate_catalog.o distributions.o generate_random_values.o random_gaussian.o

 
INCL   = *.h Makefile
LIB_OPT = -lgsl -lgslcblas -lm

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


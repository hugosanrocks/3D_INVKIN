#include ../../Makefile.inc
include ../Makefile.inc

CXX = $(FC)
FLAGS = $(OPTF)

##SOURCES 
COM = ../COMMON/src/
KER = ./kernel/src/
LIB = ./lib/
CXX_SOURCES_kernel = $(KER)TRN.f90 $(KER)init_TRN.f90 $(KER)descent_TRN.f90 $(KER)forcing_term_TRN.f90 $(KER)print_info_TRN.f90 $(KER)finalize_TRN.f90
CXX_SOURCES_common = $(COM)normL2.f90 $(COM)scalL2.f90 $(COM)std_linesearch.f90 $(COM)std_test_conv.f90 $(COM)project.f90

CXX_SOURCES= $(CXX_SOURCES_common) $(CXX_SOURCES_kernel)  
CXX_OBJECTS = $(CXX_SOURCES:.f90=.o) 

INCLUDE = -I ../COMMON/include

%.o : %.f90
	$(CXX) $(FLAGS) $(INCLUDE) -c -o $@ $<

donothing : # for giving the information when typing make
	@echo please type -make lib- for compiling the library

lib : $(CXX_OBJECTS) 
	$(AR) $(ARFUNCT) $(LIB)libTRN.a $(CXX_OBJECTS)
	cp $(LIB)libTRN.a ../lib

clean :	
	rm -f $(KER)*.o $(KER)*~ $(KER)*backup $(KER)*.mod $(KER)*genmod* 
	rm -f $(COM)*.o $(COM)*~ $(COM)*backup $(COM)*.mod $(COM)*genmod* 
	rm -f *.mod *genmod* 
	rm ./lib/*.a	
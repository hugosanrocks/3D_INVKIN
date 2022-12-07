#include ../../Makefile.inc
include ../Makefile.inc

CXX = $(FC)
FLAGS = $(OPTF)

##SOURCES 
COM = ../COMMON/src/
KER_TRN = ../TRN/kernel/src/
KER = ./kernel/src/

LIB = ./lib/
CXX_SOURCES_kernel = $(KER)PTRN.f90 $(KER)descent_PTRN.f90 $(KER)init_PTRN.f90 $(KER)print_info_PTRN.f90 $(KER_TRN)forcing_term_TRN.f90 $(KER)finalize_PTRN.f90
CXX_SOURCES_common = $(COM)normL2.f90 $(COM)scalL2.f90 $(COM)std_linesearch.f90 $(COM)std_test_conv.f90 $(COM)project.f90

CXX_SOURCES= $(CXX_SOURCES_common) $(CXX_SOURCES_kernel)  
CXX_OBJECTS = $(CXX_SOURCES:.f90=.o) 

INCLUDE = -I ../COMMON/include

%.o : %.f90
	$(CXX) $(FLAGS) $(INCLUDE) -c -o $@ $<

donothing : # for giving the information when typing make
	@echo please type -make lib- for compiling the library

lib : $(CXX_OBJECTS) 
	$(AR) $(ARFUNCT) $(LIB)libPTRN.a $(CXX_OBJECTS)
	cp $(LIB)libPTRN.a ../lib

clean :	
	rm -f $(KER)*.o $(KER)*~ $(KER)*backup $(KER)*.mod $(KER)*genmod* 
	rm -f $(COM)*.o $(COM)*~ $(COM)*backup $(COM)*.mod $(COM)*genmod* 
	rm -f *.mod *genmod* 
	rm ./lib/*.a	
include ../../Makefile.inc

## SOURCES 
COM = ../COMMON/src
KER = ./kernel/src
KER_LB = ../LBFGS/kernel/src
LIB = ./lib
CXX_SOURCES_kernel = $(KER)/PLBFGS.f90 $(KER)/init_PLBFGS.f90 $(KER)/descent_PLBFGS.f90 $(KER_LB)/save_LBFGS.f90 $(KER_LB)/update_LBFGS.f90 $(KER)/finalize_PLBFGS.f90
CXX_SOURCES_common = $(COM)/normL2.f90 $(COM)/scalL2.f90 $(COM)/std_linesearch.f90 $(COM)/std_test_conv.f90 $(COM)/print_info.f90 $(COM)/project.f90

CXX_SOURCES= $(CXX_SOURCES_common) $(CXX_SOURCES_kernel)  
CXX_OBJECTS = $(CXX_SOURCES:.f90=.o) 

## TEST
TEST_DIR = ./test
TEST_SRC = $(TEST_DIR)/src
TEST_BIN = $(TEST_DIR)/bin
TEST_RUN = $(TEST_DIR)/run
TEST_COMMON_SRC = ../COMMON/test
TEST_SOURCES = $(TEST_COMMON_SRC)/rosenbrock.f90 $(TEST_SRC)/test_PLBFGS.f90
TEST_OBJECTS = $(TEST_SOURCES:.f90=.o) 
TEST_TARGET = $(TEST_BIN)/test.bin
TEST_LIB = -L../lib -lSEISCOPE_OPTIM
## Common to compil and test
INCLUDE = -I ../COMMON/include

%.o : %.f90
	$(FC) $(FFLAGS) $(INCLUDE) -c -o $@ $<

all : lib test 

lib : $(CXX_OBJECTS) 
	$(AR) $(ARFUNCT) $(LIB)/libPLBFGS.a $(CXX_OBJECTS)
	cp $(LIB)/libPLBFGS.a ../lib

clean :	
	$(RM) $(KER)/*.o $(KER)/*~ $(KER)/*backup $(KER)/*.mod $(KER)/*genmod* 
	$(RM) $(COM)/*.o $(COM)/*~ $(COM)/*backup $(COM)/*.mod $(COM)/*genmod* 
	$(RM) *.mod *genmod* $(LIB)/*.a

test : $(TEST_TARGET)
	cd $(TEST_RUN) ;  ../bin/test.bin

$(TEST_TARGET) : $(TEST_OBJECTS) 
	$(FC) $(FFLAGS) $(TEST_OBJECTS) -o $(TEST_TARGET) $(TEST_LIB)

testclean :
	$(RM) $(TEST_TARGET) $(TEST_RUN)/*dat $(TEST_SRC)/*.o $(TEST_COMMON_SRC)/*.o

distclean : clean testclean

include ../../../../Makefile.inc

##SOURCES 

COMMON_SRC = ../../../COMMON/test/
CXX_SOURCES = $(COMMON_SRC)rosenbrock.f90 test_PSTD.f90 
CXX_OBJECTS = $(CXX_SOURCES:.f90=.o) 

INCLUDE = -I ../../../COMMON/include
TARGET = ../bin/test.bin 

LIB =  -L ../../../lib -lSEISCOPE_OPTIM

$(TARGET) : $(CXX_OBJECTS) 
	$(FC) $(FFLAGS) $(CXX_OBJECTS) -o $(TARGET) $(LIB)

%.o : %.f90
	$(FC) $(FFLAGS) $(INCLUDE) -c -o $@ $<

clean :
	$(RM) *.o *~ *backup *.mod *genmod*
	$(RM) $(COMMON_SRC)*.o $(COMMON_SRC)*~ $(COMMON_SRC)*backup $(COMMON_SRC)*.mod $(COMMON_SRC)*genmod*


test:
	$(RM) *.o ../run/*.dat $(TARGET) ; $(MAKE) ; cd ../run ; $(TARGET)


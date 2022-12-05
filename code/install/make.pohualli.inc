#
# libraries
#

# Library variables
#LTOOLS_BOX = ../../TOOLS_BOX
LTOOLS_BOX =~/svn/TOOLS_BOX/trunk


LFFT = ~/fftw3/lib
FFTW3 = $(LFFT) -lfftw3f


##OPTIMIZATION LIBRARIES
LPSTD = $(LTOOLS_BOX)/OPTIMIZATION/lib -lPSTD
LPLBFGS = $(LTOOLS_BOX)/OPTIMIZATION/lib -lPLBFGS
LLBFGS = $(LTOOLS_BOX)/OPTIMIZATION/lib -lLBFGS
LPNLCG = $(LTOOLS_BOX)/OPTIMIZATION/lib -lPNLCG

IOPTIM  = $(LTOOLS_BOX)/OPTIMIZATION/COMMON/include

#
# compiler options
#
#FOPTION= -O3 -cpp
FOPTION=-O3 -cpp 
##-check bounds -check all -traceback
##-fbounds-check -DDOUBLE_PRECISION -g -fcheck=all -Wall -fbacktrace
CC_OPTION= -O3

#
# compilers
#
#F77=gfortran
#F90=gfortran
#CC=gcc
F77=ifort
F90=ifort
CC=icc



# if DOUBLE PRECISION is desired for FFTW3
# LIBS2 = -lfftw3 -lm -lc
# and change all the sfft subroutines call in forward folder

# libs
LIBS = -L $(LPSTD) -L $(LPLBFGS) -L $(LLBFGS) -L $(LPNLCG)
LIBS2 = -L $(FFTW3) -lm -lc


# libs mkl
LMKL_DIR= /opt/intel/compilers_and_libraries_2016.0.109/linux/mkl/lib/intel64
LIBS3= -L$(LMKL_DIR) -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -openmp













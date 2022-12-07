#
# libraries
#

# Library variables for inversion (SEISCOPE TOOLS)
LTOOLS_BOX = /home/sanchezh/Tools/3D_INVKIN/code/other_tools/TOOLS_BOX/trunk

# FFT Dir and libraries
fftw_DIR = /home/sanchezh/Tools/3D_INVKIN/code/other_tools/fftw-3.3.10
LFFT = $(fftw_DIR)/fftw-3.3.10
FFTW3 = $(LFFT) -lfftw3f


##OPTIMIZATION LIBRARIES INVERSION

LPSTD = $(LTOOLS_BOX)/OPTIMIZATION/lib -lPSTD
LPLBFGS = $(LTOOLS_BOX)/OPTIMIZATION/lib -lPLBFGS
LLBFGS = $(LTOOLS_BOX)/OPTIMIZATION/lib -lLBFGS
LPNLCG = $(LTOOLS_BOX)/OPTIMIZATION/lib -lPNLCG

IOPTIM  = $(LTOOLS_BOX)/OPTIMIZATION/COMMON/include

#
# compiler options
#
FOPTION= -O3 -cpp -assume byterecl
#FOPTION=-O3 -cpp -msse2 -check bounds -check all -traceback

##DEBUG OPTIONS
##FOPTION=-cpp -C -g -DEBUG -check -assume byterecl -check all -gen_interfaces -traceback -fpe0 -warn all
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



# libs fft and inversion
LIBS = -L $(LPSTD) -L $(LPLBFGS) -L $(LLBFGS) -L $(LPNLCG)
LIBS2 = -L $(FFTW3) -lm -lc

# libs mkl
LMKL_DIR=/home/sanchezh/intel/oneapi/mkl/latest/lib/intel64/
LIBS3= -L$(LMKL_DIR) -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -qopenmp


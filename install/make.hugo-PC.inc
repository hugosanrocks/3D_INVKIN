#
# Path to libraries
#

# Optimization tool box
LTOOLS_BOX =~/svn/TOOLS_BOX/trunk

# OPTIMIZATION LIBRARIES
LPSTD = $(LTOOLS_BOX)/OPTIMIZATION/lib -lPSTD
LPLBFGS = $(LTOOLS_BOX)/OPTIMIZATION/lib -lPLBFGS
LLBFGS = $(LTOOLS_BOX)/OPTIMIZATION/lib -lLBFGS
LPNLCG = $(LTOOLS_BOX)/OPTIMIZATION/lib -lPNLCG

# Common file
IOPTIM  = $(LTOOLS_BOX)/OPTIMIZATION/COMMON/include

####################
# compiler options
####################
FOPTION=-O3 -cpp -assume byterecl 
##-msse2 -check bounds -check all -traceback
##DEBUG
##FOPTION=-cpp -C -g -DEBUG -check -assume byterecl -check all -gen_interfaces -traceback -fpe0 -warn all 
#for debugging with ifort

#for debugging with gfortran 
##-fbounds-check -DDOUBLE_PRECISION -g -fcheck=all -Wall -fbacktrace

#C options
CC_OPTION= -O3

##################
# compilers
##################
#CODES DOES NOT WORK
#WITH GFORTRAN AS
#SOME MATRIX MULTIPLICATIONS
#REQUIRE INTEL MKL
#LIBRARIES

#gfortran
#F77=gfortran
#F90=gfortran
#CC=gcc

#ifort
F77=ifort
F90=ifort
CC=icc


# All optimization libraries
LIBS = -L $(LPSTD) -L $(LPLBFGS) -L $(LLBFGS) -L $(LPNLCG)
# fft library
fftw_DIR=~/fftw3
LIBS2 = -lfftw3f -lm -lc

# library mkl for bigmatrix multiplication
LMKL_DIR=/opt/intel/mkl/lib/intel64/
LIBS3= -L$(LMKL_DIR) -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -qopenmp

# Covariance tools_box
##COV_TOOL_BOX_DIR=~/svn/COVARIANCE_TOOLBOX/trunk
##LIB_COV_TOOL_BOX= -L $(COV_TOOL_BOX_DIR)/lib -lLAP_COV
##INC_CONV_TOOLS_BOX = -I $(COV_TOOL_BOX_DIR)/include -I



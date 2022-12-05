#
# libraries
#

# Library variables
LTOOLS_BOX = /data/cycle/sanchezh/TOOLS_BOX/trunk

fftw_DIR = /data/cycle/sanchezh/3D_INVKIN/branches/predict_prior
LFFT = /data/cycle/sanchezh/fftw-3.3.10/lib
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
FOPTION= -O3 -cpp -assume byterecl
#FOPTION=-O3 -cpp -msse2 -check bounds -check all -traceback
##DEBUG
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



# if DOUBLE PRECISION is desired for FFTW3
# LIBS2 = -lfftw3 -lm -lc
# and change all the sfft subroutines call in forward folder

# libs
LIBS = -L $(LPSTD) -L $(LPLBFGS) -L $(LLBFGS) -L $(LPNLCG)
LIBS2 = -L $(FFTW3) -lm -lc


# libs mkl
LMKL_DIR= /soft/intel-2018/compilers_and_libraries/linux/mkl/lib/intel64
LIBS3= -L$(LMKL_DIR) -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -qopenmp

# Covariance tools_box
#COV_TOOL_BOX_DIR=~/svn/COVARIANCE_TOOLBOX/trunk
#LIB_COV_TOOL_BOX= -L $(COV_TOOL_BOX_DIR)/lib -lLAP_COV
#INC_CONV_TOOLS_BOX = -I $(COV_TOOL_BOX_DIR)/include -I












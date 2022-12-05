#
# libraries
#

INC = 
LIB = 

#
# compiler options
#
FOPTION = -O3 -cpp
# FOPTION = -C -g -cpp -fbacktrace -fbounds-check 

#
# compilers
#
FC = mpif90
CC = mpicc
# FC = gfortran
# CC = gcc

LIBS = -L $(LIB)

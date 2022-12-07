#Options for all Makefile defined here

# ---------------------------------------------------------- #
# --- Intel compiler by default (make -j 16)
# ---------------------------------------------------------- #
FC = ifort
CC =  icc
# ---------------------------------------------------------- #
# --- Preprocessing, opt. &  flags
# ---------------------------------------------------------- #
OPTPRE = -qopenmp -fpp #-DTIME_PROFILING
OPTF = -O3 -xHost -inline $(OPTPRE)
FFLAGS = $(OPTF) -assume byterecl -heap-arrays 50
OPTC =  -O3
CFLAGS = $(OPTC)
ifeq ($(DBG),1) # Debug
  OPTF = -C -g -warn -DEBUG -check -assume byterecl -heap-arrays 50 -check all -gen_interfaces -traceback -check bounds -fpe0 -warn all $(OPTPRE) # += -ipo
  FFLAGS = $(OPTF)
  OPTC = -g
  CFLAGS = $(OPTC)
endif
# ---------------------------------------------------------- #
# --- variables for the libraries
# ---------------------------------------------------------- #
AR = xiar
ARFUNCT = cruvfs


# ---------------------------------------------------------- #
# --- Overload some variables for the case of gnu compilers
# --- (make GNU=1 -j 16)
# ---------------------------------------------------------- #
ifeq ($(GNU),1)
  FC = gfortran
  CC = gcc
  # --- Preprocessing, opt. &  flags
  OPTPRE = -fopenmp -cpp -ffree-line-length-0 #-DTIME_PROFILING
  OPTF = -O3 $(OPTPRE)
  FFLAGS = $(OPTF)
  OPTC = -O3 $(OPTPRE)
  CFLAGS = $(OPTC)
  ifeq ($(DBG),1) # Debug -ffpe-trap=invalid,overflow,zero
    OPTF = -g -fbounds-check -ffpe-trap=invalid,overflow -Wtabs $(OPTPRE)
    FFLAGS = $(OPTF) -pedantic
    OPTC = -g $(OPTPRE)
    CFLAGS = $(OPTC)
   endif
   # --- variables for the libraries
   AR = ar
   UNAME_S := $(shell uname -s)
   ifeq ($(UNAME_S),Darwin)
     ARFUNCT = cruvTs
   endif
endif


# ---------------------------------------------------------- #
# --- overload some variables for the case of IBM compilers
# --- (make XL=1 -j 16)
# ---------------------------------------------------------- #
ifeq ($(XL),1)
  FC = xlf90_r
  CC = xlc_r
  # --- Preprocessing, opt. &  flags
  OPTPRE = -qsuffix=cpp=f90 #-WF,-DTIME_PROFILING
  OPTF = -O3 -qstrict -qsmp=omp -shared -qalign=struct=natural $(OPTPRE)
  FFLAGS = $(OPTF)
  OPTC = -O3 $(OPTPRE)
  CFLAGS = $(OPTC)
  ifeq ($(DBG),1) # Debug
    OPTF = -g $(OPTPRE)
    FFLAGS = $(OPTF) -qsmp=omp -shared -qalign=struct=natural
    OPTC = -g $(OPTPRE)
    CFLAGS = $(OPTC) -qsmp=omp -shared -qalign=struct=natural
  endif
  # --- variables for the libraries
  AR = ar
endif

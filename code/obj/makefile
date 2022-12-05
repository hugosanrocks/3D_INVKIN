
########################################
#  INV3DKIN Makefile
#
#  Requires IFORT to compile
#  gfortran is not yet availabel
########################################

HOST=$(shell hostname)


#opt=1   # set generic value
#for jean-PC
ifeq ($(findstring HP-115,$(strip $(HOST))), $(strip HP-115))
include ../install/make.jean-PC.inc
opt=2
endif
#for hugo-PC
ifeq ($(findstring hugo,$(HOST)), hugo)
include ../install/make.hugo-PC.inc
opt=3
endif
#for POHULLI cluster
ifeq ($(findstring pohualli,$(HOST)), pohualli)
include ../install/make.pohualli.inc
opt=4
endif
#for IST-OAR cluster
ifeq ($(findstring ist,$(HOST)), ist)
include ../install/make.istoar.inc
opt=5
endif
#for froggy cluster
ifeq ($(findstring frog,$(HOST)), frog)
include ../install/make.froggy.inc
opt=6
endif
#generic option. Modify your ../install/make.generic.inc
#if necessary
ifeq ($(opt), 1)
include ../install/make.generic.inc
endif
#for IST-OAR cluster
ifeq ($(findstring ist-159,$(strip $(HOST))), $(strip ist-159))
include ../install/make.sanchezh.inc
opt=7
endif


#options of make
all : INV3DKIN PREPROCESS FORWARD SENSITIVITY HESSIAN


#directories
BASE_DIR        =      ..
INV_BASE_DIR    =      $(BASE_DIR)
PREPROC_DIR     =      $(INV_BASE_DIR)/preprocess
INV_GENERAL_DIR =      $(INV_BASE_DIR)/general
FILTER_DIR      =      $(INV_BASE_DIR)/filtro/test
FWI_INV_DIR     =      $(INV_BASE_DIR)/inversion
COV_INV_DIR     =      $(INV_BASE_DIR)/covariance
FOR_INV_DIR     =      $(INV_BASE_DIR)/forward
ADJ_INV_DIR     =      $(INV_BASE_DIR)/adjoint
SENSITI_DIR     =      $(INV_BASE_DIR)/sensitivity
HESSIAN_DIR     =      $(INV_BASE_DIR)/hessian
INV_INC_DIR     =      $(INV_BASE_DIR)/include
INV_EXE_DIR     =      $(INV_BASE_DIR)/bin
PRIOR___DIR     =      $(INV_BASE_DIR)/prior
DELAUN__DIR     =      $(INV_BASE_DIR)/delaunay

VPATH = $(BASE_DIR):$(INV_BASE_DIR):$(HESSIAN_DIR):$(PREPROC_DIR):$(FILTER_DIR):$(INV_GENERAL_DIR):$(INV_EXE_DIR):$(COV_INV_DIR):$(FWI_INV_DIR):$(FOR_INV_DIR):$(ADJ_INV_DIR):$(PRIOR___DIR):$(SENSITI_DIR):$(DELAUN__DIR)

FFLAGS    = -I. -I $(INV_INC_DIR) -I $(IOPTIM) -I $(fftw_DIR)/include $(FOPTION)
FFLAGS2   = -qopenmp
#FFLAGS2   = -openmp     The Hugo flag is the correct one: should stay like that

#objects

COV_OBJ     =  covariance.o exp_covar.o edge.o covedge.o modeledge.o timemask.o covtimer.o modeltimer.o laplacian.o covtiko.o modeltiko.o modelpri.o covpri.o

FILTER_OBJ  =  init_butterworth.o filtfilt_butterworth.o

CAUSAL_OBJ  =  synwindow.o

GENERAL_OBJ = lap_filter.o model_norm.o orderinv.o linreg.o conv.o ortogonal.o rotate.o rakeout.o modeldecomp.o inv_mat.o rakeuni.o rakeangle.o normal_vector.o

PREPROC_OBJ = preprocess.o readinfop.o mul_src.o stress_axi.o interpolation.o stress_dwn.o interp_axi.o interp_dwn.o traction.o normvec.o tract_fft.o medium.o

OPTIM_OBJ = laplace_matrix.o fcost.o fwi_option.o fwi_pstd.o fwi_lbfgs.o fwi_plbfgs.o fwi_pnlcg.o invoption.o prediction.o detect.o fwi_precondition.o depth_preco.o

INIT_OBJ = initialize.o monitor_alloc.o modelcheck.o initwin.o

MAIN_OBJ  = main.o

SENSITIVE_OBJ = lap_filter.o readinfo.o normal_vector.o coortrans.o write_syn.o initialize.o monitor_alloc.o residual.o adjtraction.o initwin.o syn_velocity.o conv.o readgrad.o grad_time_interp.o modeldecomp.o orderinv.o synwindow.o writegrad.o write_residual.o rotate.o

DELAUN_OBJ  =  interplib.o triangles.o delaunay_mod.o interpolation_mod.o

FORWARD_OBJ = laplace_filter.o readinfo.o coortrans.o forward.o time_interp.o syn_velocity.o write_syn.o progress.o

ADJOINT_OBJ = adjoint.o residual.o write_residual.o adjtraction.o grad_time_interp.o read_obs2.o smooth_grad.o

MODEL_GRAD_OBJ = writemodel.o newgrad.o readgrad.o

HESS_OBJ = $(DELAUN_OBJ) $(FILTER_OBJ) $(ADJOINT_OBJ) $(INIT_OBJ) lap_filter.o normal_vector.o orderinv.o conv.o rotate.o modeldecomp.o $(FORWARD_OBJ) writemodel.o readgrad.o synwindow.o writegrad.o

PRIOR_OBJ = prior_past.o

ALL_OBJ = $(DELAUN_OBJ) $(INIT_OBJ) $(FILTER_OBJ) $(PREPROC_OBJ) $(GENERAL_OBJ) $(OPTIM_OBJ) $(MAIN_OBJ) $(FORWARD_OBJ) $(MODEL_GRAD_OBJ) $(ADJOINT_OBJ) $(COV_OBJ) $(CAUSAL_OBJ) $(PRIOR_OBJ)

#make options:

FOR_OBJ = $(FILTER_OBJ) $(DELAUN_OBJ) $(PREPROC_OBJ) $(CAUSAL_OBJ) residual.o fcost.o mainforward.o write_residual.o writemodel.o $(INIT_OBJ) $(FORWARD_OBJ) $(GENERAL_OBJ)

SENSITIVITY_OBJ = $(DELAUN_OBJ) $(FILTER_OBJ) $(SENSITIVE_OBJ) smooth_grad.o sensitive.o

HESSIAN_OBJ = $(HESS_OBJ) main_hessian.o write_hessian.o

PRE_OBJ = $(FILTER_OBJ) $(PREPROC_OBJ) detect_focal.o linreg.o main_preprocess.o rotate.o

.SUFFIXES: .f90 .o
.f90.o:
	$(F90) $(FFLAGS2) -c  $(FFLAGS) $<

.SUFFIXES: .c .o
.c.o:
	$(CC) -c $(C_OPTION) $(FFLAGS) $<

PREPROCESS : $(PRE_OBJ)
	$(F90) $(FFLAGS) -o $(INV_EXE_DIR)/$@ $(PRE_OBJ) $(LIBS2) $(FFLAGS2)

FORWARD : $(FOR_OBJ)
	$(F90) $(FFLAGS) -o $(INV_EXE_DIR)/$@ $(FOR_OBJ) $(LIBS2) $(LIBS) $(LIBS3) $(FFLAGS2)

INV3DKIN : $(ALL_OBJ)
	$(F90) $(FFLAGS) -o $(INV_EXE_DIR)/$@ $(ALL_OBJ) $(LIBS) $(LIBS2) $(LIBS3)

HESSIAN : $(HESSIAN_OBJ)
	$(F90) $(FFLAGS) -o $(INV_EXE_DIR)/$@ $(HESSIAN_OBJ) $(LIBS) $(LIBS2) $(LIBS3)

SENSITIVITY : $(SENSITIVITY_OBJ)
	$(F90) $(FFLAGS) -o $(INV_EXE_DIR)/$@ $(SENSITIVITY_OBJ) $(LIBS) $(LIBS2) $(LIBS3)

objs : $(ALL_OBJ)

$(ALL_OBJ): $(INV_INC_DIR)/precis.h $(INV_INC_DIR)/green.h $(fftw_DIR)/include/fftw3.f

clean:
	rm -f *.o *.f90 *~ ../bin/* *.mod


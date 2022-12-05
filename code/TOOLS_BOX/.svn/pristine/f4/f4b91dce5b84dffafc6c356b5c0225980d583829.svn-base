!########################################################
!# DATA STRUCTURE optim BUILT FOR REVERSE COMMUNICATION #
!#       IN ALL THE ALGORITHM OF THE TOOLBOX            #
!########################################################

TYPE optim_type
                                                                      
SEQUENCE

!#### DEBUG OPTION FOR USER
!# when debug is false: no information printed
!# when debug is true: additional information printed
LOGICAL :: debug

!#### BOUND CONSTRAINTS OPTION FOR USER
!# when bound is 0: no bound constraints
!# when bound is 1: bound constraints activated
integer :: bound
real :: threshold !#tolerance on bound constraints satisfaction

!#### PRINTING FLAG FOR MPI APPLICATION
INTEGER :: print_flag

!#### LINESEARCH PARAMETERS
LOGICAL :: first_ls
CHARACTER*8 :: task
INTEGER :: nls_max,cpt_ls,nfwd_pb,cpt_iter,niter_max
REAL    :: f0,fk,conv
REAL    :: m1,m2,mult_factor,alpha_L,alpha_R,alpha
REAL    :: q0, q
REAL,allocatable,dimension(:) :: ub,lb

!#### STEEPEST DESCENT PARAMETERS
REAL,allocatable,dimension(:) :: xk,grad 

!#### NONLINEAR CONJUGATE GRADIENT PARAMETERS
REAL,allocatable,dimension(:) :: grad_prev,descent,descent_prev

!### LBFGS PARAMETERS 
INTEGER :: cpt_lbfgs,l
REAL,allocatable,dimension(:,:) :: sk,yk

!### PLBFGS PARAMETERS
REAL,allocatable,dimension(:) :: q_plb,alpha_plb,rho_plb

!### TRUNCATED NEWTON PARAMETERS
LOGICAL :: conv_CG
CHARACTER*4 :: CG_phase,comm
INTEGER :: cpt_iter_CG,niter_max_CG,nhess
REAL :: qk_CG,qkm1_CG,hessian_term,eta,norm_grad,norm_grad_m1,norm_residual
REAL,dimension(:),allocatable :: residual,d,Hd,eisenvect

!### PRECONDITIONED TRUNCATED NEWTON PARAMETERS
REAL :: dHd,res_scal_respreco,alpha_CG
REAL,dimension(:),allocatable :: residual_preco



END TYPE optim_type

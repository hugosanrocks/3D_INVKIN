       TYPE mesh

       SEQUENCE

       CHARACTER (LEN=4) dat, out                     	! Working directories
       REAL*4, DIMENSION(:,:), POINTER ::  slipmod      ! Slip rate modulus
       
       !Arrays: Synthetics, adjoint traction, residuals
       ! observations, slip-rate functions
       REAL*4, DIMENSION(:,:),POINTER :: syn, synt, res
       REAL*4, DIMENSION(:,:),POINTER :: obs, slip

       !Focal mechanism (local variables)
       !hypocenter location(x,y,z), rupture velocity
       REAL*4 stk, dip, rak, hypo(3), vrup
       REAL, DIMENSION(:), POINTER :: stk_i, dip_i

       !Length along strike & dip, lengths of subfaults
       REAL*4 lenstk, lendip, stk_s, dip_s

       !Traction vector
       REAL*4 traction(3)
 
       !Arrays used for local unitary vectors:
       !normal, 3D slip, 2D slip, stike, dip, vectors 
       REAL*4, DIMENSION(:,:),POINTER :: vnorm, vslip, vslip2     ! Vectors along the fault
       REAL*4, DIMENSION(:,:),POINTER :: vstk, vdip               ! Vectors along the fault
       REAL*4, DIMENSION(:,:,:),POINTER :: slipm                  ! Vectors along the fault

       !Coordinates of subfaults (x,y,z)
       REAL*4, DIMENSION(:,:),POINTER :: fault, distcorr                    ! Position of each subfault

       !Time variables:
       !simt: final time of Green's function files
       !simdt: dt of Green's functions
       !slipt: final time of slip-rate functions
       !slipdt: dt of slip-rate functions
       REAL*4 simt, simdt, slipt, slipdt, intdt

       !Array to compute cost and matrices of penalty
       REAL*4, DIMENSION(:),POINTER :: cost, diag, rup, diag2             	  ! Cost function

       !Final cost (costa), model cost (costm)
       !hyperparameters: lam1, lam2, lam3, lam4
       !lam1: Rupture time penalty
       !lam2: Edge slip penalty
       !lam3: Tikhonov term
       !lam4: Prior model term
       REAL*4 costa, costm, costdini, costd, conv, lam1, lam2, lam3, lam4, lam5

       !seismic moment
       !Variables to remove soon rt, ot, syn_sec
       REAL*4 moment, rt, ot, syn_sec                     	  
       !seismic moment, shear modulus, rise time, origin time

       !Num stations, Num components, Num subfaults
       INTEGER*4 nsta, ncomp, msub, lensyn, lensynf, lenobs ! Number of: stations, components and subfaults, length of traces
       INTEGER*8 stcomp, sta_i, comp_i, mjump            	  ! Counter on stations, components and staXcomp, jump inside files

       INTEGER*4 simsam, slipsam, interp_i, ntint, trac_i, mod_i
       INTEGER*4 interpadj_i      	  			  ! Number of samples of: simulation, slip, forward interpolation, adjoint records
       INTEGER*4 delays, iter                              	  ! Number of samples of delay (due to origin time), number of  iteration
       INTEGER*4 modelsize, modelsize2, sizeint                      	  ! number of samples in the model to optimize


       INTEGER*4, DIMENSION(:), POINTER :: synsam, nodes, samp4win  ! rupture sample time
       REAL*4, DIMENSION(:,:), POINTER :: twin

       CHARACTER (LEN=3) :: sta, sub                      	  ! Variables to write file names
       CHARACTER (LEN=1) :: comp                         	  ! Variables to write file names
       INTEGER :: iter_i

       REAL*4, DIMENSION(:,:), POINTER :: slipr, tracadj, tottrac 
       REAL*4, DIMENSION(:,:), POINTER :: tractionvec, tracint 
								  ! Arrays to save current sliprate and gradient

       REAL*4, DIMENSION(:), POINTER :: model,modelp,model2
       REAL*4, DIMENSION(:), POINTER :: model2p,grad, grad2
       REAL*4, DIMENSION(:), POINTER :: gradad, rtimes            ! 1D arrays used by TOOLBOX to optmize
       REAL*4, DIMENSION(:,:), POINTER :: slipr2, tracad          ! 1D arrays used by TOOLBOX to optmize
       REAL*4, DIMENSION(:), POINTER :: modelint
       INTEGER*4, DIMENSION(:), POINTER :: rsamp                  ! rupture sample time

       INTEGER*4, DIMENSION(:,:), POINTER :: samwin               ! samples of time window


       !DOUBLE COMPLEX, DIMENSION (:,:), POINTER :: tracf, slipf  ! Arrays for fft conversion (slip and traction)
       COMPLEX, DIMENSION (:,:), POINTER :: tracf, slipf, resif
       INTEGER tfft_i                                             ! Counter used by fftw3

       INTEGER debug_i, stress_opt, for_opt, optm, mext, mxint    ! Debug options, origin of stress files, forward option
       LOGICAL debug         


       REAL*4, DIMENSION(:,:),POINTER :: cd, cm, ce, ct, la       ! Arrays to covariance matrices
       INTEGER weig
       REAL*4, DIMENSION(:),POINTER :: tseries                    ! Subfault's positions (x,y,z)
       INTEGER wininv, dowin, synwin                                             ! Number of the time window to forward
       INTEGER*4, DIMENSION(:), POINTER :: idsub, map
       INTEGER*4, DIMENSION(:,:), POINTER :: win
       INTEGER nsstk, nsdip
       CHARACTER*4 :: flag_d

       REAL*4, DIMENSION(:),POINTER :: model1, grad1, gradad1    ! Model, gradient, additional gradient rake fixed

       REAL*4, DIMENSION(:),POINTER :: target                  ! Array for model norm

       INTEGER modelsize1
       INTEGER rake_opt, hess_opt, hess_cont, lap_opt
       REAL*4, DIMENSION(:,:),POINTER :: precon                  ! Matrix to save H^-1
       REAL*4, DIMENSION(:),POINTER :: q_dw, q_up
       INTEGER coarse1d                                          ! Size of H^-1
       REAL*4 ub, lb                                             ! upper and lower bound
       INTEGER dir_n, dir_i, subf_i                              ! Number of directions, counter direction
       INTEGER, DIMENSION(:),POINTER :: vnorm_i                 ! array to point to the corresponding normal vector

       !Limits rake direction
       REAL*4, DIMENSION(:,:), POINTER :: rake_lim

       !LAPLACIAN FILTER
       REAL*4, DIMENSION(:,:), POINTER :: vector_in, vector_out
       REAL*4, DIMENSION(:,:), POINTER :: lx, lz, lt, phi, theta
       REAL*4 h_lap, tol_lap
       INTEGER cont_i

       !LAPLACIAN TIME FILTER
       REAL*4, DIMENSION(:), POINTER :: tseries_in, tseries_out
       REAL*4, DIMENSION(:,:,:), POINTER :: filter, resfil

       !Depth preconditioner
       INTEGER depth_opt
       REAL depth_coef

       REAL*4, DIMENSION(:), POINTER :: rake_i, rake_p
       REAL*4 coef, coef2

       !Variables used to filter
       REAL*4 fc_obs, fc_syn, fc_grad
       !filter options 1=filter, 2=no filter, order of filter
       INTEGER f_obs, f_syn, f_grad, order
       !2D Gaussian filter:
       !r1: radius along dip, r2: along strike
       !repeat: times it is passed
       !filtgrad option: 1=filter, 2=no filter
       INTEGER r1, r2, repeat, filtgrad

       !Prior model information
       REAL*4, DIMENSION(:,:), POINTER :: p_model1d, p_model2d
       REAL balance1, balance2, balance3, balance4
       REAL quota1, quota2, quota3, quota4

       REAL progdt, progot
       INTEGER rak_case
       INTEGER niter_max, fwiopt

       !Triangular mesh
       REAL*4, DIMENSION(:), POINTER :: stk_ref, dip_ref
       REAL*4, DIMENSION(:), POINTER :: stk_coor, dip_coor
       REAL*4 dstk, ddip
       REAL*4, DIMENSION(:,:), POINTER :: values, xw, nodes_coor
       REAL*4, DIMENSION(:,:), POINTER :: tractionvec_ref
       INTEGER :: nelem, n_gauss, siz_gauss, opt_mesh, nnodes
       INTEGER, DIMENSION(:,:), POINTER :: corners

       !PIS options added
       REAL*4 :: percent_win
       REAL*4 :: percent_pri
       REAL*4 :: new_weight
       INTEGER :: wintochange
       INTEGER :: xhyp, zhyp         !for local laplacian filter


       END TYPE mesh

!---------------------
! Butterworth filter
!---------------------
TYPE butter

SEQUENCE

   ! Order of the butterworth filter
   INTEGER order

   ! Cutoff frequency
   REAL fc

   ! Filter coefficients
   REAL C(5,10)

   ! Coefficients defining the filter memory
   REAL D(2,10)

   ! Group delay in seconds
   REAL Tg

   ! Number of required 2nd order sections
   INTEGER NSections


END TYPE butter


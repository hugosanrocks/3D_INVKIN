       TYPE mesh

       SEQUENCE

       CHARACTER (LEN=4) dat, out                       	! Working directories dat=input, out=output

       !Input data coming from *.info files
       INTEGER*4 nsta, ncomp, msub      			! Number of stations, components and subfaults
       INTEGER*4 stress_opt					! Option of input data (0=GEODG3D, 1=DWN)
       INTEGER*8 stcomp, sta_i, comp_i, mjump            	! Total number of traces, counter on stations, components and staXcomp, jump file
       INTEGER*8 simsam, slipsam, interp_i, trac_i, interpadj_i ! Samples of: simulation, slip, forward interpolation, adjoint simulation
       INTEGER*4 delays                              		! Number of samples of delay (given by origin time of the force appli
								! to compute Green's functions)
       INTEGER tfft_i, nsdip, nsstk                                           ! Samples of traction array in frequency domain
       
       REAL*4 stk, dip, rak, dip_s, stk_s                                     ! Focal mechanism
       REAL*4 simt, simdt, slipdt, slipt 	              	! Simulation time, simulation dt, slip final time, slip dt, rise time
       REAL*4 moment, mu, rt, ot                                ! seismic moment, rise time, origin time

       !Data calculated and used for preprocessing
       !REAL*4 vnorm(3), vslip(3)                              	! Normal and slip unitary vectors
       REAL*4, DIMENSION(:,:), POINTER :: vnorm, vslip          ! Normal and slip unitary vectors
       REAL*4 stso(3,3), traction(3)                          	! Stress tensor and traction vector
         
       REAL*4, DIMENSION(:,:), POINTER :: stsi, stinterp	! Input stress tensor, input stress tensor interpolated
       REAL*4, DIMENSION(:,:), POINTER :: tractionv             ! Array to keep tractions (pseudo Green's functions)

       REAL*4, DIMENSION(:,:),POINTER :: fault
       REAL*4, DIMENSION(:),POINTER :: mus, tseries             ! Subfault's positions (x,y,z)


       CHARACTER (LEN=3) :: sta, sub                            ! Variables to compose file names
       CHARACTER (LEN=1) :: comp
       CHARACTER (LEN=20) :: FILE_S, FILE_OPT
       INTEGER dir_n, dir_i, subf_i                                     ! Number of directions, counter direction
       INTEGER, DIMENSION(:),POINTER :: vnorm_i                 ! array to point to the corresponding normal vector

       !Filter options
       INTEGER f_pre, order
       REAL*4 fc_pre


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


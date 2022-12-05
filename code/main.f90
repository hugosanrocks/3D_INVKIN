!----------------------------------------------------------------------
!     INV3DKIN
!     3D kinematic source inversion
!     by an adjoint-state method approach
!
!     ------------- Procedure info ---------------------------
!
!     This program:
!
!     1) Initializes all regularization and preconditioning
!        strategies
!
!     2) Forward modeling:
!
!      2.1) Using the sliprate model and unitary traction
!           vectors, computes the synthetic seismograms.
!
!    3) Adjoint:
!     
!     3.1) Estimate the residuals (res=syn-obs) at each station 
!          and compnent. Both, synthetics and observed seismograms
!          are reversed in time t'=T-t
!     3.2) The residuals at the recievers are projected to
!          the fault through the unitary traction used during the 
!          forward modeling and the total traction at the fault is 
!          estimated by convolving the untary traction vector
!          with the residuals (adjoint formulation).
!
!    4) fwi_option:
!
!     4.1) Searchs for a sliprate model that minimizes the 
!          misfit function (L2 Norm). Several methods are available
!          (PSTD, CG, LBFGS) to solve the optimization problem.
!----------------------------------------------------------------------
!    Hugo S. Sanchez Reyes 20/05/2018
!    Universite de Grenoble Alpes "UGA"
!    Supervisors: J. Tago, V.M. Cruz-Atienza, L. Metivier and J. Virieux
!---------------------------------------------------------------------
!#include "precis.h"

      program INV3DKIN

!----------------------Definition of variables--------!
      ! COMMON VARIABLES
      IMPLICIT NONE
      ! Define all the GLOBAL variables needed
      ! variables inside include/green.h
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh
      ! Optimization TOOLS_BOX
      include 'optim_type.h'
      type (optim_type) :: optim

      integer :: i

      !Local variables
!------------------------------------------------------!
      ! Timming variables
      real :: start, fini
!------------------------------------------------------!

      WRITE(6, *) '================================================='
      WRITE(6, *) '     3D Kinematic Source Inversion INV3DKIN      '
      WRITE(6, *) '================================================='


!=====================================================!
      ! ONLY IF YOU WANT TO PERFORM THE PREPROCESS
      ! STEPS UNCOMMENT THE NEXT LINE
!=====================================================!
      ! Prepare pseudo Green functions as needed
      !call preprocess(green_mesh)
!===================================================!

!===================================================!
!     READ INPUT INFO AND INITIALIZE ARRAYS
!===================================================!
      ! Read information, time samples, focal info
      ! reciever's geometry, optim. strategy, etc.
      call read_info(green_mesh)
      ! Initialize arrays
      call initialize(green_mesh)
!===================================================!


!===================================================!
!     Prepare regularization matrices
!===================================================!
      IF (green_mesh%msub .EQ. 1) THEN
       !Only one subfault = Point source approximation
      ELSE
       !Initialize regularizing terms
       call init_laplacian(green_mesh)
       !Penalize slip at the borders of the fault
       call edge(green_mesh)
      ENDIF
!===================================================!

!===================================================!
!     Read pseudo Green's functions                 !
!===================================================!
      ! Read tractions from file (time domain) 
      call read_time(green_mesh)
!===================================================!

!===================================================!
!     Assuming a REALTIME picking detect            !
!     time window for inversion. Now it only reads  !
!     the time windows from a file                  !
!===================================================!
      call windows(green_mesh)
!===================================================!

!===================================================!
!     INVERSION STRATEGIES                          !         !
!===================================================!

      !hess_opt = 1 Hessian   hess_opt = 2 inversion
      green_mesh%hess_opt = 3

      !time mark
      !call cpu_time(start)

      !Choose the inversion strategy SIS or PIS
      call invoption(green_mesh)

      !If PIS is carried out
      if ( green_mesh%for_opt .eq. 2 ) then
       stop
      endif

      !If SIS is going to be carried out
!===================================================!
!     FIRST FORWARD AND ADJOINT
!===================================================!
      call cpu_time(start)
      call forward(green_mesh)
      call cpu_time(fini)
      WRITE(6, *) '================================================='
      WRITE(6,*)  ' Time used for forward: ', fini-start,'seconds'
      WRITE(6, *) '================================================='

!print *, green_mesh%syn(1:5,1), 'syn'

      ! Estimate the tractions (adjoint problem) using as 
      ! forces the residuals
      call cpu_time(start)
      call adjoint(green_mesh)
      call cpu_time(fini)
      WRITE(6, *) '================================================='
      WRITE(6,*)  ' Time used for adjoint: ', fini-start,'seconds'
      WRITE(6, *) '================================================='
!print *, green_mesh%res(1:5,1), 'res'
!print *, green_mesh%tottrac(1:5,1), 'trac'

!===================================================!
      !initial misfit for the whole seismograms
      green_mesh%costdini = green_mesh%costd
!===================================================!

       !Arrange the graient in a 1D vector
       call read_grad(green_mesh)

!stop
       !Add penaty terms (regularization
       call penalty_terms(green_mesh)

!===================================================!
!     OPTIMIZATION STRATEGY 
!===================================================!
      ! Use any optimization strategy to follow the -gradient 
      ! direction to find a model (solution)
      call cpu_time(start)
      call fwi_option(green_mesh,optim)
      call cpu_time(fini)
      WRITE(6, *) '================================================='
      WRITE(6,*)  ' Time used for inversion: ', fini-start,'(sec)'
      WRITE(6, *) '================================================='
      !===========================================================!

      !write rake angles
      call rake_uni(green_mesh)

      !NOT WORKING ANYMORE
      ! Uncomment this line to check orthogonality of
      ! slip and normal vector2
      !call ortogonal(green_mesh)

      ! Destroy all the memory used
      call destroy_arrays(green_mesh)


      end program INV3DKIN


!----------------------------------------------------------------------
!
!     Inverse problem for 3D kinematic source inversion
!     by an adjoint-state method approach
!
!     ------------- Procedure info ---------------------------
!
!     This program:
!
!     1) Preprocess:
!
!      1.1) Reads the: 
!              - focal mechanism information.
!              - simulation (GEODG3D) information.
!              - slip rate model information (AXITRA).
!      1.2) Using this information, the program 
!           performs a coordinate transformation
!           slip rate (x,y,z).
!      1.3) Prepares the stress tensor by:
!             - Reading files P_CX TAU_*_CX SIGMA_**_CX.
!             - Rearranging this files into a unique file
!               for each station and component.
!             - Interpolating at the stress tensor at 
!               sliprate sampling rate.
!      1.4) Estimating the unitary traction vector at each 
!           staion and component.
!   
!     2) Forward:
!
!      2.1) Using the sliprate model and unitary traction
!           vector, computes the synthetic seismograms.
!
!    3) Adjoint:
!     
!     3.1) First, read the adjoint information (samples, etc.)
!     3.2) Estimate the residuals (r=syn-obs) at each station 
!          and compnent. Both, synthetics and observed seismograms
!          are reversed in time t'=T-t
!     3.4) The residuals at the recievers are retropropagated to
!          the fault and the total traction at the fault is 
!          estimated by convolving the untary traction vector
!          with the residuals (adjoint formulation).
!
!    4) fwi_pstd:
!
!     4.1) Searchs for a sliprate model that minimizes the 
!          cost function (L2). The steepest descent method
!          (PSTD) is used inside a loop, where the updating 
!          parameter to compute is the gradient (total
!          traction of adjoint problem).
!----------------------------------------------------------------------
!    Hugo S. Sanchez Reyes 29/05/15
!    Universite de Grenoble 1 "Joseph Fourier"
!    Supervisors: J. Tago, V.M. Cruz-Atienza, V. Hjorleifsdottir, J. Virieux
!---------------------------------------------------------------------
#include "precis.h"

      program INV3DKIN

!----------------------Definition of variables--------!
      ! COMMON VARIABLES
      IMPLICIT NONE
      ! Define all the GLOBAL variables needed
      ! variables inside include/green.h
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh
      ! Optimization TOOL_BOX
      include 'optim_type.h'
      type (optim_type) :: optim
      !INCLUDE 'varprec.h'
      !TYPE (var_type) :: vartype
      !MYFLOAT var1
      !Variables needed only by this program
      REAL pi/3.1415926535897932/
      integer i, j, k, iunit, fwiopt, m
!------------------------------------------------------!
      ! Timming variables
      real :: start, fini
      character*4 :: flag_d

      WRITE(6, *) '================================================='
      WRITE(6, *) '     3D Kinematic Source Inversion INV3DKIN      '
      WRITE(6, *) '================================================='


      ! Set working directories
      ! Directory containing input data files P_C* TUA_C* SIGMA_**_C*
      green_mesh%dat='dat/'
      ! Directory containing output data files
      green_mesh%out='out/'

      !CHANGE TO DOUBLE PRECISSION ==== NOT YET IMPLEMENTED
      !var_type%var1=1.d0
      !print *, var_type%var1

!=====================================================!
!=====================================================!
      ! ONLY IF YOU WANT TO PERFORM THE PREPROCESS
      ! STEPS UNCOMMENT THE NEXT LINE

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
!     INITIALIZE REGULARIZING MATRICES
!===================================================!
      IF (green_mesh%msub .EQ. 1) THEN
       !Only one subfault = Point source
      ELSE
      ! Initialize regularizing terms
       !call model_pri(green_mesh)
       !call laplacian(green_mesh)
      !Rupture time regularization
!      call exp_covar(green_mesh)
       !Correlation in time
       !call time_corr(green_mesh)
       !Penalize slip at the borders of the fault
       !call edge(green_mesh)
       !Penalize difference at hypocenter location
       !call hypo(green_mesh)
      ENDIF
!===================================================!

!===================================================!
!     Read pseudo Green's functions                 !
!===================================================!
      ! Read tractions from file (time domain) 
      call read_time(green_mesh)
!===================================================!

!===================================================!
!     Assuming a REALTIME picking (STALTA) detect   !
!     time window for inversion. Now it only reads  !
!     the time windows from a file                  !
!===================================================!
      call windows(green_mesh)
!===================================================!

!===================================================!
!     INVERSION STARTS HERE                         !
!===================================================!

      iunit = 33
      open(iunit,file=green_mesh%dat//'hessian_opt.dat',status='old',&
  &        action='read',form='formatted')
      read(iunit,*) green_mesh%hess_opt
      read(iunit,*) green_mesh%hess_cont
      close(iunit)

         if (green_mesh%for_opt .eq. 1) then
           write(*,*) ' Inversion option 1 (initial model 0) '
           call initwin(green_mesh)
         else
         endif

      !Beggining estimation
      call cpu_time(start)
     if ( green_mesh%hess_opt .eq. 1 ) then
       write(*,*) ' Hessian option 1, Hessian columns'
       k=green_mesh%modelsize1
     elseif ( green_mesh%hess_opt .eq. 2 ) then
       write(*,*) ' Hessian option 2, A^Tb estimation'
       k=1
     endif


     print *, ' Total number of parameters: ', k
      do i=1,k
       green_mesh%hess_cont = i
       print *, i

      !Rake is fixed
      if (green_mesh%rake_opt .eq. 1) then
         write(*,*) ' Rake is not an unknown rake_opt=1'
         if (green_mesh%for_opt .eq. 1) then
           call coor_trans(green_mesh)
         else
         endif
       endif

       print *, green_mesh%model1(green_mesh%hess_cont), green_mesh%model1(i+1), i, i+1

!===================================================!
!     FIRST FORWARD AND ADJOINT
!===================================================!
!      call cpu_time(start)
      call forward(green_mesh)
!      call cpu_time(fini)
!      WRITE(6, *) '================================================='
!      WRITE(6,*)  ' Time used for forward: ', fini-start,'seconds'
!      WRITE(6, *) '================================================='

      ! Estimate the tractions (adjoint problem) using as 
      ! forces the residuals
!      call cpu_time(start)
      call adjoint(green_mesh)
!      call cpu_time(fini)
!      WRITE(6, *) '================================================='
!      WRITE(6,*)  ' Time used for adjoint: ', fini-start,'seconds'
!      WRITE(6, *) '================================================='
!===================================================!

       call read_grad(green_mesh)

       call write_hessian(green_mesh)

       !green_mesh%rak=green_mesh%rak-pi
      enddo
 
      call cpu_time(fini)
      WRITE(6, *) '================================================='
      WRITE(6,*)  ' Time used for for estimation: ', fini-start,'seconds'
      WRITE(6, *) '================================================='


      call destroy_arrays(green_mesh)


      end program INV3DKIN




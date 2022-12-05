!--------------------------------------------------------!
!     3D Kinematic Seismic Source Inversion              !
!          by an adjoint-state method                    !
!--------------------------------------------------------!
!------- MAIN FOR SENSITIVITY ANALYSIS ------------------!
!
!     This program:
!
!
!----------------------------------------------------------------------
!    Author:
!    Hugo S. Sanchez Reyes 06/09/16
!    Universite de Grenoble 1 "Joseph Fourier"
!    Supervisors: J. Tago, V.M. Cruz-Atienza, V. Hjorleifsdottir, 
!                 L. Metivier, J. Virieux
!---------------------------------------------------------------------

      program sensitivity

!----------------------Definition of variables--------------------------------------

      implicit none
      !COMMON VARIABLES
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      !local variables
      real :: sample


      ! Set working directories
      ! Directory containing input data files P_C* TUA_C* SIGMA_**_C*
      green_mesh%dat='dat/'
      ! Directory containing output data files
      green_mesh%out='out/'

      !CHANGE TO DOUBLE PRECISSION ==== NOT YET IMPLEMENTED
      !var_type%var1=1.d0
      !print *, var_type%var1


!===================================================!
!     READ INPUT INFO AND INITIALIZE ARRAYS
!===================================================!
      ! Read information, time samples, focal info
      ! reciever's geometry, optim. strategy, etc.
      call read_info(green_mesh)
      ! Initialize arrays
      call initialize(green_mesh)
!===================================================!

      call read_time(green_mesh)

      call windows(green_mesh)

      ! 3 = start from a given initial model and invert the whole history !
      !===================================================================!

      !Rake is fixed
      if (green_mesh%rake_opt .eq. 1) then
         write(*,*) ' Rake is not an unknown rake_opt=1'
         if (green_mesh%for_opt .eq. 1) then
           write(*,*) ' Inversion option 1 (initial model 0) '
           call initwin(green_mesh)
           call coor_trans(green_mesh)
         endif
       endif

      green_mesh%syn(:,:) = 0.
      !Read observed seismograms at stations
      call read_obs(green_mesh)
      call residual(green_mesh)
      sample = green_mesh%res(40,1)
      green_mesh%res(:,:) = 0.
      green_mesh%res(40,1) = sample
      call write_residual(green_mesh)

      !Flush the array to store total traction
      green_mesh%tottrac(:,:)=0.

      call conadjtime(green_mesh)

      call read_grad(green_mesh)

      call write_grad(green_mesh)

      !call destroywin(green_mesh)
      call destroy_arrays(green_mesh)
      endprogram sensitivity

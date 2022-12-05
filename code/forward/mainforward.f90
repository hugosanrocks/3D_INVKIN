!--------------------------------------------------------
!     3D Kinematic Seismic Source Inversion
!     by an adjoint-state method
!--------------------------------------------------------
!------------- MAIN FORWARD MODELING PROGRAM ------------
!
!     This program:
!
!     2) Forward:
!
!      2.1) Using the sliprate module for each subfault -> green_mesh%slipr
!           and the pseudo-Green's functions            -> green_mesh%traction
!           computes the synthetic seismograms.
!
!----------------------------------------------------------------------
!    Author:
!    Hugo S. Sanchez Reyes 29/06/15
!    Universite de Grenoble 1 "Joseph Fourier"
!    Supervisors: J. Tago, V.M. Cruz-Atienza, V. Hjorleifsdottir, 
!                 L. Metivier, J. Virieux
!---------------------------------------------------------------------

      program main_forward

!----------------------Definition of variables--------------------------------------
      !COMMON VARIABLES
      implicit none
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

!     Variables needed only by this program
      integer iunit
      real start, fini
!------------------------------------------------------

!     Set working directories
!     Directory containing input data files P_C* TUA_C* SIGMA_**_C*
      green_mesh%dat='dat/'
!     Directory containing output data files
      green_mesh%out='out/'

      !Read general information and prepare stress tensor files
      call read_info(green_mesh)

      call initialize(green_mesh)

      call initwin(green_mesh)


      if (green_mesh%rake_opt .eq. 1) then
        write(*,*) ' Rake is not an unknown rake_opt=1'
        if (green_mesh%for_opt .eq. 1) then
          write(*,*) ' Forward option 1 (dat/vitesse.out) '
          call coor_trans(green_mesh)
        elseif (green_mesh%for_opt .eq. 2) then
          write(*,*) ' Forward option 2 (dat/modelpri1d.dat) '
        endif
      elseif(green_mesh%rake_opt .eq. 2) then
        write(*,*) ' Rake is an unknown rake_opt=1'
         if (green_mesh%for_opt .eq. 1) then
           call coor_trans(green_mesh)
           print *, ' Forward option 1 (dat/vitesse.out) '
         elseif (green_mesh%for_opt .eq. 99) then
           !<F3>call read_model(green_mesh)
           print *, ' Forward option 2 (dat/model.out) '
         elseif (green_mesh%for_opt .eq. 2) then
           call read_modelpri(green_mesh)
           print *, ' Forward option 3 ( dat/modelpri.dat)'
         endif
       endif
!==================================================!


      !Read Green's functions
      call read_time(green_mesh)

      !call green_interp(green_mesh)

         !Assuming a REALTIME picking (STALTA) detect time window for inversion
         !======Now it only reads from a file==================================
         call windows(green_mesh)


      !Compute the forward problem and compute first synthetics
      !asociated with the first model (0 ZEROS)
      !lensyn = dimension of synthetic vectors
       call cpu_time(start)
        call forward(green_mesh)
       call cpu_time(fini)
       write(*,*) '****************************************************'
       print *, 'Forward modeling done in:', fini-start, 'seconds'
       write(*,*) '****************************************************'
      call write_syn(green_mesh)

call rake_uni(green_mesh)
      call write_model(green_mesh)
      !call rake_angle(green_mesh)

      call destroywin(green_mesh)
      !call destroy_arrays(green_mesh)


!      deallocate(green_mesh%fault)
!      deallocate(green_mesh%slipmod)
!      deallocate(green_mesh%model,green_mesh%model2,green_mesh%grad2)
!      deallocate(green_mesh%tractionvec,green_mesh%syn,green_mesh%slipr)
!      deallocate(green_mesh%slip)
!      deallocate(green_mesh%gradad,green_mesh%slipr2)
!      deallocate(green_mesh%obs,green_mesh%cd,green_mesh%cm)
!      deallocate(green_mesh%ce,green_mesh%ct,green_mesh%rtimes)
!      deallocate(green_mesh%rsamp,green_mesh%diag2)
!      deallocate(green_mesh%diag)
!      deallocate(green_mesh%la)
!      deallocate(green_mesh%modelp)
!      deallocate(green_mesh%model2p)
!      deallocate(green_mesh%tseries)
!      deallocate(green_mesh%samwin)
!      deallocate(green_mesh%synsam)
!      deallocate(green_mesh%win)


      !Frequency domain
!      deallocate(green_mesh%tracf,green_mesh%slipf)


      end program main_forward

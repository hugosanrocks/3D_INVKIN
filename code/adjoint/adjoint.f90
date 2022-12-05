!--------------------------------------------------------
!     3D Kinematic Seismic Source Inversion
!     by an adjoint-state method
!--------------------------------------------------------
!------------- ADJOINT SUBROUTINE -----------------------
!
!     This program:
!
!    3) Adjoint:
!
!     3.2) Estimate the residuals (r=syn-obs) at each station 
!          and compnent. Both, synthetics and observed seismograms
!          are reversed in time t'=T-t
!     3.4) The residuals at the recievers are retropropagated to
!          the fault and the total traction at the fault is 
!          estimated by convolving the pseudo-Green's functions
!          with the residuals (adjoint formulation).
!
!----------------------------------------------------------------------
!    Author:
!    Hugo S. Sanchez Reyes 29/06/15
!    Universite de Grenoble 1 "Joseph Fourier"
!    Supervisors: J. Tago, V.M. Cruz-Atienza, V. Hjorleifsdottir, 
!                 L. Metivier, J. Virieux
!---------------------------------------------------------------------


      subroutine adjoint(green_mesh)

!----------------------Definition of variables--------------------------------------

      implicit none

      !COMMON VARIABLES
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      !local variables
      integer i
      real start, fini

         !Read observed seismograms at stations
         call read_obs(green_mesh)

         if (green_mesh%hess_opt .eq. 1) then
            print *, ' Hessian option on, observations set to zero'
            green_mesh%obs(:,:) = 0.
         elseif (green_mesh%hess_opt .eq. 2) then
            print *, ' Hessian option on, synthetics set to zero'
            green_mesh%syn(:,:) = 0.
         else
            print *, ' Hessian option off, observations read'
         endif


         !Compute residuals for all stations and components
        call residual(green_mesh)

         !Flush the array to store total traction
         green_mesh%tottrac(:,:)=0.

         !Perform convolution in the time domain
         call conadjtime(green_mesh)

!==================================================!
!       ATTENTION: UNCOMMENT THE NEXT SECTION
!  IF CONVOLUTIONS ARE PERFORM IN FREQUENCY DOMAIN
!==================================================!
!         green_mesh%tottrac(:,:)=0.

!         call cpu_time(start)
!         do i = 1,green_mesh%msub
!           !Jump inside traction files for m-subfaults
!           green_mesh%tfft_i = i
!           call adj_trac(green_mesh)
!         enddo
!         call cpu_time(fini)
!         print *, 'Time for adjoint:', fini-start
!=================================================!



         
      end subroutine adjoint

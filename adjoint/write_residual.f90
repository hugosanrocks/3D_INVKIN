!--------------------------------------------------------
!     3D Kinematic Seismic Source Inversion
!     by an adjoint-state method
!--------------------------------------------------------
!------------- WRITE VELOCITY RESIDUALS SUBROUTINE ------
!
!     This program:
!
!     Writes into ASCII files the corresponding velocity
!     residuals at each station and component
!
!     THIS SUBROUTINE IS INDEPENDENT OF THE INVERSION
!----------------------------------------------------------------------
!    Author:
!    Hugo S. Sanchez Reyes 29/06/15
!    Universite de Grenoble 1 "Joseph Fourier"
!    Supervisors: J. Tago, V.M. Cruz-Atienza, V. Hjorleifsdottir, 
!                 L. Metivier, J. Virieux
!---------------------------------------------------------------------


       subroutine write_residual(green_mesh)

       !COMMON VARIABLES
       IMPLICIT NONE
       INCLUDE 'green.h'
       TYPE (mesh) :: green_mesh

       !Local variables
       INTEGER i, ii, jj, k, iunit, lenr

       !length of one residual element
       INQUIRE(iolength=lenr) green_mesh%res(1,1)

       !Unit to open files
       iunit=12

       k=1
       do ii=1,green_mesh%nsta
        do jj=1,green_mesh%ncomp
         write(green_mesh%sta,'(I3.3)') ii
         write(green_mesh%comp,'(I1.1)') jj
         !Write residuals into a binary file
         OPEN(iunit,FILE=green_mesh%out//'res_S'//green_mesh%sta//'_C'//green_mesh%comp//'.bin',&
    &    status='unknown',&
    &    FORM='UNFORMATTED',&
    &    ACCESS='DIRECT',recl=lenr)
         !ASCII file with residuals (used to check)
         OPEN(99,FILE=green_mesh%out//'res_S'//green_mesh%sta//'_C'//green_mesh%comp//'.ascii',&
    &    status='unknown')
         do i=1,green_mesh%samwin(ii,green_mesh%dowin)     !syn_i = interp_i
          write(iunit,rec=i) green_mesh%res(i,k)
          write(99,*) green_mesh%res(i,k)
         enddo
         close(iunit)
         close(99)
        k=k+1
        enddo
       enddo

       end subroutine write_residual

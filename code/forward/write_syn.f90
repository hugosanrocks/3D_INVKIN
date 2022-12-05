!--------------------------------------------------------
!     3D Kinematic Seismic Source Inversion
!     by an adjoint-state method
!--------------------------------------------------------
!------------- WRITE SYNTHETIC SEISMOGRAMS SUBROUTINE ---
!
!     This program:
!
!     Writes into ASCII files the corresponding synthetic
!     seismograms to each station and every component
!
!     THIS SUBROUTINE IS INDEPENDENT OF THE INVERSION
!----------------------------------------------------------------------
!    Author:
!    Hugo S. Sanchez Reyes 29/06/15
!    Universite de Grenoble 1 "Joseph Fourier"
!    Supervisors: J. Tago, V.M. Cruz-Atienza, V. Hjorleifsdottir, 
!                 L. Metivier, J. Virieux
!---------------------------------------------------------------------

       subroutine write_syn(green_mesh)

       !COMMON VARIABLES
       IMPLICIT NONE
       INCLUDE 'green.h'
       TYPE (mesh) :: green_mesh

       !local variables
       INTEGER ii, jj, k, iunit, iunit2, iunit3, iunit4
       INTEGER start, j
       REAL    t, dt     
       REAL    maxamp, maxamp2

      call read_obs(green_mesh)
      

      !Unit to open files
      iunit=22
      iunit2=23
      iunit3=24
      iunit4=25

      dt = green_mesh%intdt
      start= 1

      k=1
      do ii=1,green_mesh%nsta
       write(green_mesh%sta,'(I3.3)') ii
       open(iunit3,file=green_mesh%out//'syn_S'//green_mesh%sta//'.win',&
    &  status='unknown')
       write(iunit3,*) '# NUM_WIN = ', green_mesh%wininv-1
       !do jj=1,green_mesh%wininv-1
       write(iunit3,*) 1, real(green_mesh%samwin(k,1))*dt, green_mesh%twin(ii,green_mesh%synwin), '0'
       !enddo
       close(iunit3)
       do jj=1,green_mesh%ncomp
         t = 0.
!        write(green_mesh%sta,'(I3.3)') ii
        write(green_mesh%comp,'(I1.1)') jj
        open(iunit2,file=green_mesh%out//'syn_S'//green_mesh%sta//'_C'//green_mesh%comp//'.head',&
    &   status='unknown')
        write(iunit2,*) '# NPTS = ', 875
        open(iunit3,file=green_mesh%out//'obs_S'//green_mesh%sta//'_C'//green_mesh%comp//'.head',&
    &   status='unknown')
        write(iunit3,*) '# NPTS = ', 875  !changed 
       !Write synthetic seismogram ASCII file (used to check)
        open(iunit,file=green_mesh%out//'syn_S'//green_mesh%sta//'_C'//green_mesh%comp//'.ascii',&
    &   status='unknown')
        open(iunit4,file=green_mesh%out//'obs_S'//green_mesh%sta//'_C'//green_mesh%comp//'.a',&
    &   status='unknown')
        do j=start,green_mesh%samwin(ii,green_mesh%synwin)
         write(iunit,*) t, green_mesh%syn(j,k)
         write(iunit4,*) t, green_mesh%obs(j,k)
         t = t + dt
        enddo
       close(iunit) !Close writing units
       close(iunit4)
       !Maximum values in the plot
       maxamp  = maxval(abs(green_mesh%syn(1:green_mesh%samwin(ii,green_mesh%wininv),k)))
       maxamp2 = maxval(abs(green_mesh%obs(:,k)))
       maxamp = max(maxamp,maxamp2)
       write(iunit3,*) '# PLOT_MAX = ', maxamp
       write(iunit3,*) '# T_START = ',  0.d0
       write(iunit3,*) '# T_END = ',    real(green_mesh%lenobs-1)*green_mesh%intdt
       close(iunit3)
       write(iunit2,*) '# PLOT_MAX = ', maxamp
       write(iunit2,*) '# T_START = ',  0.d0
       write(iunit2,*) '# T_END = ',    real(green_mesh%lenobs-1)*green_mesh%intdt
       close(iunit2)

       k=k+1
       enddo
      enddo
      


      end subroutine write_syn






       subroutine writesyn_std(green_mesh)

       !COMMON VARIABLES
       IMPLICIT NONE
       INCLUDE 'green.h'
       TYPE (mesh) :: green_mesh

       !local variables
       INTEGER ii, k, iunit
       INTEGER start, j


       !call read_obs(green_mesh)

       !Unit to open files
       iunit =22


       k=1
       do ii=1,green_mesh%nsta
        write(green_mesh%sta,'(I3.3)') ii
        open(iunit,file=green_mesh%out//'syn_S'//green_mesh%sta//'.ascii',&
    &    status='unknown')
         write(iunit,*) '# NPTS = ', green_mesh%trac_i  !changed 
         write(iunit,*) 'sanchez_2016 ', 875  !changed 
         write(iunit,*) '# NPTS = ', 875  !changed 
         !Write synthetic seismogram ASCII file (used to check)
         do j=start,green_mesh%samwin(ii,green_mesh%synwin)
          write(iunit,*) green_mesh%syn(j,k:k+2)
         enddo
         k=k+3
       enddo
      


      end subroutine writesyn_std

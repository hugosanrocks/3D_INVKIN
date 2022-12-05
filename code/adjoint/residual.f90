!--------------------------------------------------------
!     3D Kinematic Seismic Source Inversion
!     by an adjoint-state method
!--------------------------------------------------------
!------------- RESIDUALS SUBROUTINE -----------------------
!
!     This program:
!
!    3) Adjoint:
!
!     3.2) Estimate the residuals (r=syn-obs) at each station 
!          and compnent. Both, synthetics and observed seismograms
!          are reversed in time t'=T-t
!
!----------------------------------------------------------------------
!    Author:
!    Hugo S. Sanchez Reyes 29/06/15
!    Universite de Grenoble 1 "Joseph Fourier"
!    Supervisors: J. Tago, V.M. Cruz-Atienza, V. Hjorleifsdottir, 
!                 L. Metivier, J. Virieux
!---------------------------------------------------------------------

      subroutine residual(green_mesh)

      !COMMON VARIABLES
      implicit none
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      integer i, j, k

      !Flush residual array
      green_mesh%res(:,:)=0.d0
      
      !Flush the cost function before computation
      green_mesh%costa=0.d0

      !ALL residuals have same length      
!      green_mesh%res(:,:) = green_mesh%syn(:,:) - green_mesh%obs(:,:)

      
      !Invert for different time windoes
      k = 1
      do i=1,green_mesh%nsta
        do j=1,green_mesh%ncomp
          green_mesh%res(1:green_mesh%samwin(i,green_mesh%synwin),k) = &
  &       green_mesh%syn(1:green_mesh%samwin(i,green_mesh%synwin),k) - green_mesh%obs(1:green_mesh%samwin(i,green_mesh%synwin),k)
          k = k + 1
        enddo
      enddo
print *, green_mesh%obs(1:5,1), 'obs'


      !Multiply by data covariance matrix (Weighting factors)
      call prod_res(green_mesh)
print *, green_mesh%costa
      !Print only if you want to check
      !call write_residual(green_mesh)

      !Transform residuals to Frequency domain
      !turn off this if convolutions are done in tme
      !call resi_fft(green_mesh)

      end subroutine residual



      subroutine filter_grad_time(green_mesh)

      USE INIT_BUTTERWORTH_MOD
      USE FILTFILT_BUTTERWORTH_MOD

      !COMMON VARIABLES
      implicit none
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh
      TYPE (butter) :: butt

      integer ii
      integer nrec, nstep
      real dt
      real,dimension(:),pointer :: tseries

          !Butterworth definition
          butt%order = green_mesh%order
          butt%fc = green_mesh%fc_grad
          nstep = green_mesh%interp_i
          nrec  = green_mesh%ncomp*green_mesh%msub
          dt    = green_mesh%slipdt
          allocate(tseries(nstep))
          !Initialize the butterworth filter
          CALL INIT_BUTTERWORTH(butt, dt)

          !-------------------------------------------------------------------
          !Parameters needed
       do ii=1,nrec
          !Asign the time series to be filtered
          tseries(:) = green_mesh%tottrac(:,ii)
          !PRINT ONLY TO CHECK THE FILTER
          !do jj=1,nstep
          ! write(54,*) green_mesh%tottrac(jj,ii)
          !enddo
          ! Apply the butterworth filter
          CALL FILTFILT_BUTTERWORTH(tseries, butt, nstep)
          !copy_tseries(:,2) = tseries(:)
          green_mesh%tottrac(:,ii) = tseries(:)
          !PRINT TO CHECK THE FILTER
          !do jj=1,nstep
          ! write(55,*) green_mesh%tottrac(jj,ii)
          !enddo
          ! Print the resultsprint input and filtered signals to output file
          !-----------------------------------------------------------------
          ! END OF OPTIONAL FILTER APPLIED TO GREEN'S FUNCTIONS
          !-----------------------------------------------------------------
      enddo

      deallocate(tseries)
      endsubroutine filter_grad_time



      subroutine filter_syn(green_mesh)

      USE INIT_BUTTERWORTH_MOD
      USE FILTFILT_BUTTERWORTH_MOD

      !COMMON VARIABLES
      implicit none
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh
      TYPE (butter) :: butt

      integer ii
      integer nrec, nstep
      real dt
      real,dimension(:),pointer :: tseries

          !Butterworth definition
          butt%order = green_mesh%order
          butt%fc = green_mesh%fc_syn
          nstep = green_mesh%lenobs
          nrec  = green_mesh%stcomp
          dt    = green_mesh%slipdt
          allocate(tseries(nstep))
          !Initialize the butterworth filter
          CALL INIT_BUTTERWORTH(butt, dt)

          !-------------------------------------------------------------------
          !Parameters needed
       do ii=1,nrec
          !Asign the time series to be filtered
          tseries(:) = green_mesh%syn(:,ii)
          !PRINT ONLY TO CHECK THE FILTER
          !do jj=1,nstep
          ! write(54,*) green_mesh%tottrac(jj,ii)
          !enddo
          ! Apply the butterworth filter
          CALL FILTFILT_BUTTERWORTH(tseries, butt, nstep)
          !copy_tseries(:,2) = tseries(:)
          green_mesh%syn(:,ii) = tseries(:)
          !PRINT TO CHECK THE FILTER
          !do jj=1,nstep
          ! write(55,*) green_mesh%tottrac(jj,ii)
          !enddo
          ! Print the resultsprint input and filtered signals to output file
          !-----------------------------------------------------------------
          ! END OF OPTIONAL FILTER APPLIED TO GREEN'S FUNCTIONS
          !-----------------------------------------------------------------
      enddo

      deallocate(tseries)
      endsubroutine filter_syn




      subroutine read_obs (green_mesh)

  USE INIT_BUTTERWORTH_MOD
  USE FILTFILT_BUTTERWORTH_MOD

      !COMMON VARIABLES
      implicit none
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh
      TYPE (butter) :: butt

      integer ii, k, iunit, p
      integer nrec, nstep
      real dt, vector_in(3)

      p=1
      iunit = 51

       do ii=1,green_mesh%nsta
         write(green_mesh%sta,'(I3.3)') ii
         !Unit to read observations
         OPEN(iunit,FILE=green_mesh%dat//'obs_S'//green_mesh%sta//'.dat',&
    &    status='old',action='read')
         do k=1,green_mesh%lenobs
          read(iunit,*) vector_in(1:3)
          green_mesh%obs(k,p:p+2) = vector_in(1:3)
         enddo
         close(iunit)
        p=p+3
       enddo  !loop over number of stations


       !Filter observations       
       if (green_mesh%f_obs .eq. 1) then

          !-------------------------------------------------------------------
          !OPTIONAL FILTER APPLIED TO GREEN'S FUNCTIONS
          butt%order = green_mesh%order
          butt%fc = green_mesh%fc_obs
          nstep = green_mesh%lenobs
          nrec  = green_mesh%stcomp
          dt    = green_mesh%slipdt
          !Initialize the butterworth filter
          CALL INIT_BUTTERWORTH(butt, dt)

          !-------------------------------------------------------------------
          !Parameters needed
       do ii=1,green_mesh%stcomp
          !Asign the time series to be filtered
          green_mesh%tseries(:) = green_mesh%obs(:,ii)
          ! Apply the butterworth filter
          CALL FILTFILT_BUTTERWORTH(green_mesh%tseries, butt, nstep)
          !copy_tseries(:,2) = tseries(:)
          green_mesh%obs(:,ii) = green_mesh%tseries(:)
          ! Print the resultsprint input and filtered signals to output file
          !-----------------------------------------------------------------
          ! END OF OPTIONAL FILTER APPLIED TO GREEN'S FUNCTIONS
          !-----------------------------------------------------------------
      enddo
      else
       print *, ' Observations not filtered'
      endif

       call weight_cov(green_mesh)

      end subroutine read_obs



      subroutine resi_fft (green_mesh)

      implicit none
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      include "fftw3.f"

      integer j, k
      integer*8 plan_forward
      !double precision in(green_mesh%lensyn)
      !double complex out(green_mesh%lensynf)
      complex out(green_mesh%lensynf)
      real in(green_mesh%lensyn)

      call sfftw_plan_dft_r2c_1d_ ( plan_forward, green_mesh%lensyn, in, out,&
     &  FFTW_ESTIMATE )

        do k = 1, green_mesh%stcomp
        in(:)=0.
          do j = 1, green_mesh%interp_i
            in(j) = green_mesh%res(j,k)
          enddo

        !do j = 1, green_mesh%lensyn
        !  write(16,*) j, in(j)
        !enddo
        !execute the plan to transform the IN data to
        !the OUT FFT coefficients.

        call sfftw_execute_ ( plan_forward )

        !do j = 1, green_mesh%lensynf
        !  write(16,*) j, out(j)
        !enddo


     !USED TO CHECK BACK FFT VALUES
     ! call sfftw_plan_dft_c2r_1d_ ( plan_backward, green_mesh%lensyn, out, in,&
     !&  FFTW_ESTIMATE )
     ! call sfftw_execute_ ( plan_backward )
     ! do j=1,green_mesh%lensyn
     ! write(16,*) j,in(j)/dble(green_mesh%lensyn)
     ! enddo
     ! call sfftw_destroy_plan_ ( plan_backward )
     ! stop

      green_mesh%resif(:,k) = out(:)

      enddo

      call sfftw_destroy_plan_ ( plan_forward )

      return
      end subroutine resi_fft






      subroutine weight_cov(green_mesh)

      implicit none
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

 !     real,dimension(:),allocatable :: x, y
      integer i,j,k
      integer iunit
      real    weights(3)


      green_mesh%cd(:,:)=0.

      !Read weights given for each recording
      k=1
      iunit=22
      open(iunit,file=green_mesh%dat//'weights.dat',status='old',action='read')
      do i=1,green_mesh%nsta
       read(iunit,*) weights
       do j=1,3
        green_mesh%cd(k,k) = weights(1+(j-1))
        k=k+1
       enddo
      enddo
      close(iunit)


      endsubroutine weight_cov



      subroutine prod_res(green_mesh)

      !COMMON VARIABLES
      implicit none
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      integer i, j, k, l, iunit, mm
      real weights(3), cost, costa, al, be
      real, dimension(:), allocatable :: resvec, vec
      real, dimension(:,:), allocatable :: weightm, resvect

      !MKL matrix multiplication coefficients
      al = 1.d0
      be = 0.d0
      costa = 0.

      !All recordings with weight from file (Cd matrix)
      if (green_mesh%weig .eq. 2) then
       iunit = 22
       open(iunit,file=green_mesh%dat//'weights.dat',status='old',action='read')
       k = 1
       do i=1,green_mesh%nsta
        mm = green_mesh%samwin(i,green_mesh%wininv)
        read(iunit,*) weights
        allocate(weightm(mm,mm),resvect(1,mm))
        allocate(resvec(mm),vec(mm))
        resvec(:) = 0.
        vec(:) = 0.
        weightm(:,:) = 0.
        do j=1,green_mesh%ncomp
            !Multiplication of matrices 
            resvec(:) = green_mesh%res(1:mm,k)
            do l=1,mm
               weightm(l,l) = weights(j)
            enddo
            call sgemm('N','N',mm,1,mm,al,        &
      &       weightm,mm,resvec,mm,be,vec,mm)
            resvect(1,:) = vec(:)
            call sgemm('N','N',1,1,mm,al,        &
      &       resvect,1,vec,mm,be,cost,mm)
            !Cumulative cost
            costa = costa + cost
            !Residuals weighted
            green_mesh%res(1:mm,k) = vec(:)
          k = k + 1
        enddo
       deallocate(weightm,resvec,resvect,vec)
       enddo
       close(iunit)
      !Option for weights != 2
      else
       !call fcost(green_mesh%res,green_mesh%stcomp,green_mesh%syn_i,costa)
       write(*,*) ' ERROR: fcost function not called'
       write(*,*) ' Value of weighting option incorrect'
       write(*,*) ' Check this value at dat/fwioption.info'
       stop
      endif

      green_mesh%costa = costa * 0.5
      green_mesh%costd = costa * 0.5


      endsubroutine prod_res




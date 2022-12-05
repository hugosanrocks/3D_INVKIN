!--------------------------------------------------------
!     3D Kinematic Seismic Source Inversion
!     by an adjoint-state method
!--------------------------------------------------------
!------------- SYNTHETIC VELOCITY RECORDINGS SUBROUTINE -
!
!     This program:
!
!      1) Uses the slip rate model (3 components) from 
!         green_mesh%slipr
!      2) Reads the pseudo-Green's functions 
!         (unitary traction) TRACT_SXXX.bin.
!      3) Convolves these two vectors for each station
!         according to the forward formulation.
!      4) Stack the corresponding components of these
!         convolutions to estimate the synthetic
!         records at each station SXXX and component CX
!      5) Writes the resulting synthetics
!
!----------------------------------------------------------------------
!    Author:
!    Hugo S. Sanchez Reyes 29/06/15
!    Universite de Grenoble 1 "Joseph Fourier"
!    Supervisors: J. Tago, V.M. Cruz-Atienza, V. Hjorleifsdottir, 
!                 L. Metivier, J. Virieux
!---------------------------------------------------------------------


      !==============================================================!
      ! Subroutine contime(green_mesh)
      ! Estimation of forward modeling with P0 approach.
      ! This subroutine computes the time convolution 
      ! between slip-rate functions and the unitary
      ! traction vector at every cartesian node. Then, 
      ! every contribution is summed over the whole mesh.
      !==============================================================!


      subroutine contime(green_mesh)

       !COMMON VARIABLES
       IMPLICIT NONE
       INCLUDE 'green.h'
       TYPE (mesh) :: green_mesh
       ! FFT Library

      integer i, j, k, ii, jj, kk, mm, cont, m
      real scalfac
      real, dimension(:), allocatable :: x, y, z

      !Finite fault scaling factor
      scalfac = green_mesh%stk_s*1000.*green_mesh%dip_s*1000. 

      !Point source scaling factor
      !green_mesh%moment/(green_mesh%msub)

      allocate(x(green_mesh%ntint),y(green_mesh%trac_i),z(green_mesh%lensyn))

      cont = 0
       do i=1,green_mesh%msub
        do j=1,green_mesh%nsta
        do m=1,green_mesh%ncomp
         do k=1,green_mesh%ncomp   !Loop over 3 traction components
          ii = m + (i-1)*3
          cont = k + (m-1)*3 + (j-1)*9 + &
  &              (i-1)*green_mesh%stcomp*green_mesh%ncomp
          jj = m + (j-1)*3
          kk = m + (k-1)*3 + (j-1)*green_mesh%ncomp*green_mesh%ncomp+ &
  &            (i-1)*green_mesh%stcomp*3
          mm = k + (j-1)*3
!          used to check elements to be convolved
          !print *, 'slip', ii, 'conv tract', kk, 'syn', mm, 'contime'
          x(:) = green_mesh%slip(:,ii)  !ii = ii
          y(:) = green_mesh%tractionvec(:,kk) !cont = kk   jj= mm
!         !Convolve slip rate model with the Green functions
          call conv(x,green_mesh%ntint,y,green_mesh%trac_i,&
  &                z,green_mesh%lensyn)
!         !Discrete convolution term
          z(:)=z(:)*green_mesh%intdt*scalfac
!         !Stack the three convolutions X, Y, Z
         green_mesh%syn(:,mm)=green_mesh%syn(:,mm)+z(green_mesh%delays+1:green_mesh%delays+green_mesh%samwin(j,green_mesh%synwin))
         !green_mesh%syn(:,mm)=green_mesh%syn(:,mm)+z(1:green_mesh%samwin(j,green_mesh%synwin))
         enddo     !Loop over 3 traction components
        enddo
       enddo
      enddo


      !Filter the forward problem
      if ( green_mesh%f_syn .eq. 1 ) then
        call filter_syn(green_mesh)
      endif

      deallocate(x,y,z)
      endsubroutine contime



      !==============================================================!
      ! Subroutine contime_p1(green_mesh)
      ! Estimation of forward modeling with P1 approach.
      ! This subroutine computes the time convolution 
      ! between slip-rate functions and the unitary
      ! traction vector at every cartesian node. Then, 
      ! the 2D space integral is computed through the 
      ! the general 2D trapezoidal rule over the cartesian
      ! grid.
      !==============================================================!


      subroutine contime_p1(green_mesh)

       !COMMON VARIABLES
       IMPLICIT NONE
       INCLUDE 'green.h'
       TYPE (mesh) :: green_mesh

       !=======================================!
       ! Forward model using Trapezoidal rule  !
       !=======================================!

       integer i, j, k, ii, kk, cont, m
       integer fi, fj, fk, count2, count1
       real scalfac
       real, dimension(:), allocatable :: x, y, z, acum
       real, dimension(:,:), allocatable :: subsyn, faultsyn

       !Variables for 2D integration
       real vx(green_mesh%nsstk), vy(green_mesh%nsdip)
       real dx, dy, result
       integer nx, ny, re

       !Size of subfaults along strike and dip
       dx = green_mesh%stk_s*1000.
       dy = green_mesh%dip_s*1000.

       !Number of subfaults along strike and dip
       nx = green_mesh%nsstk
       ny = green_mesh%nsdip

       !Build coordinate vectors for 2D integration
       vx(1) = 0.
       vy(1) = 0.
       do i=2,nx
         vx(i) = vx(i-1)+dx
       enddo
       do i=2,ny
         vy(i) = vy(i-1)+dy
       enddo

       !Scale factor according to area and moment    
       scalfac = 1.

!       re = 4*13608*green_mesh%trac_i
!       open(33,file='TRAC_4.bin',status='unknown',&
!  &        form='unformatted',ACCESS='DIRECT',recl=re)
!       write(33,rec=1) green_mesh%tractionvec(:,1:13608)
!       close(33)


       !scalfac=green_mesh%moment/(green_mesh%msub)
       allocate(x(green_mesh%ntint))
       allocate(y(green_mesh%trac_i))
       allocate(z(green_mesh%lensyn))
       allocate(acum(green_mesh%lensyn))
       allocate(subsyn(green_mesh%lensyn,green_mesh%msub))
       allocate(faultsyn(ny,nx))


      cont = 1
      do i=1,green_mesh%nsta
       do j=1,green_mesh%ncomp
        do k=1,green_mesh%msub
         acum(:) = 0.
         do m=1,green_mesh%ncomp
          ii = m + (k-1)*green_mesh%ncomp
          kk = m + (i-1)*green_mesh%ncomp**2 & 
  &          + (k-1)*green_mesh%stcomp*green_mesh%ncomp &
  &          + (j-1)*green_mesh%ncomp
!         used to check elements to be convolved
          !print *, 'slip', ii, 'conv t', kk, 'f',k,'syn',cont
          x(:) = green_mesh%slip(:,ii)  !ii = ii
          y(:) = green_mesh%tractionvec(:,kk) !cont = kk   jj= mm
!         !Convolve slip rate model with the Green functions
          call conv(x,green_mesh%ntint,y,green_mesh%trac_i,&
  &                z,green_mesh%lensyn)
          acum(:) = acum(:) + z(:)*green_mesh%intdt
!         !Discrete convolution term *dt
         enddo
         subsyn(:,k) = acum(:)
        enddo
       !Integral along fault 2D trapezoidal
       count1 = 1
       do fk = green_mesh%delays+1,green_mesh%delays+green_mesh%samwin(i,green_mesh%synwin)
        count2 = 1
        do fi = 1,ny
         do fj = 1,nx
          faultsyn(fi,fj) = subsyn(fk,count2)
          count2 = count2 + 1
         enddo
        enddo
        call trapz2D(vx,vy,nx,ny,faultsyn,result)
        green_mesh%syn(count1,cont) = result
        count1 = count1 + 1
       enddo
       !Finish 2D integration        
       cont = cont + 1
       enddo
      enddo

      !Filter the forward problem
      if ( green_mesh%f_syn .eq. 1 ) then
        call filter_syn(green_mesh)
      endif

       deallocate(x,y,z,acum)
       deallocate(subsyn,faultsyn)
      endsubroutine contime_p1



      !==============================================================!
      ! Subroutine contime_p1_trian(green_mesh)
      ! Estimation of forward modeling with P1 approach
      ! for non regular grids (triangular elements).
      ! This subroutine computes the time convolution 
      ! between slip-rate functions and the unitary
      ! traction vector at every node location. Then, 
      ! the 2D space integral is computed through the 
      ! the sum of every triangular element contribution.
      ! To estimate each contribution, we use triangles.f90 module.
      !==============================================================!

      subroutine contime_p1_trian(green_mesh)

      !load module related to triangular elements
      use triangles

       !COMMON VARIABLES
       IMPLICIT NONE
       INCLUDE 'green.h'
       TYPE (mesh) :: green_mesh

       integer i, j, k, ii, kk, cont, m
       integer fk, count1
       real, dimension(:), allocatable :: x, y, z, acum

       !Variables for 2D integration
       real :: result

       !scalfac=green_mesh%moment/(green_mesh%msub)
       allocate(x(green_mesh%ntint))
       allocate(y(green_mesh%trac_i))
       allocate(z(green_mesh%lensyn))
       allocate(acum(green_mesh%lensyn))

      cont = 1
      do i=1,green_mesh%nsta
       do j=1,green_mesh%ncomp
        do k=1,green_mesh%msub
         green_mesh%values(:,k) = 0.
         do m=1,green_mesh%ncomp
          ii = m + (k-1)*green_mesh%ncomp
          kk = m + (i-1)*green_mesh%ncomp**2 &
  &          + (k-1)*green_mesh%stcomp*green_mesh%ncomp &
  &          + (j-1)*green_mesh%ncomp
!         used to check elements to be convolved
          !print *, 'slip', ii, 'conv t', kk, 'f',k,'syn',cont
          x(:) = green_mesh%slip(:,ii)  !ii = ii
          y(:) = green_mesh%tractionvec(:,kk) !cont = kk   jj= mm
!         !Convolve slip rate model with the Green functions
          call conv(x,green_mesh%ntint,y,green_mesh%trac_i,&
  &                z,green_mesh%lensyn)
          green_mesh%values(:,k) = green_mesh%values(:,k) + &
  &                                z(:)*green_mesh%intdt
!         !Discrete convolution term *dt
         enddo
        enddo
       !Integral across the 2D fault surface
       count1 = 1
       do fk = green_mesh%delays+1,green_mesh%delays+green_mesh%samwin(i,green_mesh%synwin)
        call integrate_surface(green_mesh,fk,result)
        green_mesh%syn(count1,cont) = result
        count1 = count1 + 1
       enddo
       !Finish 2D integration        
       cont = cont + 1
       enddo
      enddo

      !Filter the forward problem
      if ( green_mesh%f_syn .eq. 1 ) then
        call filter_syn(green_mesh)
      endif

      deallocate(x,y,z,acum)
      endsubroutine contime_p1_trian


      !==============================================================!
      ! Subroutine contime_parallel(green_mesh)
      ! Estimation of forward modeling with P0 approach
      ! in parallel using OpenMP.
      ! This subroutine computes the time convolution 
      ! between slip-rate functions and the unitary
      ! traction vector in parallel at every cartesian 
      ! node. Then, every contribution is summed over 
      ! the whole mesh.
      !==============================================================!

      subroutine contime_parallel(green_mesh)

      use omp_lib
      !COMMON VARIABLES
      IMPLICIT NONE
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      integer i, td_id
      real,dimension(:,:),allocatable :: slipr, array2, array3
      integer nsstk, nsdip, nsub, nt, ntg, lensyn, cols, nsta, ncomp, stcomp
      integer nslipr
      real dt
      integer(kind=8) s1, s2, rate

      call system_clock(count_rate=rate)

      !some initial parameters
      nsta = green_mesh%nsta
      ncomp = green_mesh%ncomp
      stcomp = green_mesh%nsta*green_mesh%ncomp
      nsstk = green_mesh%nsstk
      nsdip = green_mesh%nsdip
      nsub = green_mesh%msub
      nt = green_mesh%ntint
      ntg = green_mesh%trac_i
      nslipr = green_mesh%msub*green_mesh%ncomp
      lensyn = green_mesh%trac_i+green_mesh%ntint-1
      dt = green_mesh%intdt

      cols=green_mesh%msub*green_mesh%stcomp*green_mesh%ncomp

      allocate(slipr(nt,nslipr),array2(ntg,cols),array3(lensyn,stcomp))

      slipr(:,:) = green_mesh%slip(:,:)
      array2(:,:) = green_mesh%tractionvec(:,:)
      array3(:,:) = 0.

      !parallel subroutine to compute convolutions
      call system_clock(s1)
      !$OMP PARALLEL SHARED(array2,slipr,array3) PRIVATE(td_id,i)
       td_id = OMP_get_thread_num()
       print *, 'Thread starting: ', td_id
       call convolutions(nsta,ncomp,nsub,slipr,array2,array3,nt,nslipr,ntg,cols,lensyn,stcomp,td_id)
      !$OMP END PARALLEL
      call system_clock(s2)
      write(*,*) 'Time parallel: ', real(s2-s1)/real(rate)

      !scaling factor
      green_mesh%syn(1:green_mesh%lenobs,:) = array3(11:green_mesh%lenobs+11,:)*0.1*1000.*1000.

      !Filter the forward problem
      if ( green_mesh%f_syn .eq. 1 ) then
        call filter_syn(green_mesh)
      endif


      deallocate(slipr,array2,array3)
      endsubroutine contime_parallel



      !==============================================================!
      ! Subroutine convolutions
      ! This subroutine computes the convolutions required
      ! by the forward problem in parallel (contime_parallel)
      ! using OpenMP subtoutines.
      !==============================================================!
      
      subroutine convolutions(nsta,ncomp,nsub,slipr,array2,array3,nt,nslipr,ntg,cols,lensyn,stcomp,td_id)

      use omp_lib
      implicit none
      integer i, ii, jj, kk, j, k, m
      integer, intent(inout) :: nsta, ncomp, nsub
      integer, intent(inout) :: td_id, nt, nslipr, ntg, cols, lensyn, stcomp
      real, intent(inout) :: slipr(nt,nslipr), array2(ntg,cols), array3(lensyn,stcomp)
      real x(nt), y(ntg), z(lensyn)

       !$OMP DO REDUCTION(+:array3) SCHEDULE(STATIC) PRIVATE(x,y,z,jj)
       do i = 1, nsta
        do j = 1, ncomp
         do k = 1, nsub
          do m = 1, ncomp
           ii = m + (k-1)*ncomp
           kk = m + (i-1)*ncomp**2 + (k-1)*ncomp*ncomp*nsta + (j-1)*ncomp
           y(:) = array2(:,kk)
           x(:) = slipr(:,ii)
           call conv(x,nt,y,ntg,z,lensyn)
           jj = 1+(i-1)*ncomp+(j-1)
           array3(:,jj) = array3(:,jj) + z(:)
          enddo
         enddo
        enddo
       enddo
       !$OMP END DO

      end subroutine convolutions



    !==============================================================!
    ! Subroutine trapz2D
    ! Computes the integral of a given function across a 
    ! cartesian regular grid using the general 2D trapezoidal
    ! rule.
    !==============================================================!

    subroutine trapz2D(x,y,nx,ny,mat,integral)

    integer,intent(in) :: nx, ny
    integer a, b, c, d, ix, iy, i, j
    real,intent(in) :: x(nx), y(ny), mat(ny,nx)
    real,intent(out) :: integral
    real hx, hy
    real const1, sum1, sum2, sum3, sum4, sum5, sum6

    ix = 1
    iy = 1
    hx = x(ix+1) - x(ix)
    hy = y(ix+1) - y(iy)
    a = ix
    b = nx
    c = iy
    d = ny

    const1 = 0.25*hx*hy

    !Edges
    sum1 = const1*( mat(c,a)+mat(c,b)+mat(d,a)+mat(d,b)  )

    sum2 = 0.
    sum3 = 0.
    sum4 = 0.
    sum5 = 0.
    sum6 = 0.

    do i=iy+1,d-1
      sum4 = sum4 + mat(i,a)
      sum5 = sum5 + mat(i,b)
      do j=ix+1,b-1
        sum6 = sum6 + mat(i,j)
      enddo
     enddo
    sum6 = const1*(4.*sum6)

    do i=ix+1,b-1
      sum2 = sum2 + mat(c,i)
      sum3 = sum3 + mat(d,i)
    enddo

    sum2 = const1*(2.*sum2 + 2.*sum3 +2.*sum4 +2.*sum5)

    integral = sum1 + sum2 + sum6

    endsubroutine trapz2D



      !==============================================================!
      ! Subroutine integrate_surface
      ! 1. Get gaussian points to use in the integration 
      ! 2. Transform gaussian points into the element coordinates
      ! 3. Interpolate values at required coordinates inside element
      ! 4. Integrate with interpolated values and gaussian weights
      ! 5. Add the contribution of each element to total surface
      !==============================================================!

      subroutine integrate_surface(green_mesh,cont,fval)

      use triangles
      implicit none
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh

      !time sample where the surface integral is estimated
      integer, intent(in) :: cont
      !output value of integral
      real, intent(inout) :: fval
      !counter
      integer :: i
      !interpolation points
      real :: pint(green_mesh%siz_gauss,2)

      !Set how many points to use for integration
      !check triangle_gauss_points subroutine in triangles.f90

      ! Coordinates of nodes along strike and dip
      green_mesh%stk_coor(:) = green_mesh%nodes_coor(1,:)
      green_mesh%dip_coor(:) = green_mesh%nodes_coor(2,:)


      fval = 0.
      do i = 1,green_mesh%nelem
       !Get coordinates of gauss points in triangle eleme
       pint = triangle_transform(green_mesh%stk_coor(green_mesh%corners(:,i)),&
  &           green_mesh%dip_coor(green_mesh%corners(:,i)),green_mesh%xw,&
  &           green_mesh%siz_gauss)

       !Integrate the surface contain by triangle element
       fval = fval + integrate_triangle(pint,&
  &           green_mesh%stk_coor(green_mesh%corners(:,i)),&
  &           green_mesh%dip_coor(green_mesh%corners(:,i)),&
  &           green_mesh%xw,green_mesh%values(cont,&
  &           green_mesh%corners(:,i)),green_mesh%siz_gauss)
      enddo
      
      endsubroutine integrate_surface




      !==============================================================!
      ! Subroutine read_fft
      ! Reads the Green functions saved in the frequency domain
      ! used for the forward problem .
      !==============================================================!

      subroutine read_fft(green_mesh)

       !COMMON VARIABLES
       IMPLICIT NONE
       INCLUDE 'green.h'
       TYPE (mesh) :: green_mesh

       integer bit, iunit
       integer*8 reclent

       ! Bytes used for each element of traction vector (frequency)
       INQUIRE(iolength=bit) green_mesh%tracf(1,1)
       ! Total bytes used for traction matrix (frequency)
       reclent=bit*green_mesh%ncomp*green_mesh%lensynf*green_mesh%stcomp*green_mesh%msub

       iunit=16
       OPEN(iunit,file='dat/TRACT_fft.bin',status='unknown',&
 &          form='unformatted',ACCESS='DIRECT',recl=reclent)
         read(iunit,rec=1) green_mesh%tracf(:,:)
         !read(iunit,rec=green_mesh%tfft_i) green_mesh%tracf(:,:)
         !do i=1,green_mesh%lensynf
         !   write(15,*) green_mesh%tracf(i,1)
         !enddo
       close(iunit)

      end subroutine read_fft



      !==============================================================!
      ! Subroutine confft
      ! Computes the forward problem in the frequency domain.
      !==============================================================!

      subroutine confft(green_mesh)

       !COMMON VARIABLES
       IMPLICIT NONE
       INCLUDE 'green.h'
       TYPE (mesh) :: green_mesh
       ! FFT Library
       include "fftw3.f"

!      double complex in1(green_mesh%lensynf), in2(green_mesh%lensynf), inh(green_mesh%lensynf)
!      double precision outr(green_mesh%lensyn)
      complex in1(green_mesh%lensynf), in2(green_mesh%lensynf), inh(green_mesh%lensynf)
      real outr(green_mesh%lensyn)

      integer i, j, k, ii, jj, cont

!  Set up a plan, and execute the plan to transform the IN data to
!  the OUT FFT coefficients.

       ii = green_mesh%ncomp*green_mesh%stcomp
       jj = 1 
       do i=1,green_mesh%nsta             !stations
         do j=1,green_mesh%ncomp          !slip components
          do k=1,green_mesh%ncomp         !traction components
!           slip component in1, traction component in2
            in1(:) = green_mesh%slipf(:,j)
            in2(:) = green_mesh%tracf(:,&
  &         j+(k-1)*green_mesh%ncomp+((i-1)*(green_mesh%ncomp**2))+(green_mesh%tfft_i-1)*ii) 
            inh(:) = in1(:) * in2(:)      !convolution
            cont = k+(i-1)*3
!            print *, j,  j+(k-1)*green_mesh%ncomp+((i-1)*(green_mesh%ncomp**2))+(green_mesh%tfft_i-1)*ii, cont, 'confreq'
            call backfft(green_mesh,inh,green_mesh%lensynf,&
  &                      green_mesh%lensyn,outr,cont)
          enddo
          jj = jj + 1
         enddo
       enddo

      return
      end subroutine confft


      !==============================================================!
      ! Subroutine backfft
      ! Inverse Fourier Transform using the FFTW3 package
      !==============================================================!


      subroutine backfft(green_mesh,outr,nc,n,inh,cont)

      implicit none
      ! COMMON VARIABLES
      INCLUDE 'green.h'
      TYPE (mesh) :: green_mesh
      !Dimension of vectors
      integer, intent(inout) :: nc, n, cont
      !Vectors to be convolved
!      double complex, intent(inout) :: outr(nc)
!      double precision, intent(inout) :: inh(n)
      complex, intent(inout) :: outr(nc)
      real, intent(inout) :: inh(n)

      ! FFT Library
      include "fftw3.f"

      integer*8 plan_backward
      real scalfac

      call sfftw_plan_dft_c2r_1d_ ( plan_backward, n, outr, inh,&
     &  FFTW_ESTIMATE )

      call sfftw_execute_ ( plan_backward )

      !For homogeneous media, scalfac has the same mu for all subfaults
      !evereything done in preprocess stage
      !scalfac=green_mesh%moment/(green_mesh%mu*green_mesh%msub)
      scalfac=green_mesh%moment/(green_mesh%msub)
      
      !print *, green_mesh%delays+1, green_mesh%delays+green_mesh%interp_i,'check'
      green_mesh%syn(:,cont) = green_mesh%syn(:,cont) + &
  &                 (inh(green_mesh%delays+1:green_mesh%delays+green_mesh%interp_i)*scalfac*green_mesh%slipdt)/dble(n)
      
      call sfftw_destroy_plan_ ( plan_backward )

      return
      endsubroutine backfft

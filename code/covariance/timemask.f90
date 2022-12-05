

       subroutine time_mask(green_mesh)

       !Subroutine in charge of reading the 
       !weighting for the prior model

       !In the past it was used to designed
       !and write it. Now it only reads it
       !from a file designed previously by 
       !a matlab subroutine

       !The weights for the prior model have
       !to be inside dat/weight_prior.dat
       !expected rupture times are also read from 
       !dat/time_rupture.dat

       !COMMON VARIABLES
       IMPLICIT NONE
       INCLUDE 'green.h'
       TYPE (mesh) :: green_mesh

       integer i, j, k, iunit, iunit2


       iunit = 33
       iunit2= 34

       open(iunit,file=green_mesh%dat//'weight_prior.dat',status='old',&
  &         action='read')
       open(iunit2,file=green_mesh%dat//'time_rupture.dat',status='old',&
  &         action='read')
       k=1
       do j=1,green_mesh%interp_i
       do i=1,green_mesh%msub
        read(iunit,*) green_mesh%diag(k)
        read(iunit2,*) green_mesh%rup(k)
        k=k+1
       enddo
       enddo
       close(iunit)
       close(iunit2)

       open(iunit,file=green_mesh%dat//'timemask.dat',status='unknown',&
  &         action='write')
       k=1
       do j=1,green_mesh%interp_i
       do i=1,green_mesh%msub
        write(iunit,*) green_mesh%diag(k)
        k=k+1
       enddo
       enddo      
       close(iunit)


      endsubroutine time_mask


      !======================================!
      !      GAUSSIAN DISTRIBUTION 1D        !
      !======================================!
      !
      !!This subroutine computes a
      !!normalized 1D gauss distribution
      !
      !!Imput:
      !!nt:  length of the series
      !!t:   serie to compute gaussian PDF
      !!mu:  mean of the PDF
      !!sig: sigma of the PDF
      !
      !!Output:
      !!gauss: PDF
      !======================================!
  
       subroutine gauss_1d_dist(t,gauss,nt,mu,sig)

       !This subroutine computes a
       !normalized 1D gauss distribution

       !Imput:
       !nt:  length of the series
       !t:   serie to compute gaussian PDF
       !mu:  mean of the PDF
       !sig: sigma of the PDF

       !Output:
       !gauss: PDF

       implicit none
       integer i
       integer, intent(inout) :: nt
       real, intent(in) :: t(nt), mu, sig
       real, intent(inout) :: gauss(nt)
       real dt, minv
       real pi/3.1415926535897932/

       dt = t(2) - t(1)

       do i=1, nt
         gauss(i) = 1./ (sig*sqrt(pi))* & 
  &               exp(-1.*((real(i)-1.)*dt-mu)**2./sig**2.)
       enddo

       minv = maxval(gauss)
       
       gauss(:) = gauss(:)/minv
       gauss(:) = gauss(:)*0.98
       do i=1, nt
        gauss(i) = 1. - gauss(i)
       enddo
       

       endsubroutine gauss_1d_dist



       subroutine emg_1d_dist(t,gauss,nt,h,tr,sg,ta)

       !This subroutine computes an
       !approximated exponential modified gaussian 

       !Check Kevin Lan, James W. Jorgenson (2001)

       !Imput:
       !nt:  length of the series
       !t:   serie to compute gaussian PDF
       !sg:  sigma of gaussian part
       !ta:  relax distant of exponential part
       !h:   peak magnitude
       !tr:  ordinate of peak

       !Output:
       !gauss: PDF

       implicit none
       integer i, j, samp
       integer, intent(inout) :: nt
       real, intent(in) :: t(nt), h, tr, sg, ta
       real, intent(inout) :: gauss(nt)
       real pi/3.1415926535897932/
       real c1, tini, tout(nt), tin(nt), trmap, dt

       !Controling parameters
       !A=0.08898
       !B=0.09656
       !Ordinate in time of pead
       !tr=9.345
       !Peak magnitude
       !h=0.3155
       !Sigma of Gaussian part
       !sg=0.57602
       !Relax distance exponential
       !ta=0.81019

       dt = t(2) - t(1)
       samp = floor(tr/dt)
       tini = real(samp-1)*dt
       tin(1) = tini
       do i=2, nt
        tin(i) = tin(i-1) - dt
       enddo

       !tin = t
       !trmap = tr
       trmap = 0.

       gauss(:) = 0.

       do i=1, nt
         c1 = (2.*sg**2.)+ta*(tin(i)-trmap)
         if (c1 .gt. 0.) then
            gauss(i) = 100. - h*exp(-1.*(tin(i)-trmap)**2./c1)
         else
            gauss(i) = 100. - 0.
         endif
       enddo



       endsubroutine emg_1d_dist



       !OLD CODE USED TO DESIGN THE WEIGHTS

!!       real, dimension(:,:), allocatable :: f
!!       real, dimension(:), allocatable :: t, gauss
!!       real dt, lambda
!!       real h, tr, sg, ta, mu, sig
!!       integer i, ii, j, iunit, iunit2, tsam, sub, ind(green_mesh%msub), record, reclen, k, nuc(9)


       !!dt = green_mesh%slipdt
       !!tsam = green_mesh%interp_i
       !!sub = green_mesh%msub
       !!lambda = 0.1 !0.001 !0.1

       !!allocate(t(tsam),f(tsam,sub),gauss(tsam))

       !!f(:,:) = 1.

       !!ii=1
       !!do i=1,green_mesh%nsdip
       !! do j=1,green_mesh%nsstk
       !!  do k=1,green_mesh%interp_i
       !   f(k,ii) = 1.*((i-1)*2.)            !Set maximum penalty value as 1
       !!  enddo
       !!  ii = ii+1
       !! enddo
       !!enddo
       !Time axis used for estimating the penalty 
       !at each time and subfault
       !!t(1) = 0.
       !!do i=2,tsam
       !! t(i) = t(i-1) + dt
       !!enddo

       !Gaussian distribution
       !call gauss_1d_dist(t,gauss,nt,mu,sig)

       !Ordinate in time
       !tr=9.345
       !Peak magnitude
       !!h=99.98!0.3155
       !Sigma of Gaussian part
       !!sg=3.
       !Relax distance exponential
       !!ta=7.



       !Exponential modified Gaussian
       !call emg_1d_dist(t,gauss,tsam,h,tr,sg,ta)
       !print to check gaussian subroutines
       !do i=2,tsam
       ! print *, t(i), gauss(i)
       !enddo

       !Uncomment this if you want constant rupture front
       !Constant travel time   VS_MAX
       !!ind(:) = green_mesh%rsamp(:)
       !!do i=1,green_mesh%msub
       !!  if (ind(i) .le. 0) then
       !!    ind(i) =1
       !!  endif
       !!  !write(68,*) ind(i)
       !!enddo
       !Estimated from eikonal solver
       !!open(iunit,file=green_mesh%dat//'ttimes.info',status='old',&
!!  &  !!     action='read')
       !read(iunit,*) ind(:)
       !!close(iunit)

      !If flag_d = DETC it will detect the maxima of each slip-rate
      !history and build the regularization with respect to the
      !time postion of those maximum values

      !!if (green_mesh%flag_d .ne. 'DETC') then
      !! print *, ' Using travel times from Eikonal solver'
      !!  do i=1,green_mesh%msub
!before plotting         if ((ind(i) .le. 45) .and. (i .lt.  green_mesh%nsstk*(green_mesh%nsdip-(green_mesh%nsdip/6))    )) then
      !!   if ((ind(i) .le. 35) .and. (i .lt.  green_mesh%nsstk*(green_mesh%nsdip-(green_mesh%nsdip/6))    )) then
!!          ind(i) = 1
      !!   endif
          !print *, ind(i), i, green_mesh%flag_d
      !!  enddo
      !! else
      !! print *, ' Using travel times from detected maximum peaks'
      !!  do i=1,green_mesh%msub
!!          ind(i) = green_mesh%rsamp(i)
          !print *, ind(i), i
      !!  enddo
      !! endif

      !! j = 1
      !! do ii=1,green_mesh%nsdip
      !!  do k=1,green_mesh%nsstk
      !!   do i=ind(j),tsam
      !!    f(i,j) = 1.!exp(-1. * (t(i) - t(ind(j)  )) / lambda )
          !donout
      !!    if ( t(i) .gt. t(ind(j))+ 6.) then
      !!        f(i,j) = 1.
      !!    endif
          !if ((k .le. 9) .and. (t(i) .gt. 12.))  then
          !   f(i,j)  = 1.
          !endif
      !!   enddo

          !Adding the edge condition
      !!    if (( ii .gt. 15 )  ) then
      !!     f(:,j) = 1.
      !!    elseif (( k .gt. 33)) then
      !!     f(:,j) = 1.
      !!    elseif (( k .lt. 4)) then
      !!     f(:,j) = 1.
      !!    elseif  ( ii .lt. 3)  then
      !!     f(:,j) = 1.
          !elseif ( ( k .eq. 29) ) then
          ! f(:,j) = 3.
         ! elseif (( ii .eq. 17) ) then
         !  f(:,j) = 1.5
         ! elseif (( ii .eq. 18) ) then
         !  f(:,j) = 2.
      !!    endif


          !Adding the edge condition
          !if (( ii .gt. 9 )  ) then
          ! f(:,j) = 2.
          !elseif (( k .eq. 2)) then
          ! f(:,j) = 2.
          !elseif (( k .eq. 1)) then
          ! f(:,j) = 3.
          !elseif ( ( k .eq. 28) ) then
          ! f(:,j) = 2.
          !elseif ( ( k .eq. 29) ) then
          ! f(:,j) = 3.
         ! elseif (( ii .eq. 17) ) then
         !  f(:,j) = 1.5
         ! elseif (( ii .eq. 18) ) then
         !  f(:,j) = 2.
         ! endif

      !!    j = j + 1
      !!  enddo
      !! enddo

     
       !f(:,:) = 0.
      !! j = 1
      !! do ii=1,green_mesh%nsdip
      !!  do k=1,green_mesh%nsstk
      !!   tr = real(ind(j))*green_mesh%slipdt
      !!   mu = real(ind(j))*green_mesh%slipdt
      !!   sig = 5.
!         call emg_1d_dist(t,gauss,tsam,h,tr,sg,ta)
!         call gauss_1d_dist(t,gauss,tsam,mu,sig)
!         f(:,j) = gauss(:)
      !!   j = j+1
      !!  enddo
      !! enddo

      !! deallocate(f,t,gauss)

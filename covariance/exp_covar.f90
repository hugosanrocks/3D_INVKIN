        subroutine exp_covar(green_mesh)

        !COMMON VARIABLES
        IMPLICIT NONE
        INCLUDE 'green.h'
        TYPE (mesh) :: green_mesh

        integer i, j, k, n(2), iunit
        real lambda, lambda0, sigmam, factor, hypo(1,3), t, dt
        integer hyp
        real,dimension(:,:),allocatable :: dist
        real,dimension(:),allocatable :: disthyp
        integer,dimension(:),allocatable :: front

        !Flush array to save node coordinates
        green_mesh%fault(:,:)=0.

        !File where to read the subfault positions
        iunit=10
        open(iunit,file=green_mesh%dat//'fault.pos',status='old',&
        action='read',access='DIRECT',recl=green_mesh%msub*4*green_mesh%ncomp)
        !Read subfault positions (x,y,z)
        read(iunit,rec=1) green_mesh%fault(:,:)
        close(iunit)
        !Used to check the subfault positions
        !do j=1,green_mesh%msub
        !  write(77,*) green_mesh%fault(j,:)
        !enddo

        !Estimate distances between nodes
        do i=1,green_mesh%msub
         do j=1,green_mesh%msub
          green_mesh%distcorr(i,j) = sqrt( &
  &                (green_mesh%fault(j,1)-green_mesh%fault(i,1))**2 +& 
  &                (green_mesh%fault(j,2)-green_mesh%fault(i,2))**2 +&
  &                (green_mesh%fault(j,3)-green_mesh%fault(i,3))**2  )*1000.
        enddo
       enddo

       !Read weighting model matrix
       call time_mask(green_mesh)


       endsubroutine exp_covar


       !PREVIOUS CODE FOR REGULARIZATION ESTIMATION

        !allocate(disthyp(green_mesh%msub))
        !allocate(front(green_mesh%msub))

        !THIS MUST BE CHANGED TO DETECT HYPOCENTER AUTMATICALY
        !hypo(1,:) = green_mesh%hypo(:)

        !SAME AS Mathilde Radiguet THESIS EQ (3.3) 
        !Controling parameters of covariance matrix
        !lambda = 3000.!500.     !9 km
        !lambda0 = 1500.     !3 km
        !sigmam = 0.1    !0.5 meters
        !factor = ( sigmam * (lambda / lambda0 ) )**2
        !siv1
        !hyp = 210
        !jap
        !hyp = 181


       !Distance from each subfault center to the hypocenter
       !do i=1,green_mesh%msub
       !  disthyp(i) = sqrt( (hypo(1,1)-green_mesh%fault(i,1))**2 +&
!  &                (hypo(1,2)-green_mesh%fault(i,2))**2 +&
!  &                (hypo(1,3)-green_mesh%fault(i,3))**2  )*1000.
!         !write(66,*) disthyp(i)
!       enddo

!      n(1) = 1
!      n(2) = 0
!      green_mesh%win(:,:) = 0     !flush array
!      do k = 1,2
!       t  = real(n(k))          !initial time
!       dt = 1.          !duration of rupture windows
!       do i=1,6
!        t = t + dt
!        do j=1,green_mesh%msub
!           if (disthyp(j) .lt. t * green_mesh%vrup*1. ) then !4620.*0.7 ) then
!             green_mesh%win(j,k) = green_mesh%win(j,k) + 1
!           endif
!        enddo
!       enddo
!      enddo

      !KUMAMOTO
      !n(1) = 3
      !n(2) = 2
      !green_mesh%win(:,:) = 0     !flush array
      !do k = 1,2
      ! t  = real(n(k))          !initial time
      ! dt = 2.          !duration of rupture windows
      ! do i=1,6
      !  t = t + dt
      !  do j=1,green_mesh%msub
      !     if (disthyp(j) .lt. t * green_mesh%vrup*1. ) then !4620.*0.7 ) then
      !       green_mesh%win(j,k) = green_mesh%win(j,k) + 1
      !     endif
      !  enddo
      ! enddo
      !enddo


!       !IMPORTANT NUM WIN + 1 - GREEN_MESH%WIN
!       green_mesh%win(:,:) = 7 - green_mesh%win(:,:)
!       do i=1,green_mesh%msub
!         if ( green_mesh%win(i,1) .eq. 2) then
!            green_mesh%win(i,1) = 1
!         endif
!       enddo

!       open(44,file='dat/subfaults_win.out',status='unknown')
!       do i=1,green_mesh%msub
!         write(44,*) green_mesh%win(i,1:2)
!       enddo
!       close(44)

       !Travel times
       !TO BE CHANGED FOR VARIABLE VS VELOCITY
!       do i=1,green_mesh%msub
!       green_mesh%rtimes(i) = disthyp(i) / green_mesh%vrup !4620. !4620 MAX velo siv1
!        if ( disthyp(i) .gt. 8000.) then
!        green_mesh%rtimes(i) = disthyp(i) / (green_mesh%vrup*0.7) !4620. !4620 MAX velo siv1
!        endif
!       enddo
!       do i=1,green_mesh%msub
!       green_mesh%rsamp(i) = floor( (disthyp(i) / green_mesh%vrup ) / green_mesh%slipdt )  !-2 !-8 is bestdist(hyp,i) = disthypo(i)
       !write(31,*) disthyp(i), green_mesh%rsamp(i)
!       enddo


       !deallocate(disthyp,front)


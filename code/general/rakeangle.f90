      subroutine rake_angle(green_mesh)

      implicit none
      INCLUDE 'green.h'
!     optim type strcuture
      TYPE (mesh) :: green_mesh

      ! Arrays to save rake angle and rake angle prior
      real, dimension(:), allocatable :: angle, angle_p, ang_vx, ang_vy, mean, des
      integer, dimension(:), allocatable :: ind
      real :: stk, dip, modu
      integer i, j, k, cont, n, ini, fin
      integer :: iunit, ii, cont2
      real pi2, cost_ang, trigg
      real pi/3.1415926535897932/
      pi2 = pi * 2.

        trigg = 0.1    !Trigger to measure angle and cost
        print *, ' Measuring the rake angle'
 
        !samples, subfaults
        !Arrays to store the rake angles
        allocate(angle(green_mesh%interp_i*green_mesh%msub))
        allocate(angle_p(green_mesh%interp_i*green_mesh%msub))
        allocate(ang_vx(green_mesh%interp_i*green_mesh%msub))
        allocate(ang_vy(green_mesh%interp_i*green_mesh%msub))
        allocate(ind(green_mesh%interp_i*green_mesh%msub))

        n = green_mesh%interp_i
        angle(:) = 0.
        angle_p(:) = 0.

        ini = 1
        fin = green_mesh%interp_i

        iunit=32
        open(iunit,file='angle.out',status='unknown')


        if (green_mesh%rake_opt .eq. 1) then
        call model_d(green_mesh%model,green_mesh%model2, &
  &                  green_mesh%interp_i,green_mesh%msub,&
  &                  green_mesh%slipm,green_mesh%dir_i,&
  &                  green_mesh%vnorm_i)
        endif

        green_mesh%rake_p(:) = 0.

        cost_ang = 0.
        cont2 = 0
        !Arrange slip rate model in 1D vector
        j = 1
        do i=1,green_mesh%msub
          !Select one subfault at each time
          cont=1+(i-1)*green_mesh%interp_i*2
            do k=1,green_mesh%interp_i
              stk = green_mesh%model2(cont)
              dip = green_mesh%model2(cont+green_mesh%interp_i)
              modu = sqrt(stk**2. + dip**2.)
              if (modu .gt. trigg) then
               ind(j) = 1
              else
               ind(j) = 0
              endif
               !Take the angle
               !write(777,*) stk, dip, slip1sf(k)     
               if (dip .gt. 0.) then
                 if (stk .gt. 0.) then
                  angle(j) = (pi - atan(dip/stk))
                 elseif (stk .lt. 0.) then
                  angle(j) =  (atan(-1.*dip/stk))
                 endif
               elseif (dip .lt. 0.) then
                 if (stk .gt. 0.) then
                   angle(j) = (pi + atan(-1.*dip/stk))
                 elseif (stk .lt. 0.) then
                   angle(j) = (pi*2. - atan(dip/stk))
                 endif
               elseif ((dip .eq. 0.) .and. (stk .eq. 0.)) then
                   if (j .eq. 1) then
                      angle(j) = 0.
                   else
                      angle(j) = angle(j-1)
                   endif
               endif
               angle_p(j) = (green_mesh%rake_i(green_mesh%vnorm_i(i))*(pi/180.)) - angle(j)
!print *,  (green_mesh%rake_i(green_mesh%vnorm_i(i))*(180./pi)), angle(j)*180./pi
               ang_vx(j) = ( 1. / ((dip**2. / stk) + stk ))
               ang_vy(j) = -1. * ( dip / (dip**2. + stk**2.) )
               !green_mesh%rake_p(cont) = ang_vx(j)*angle_p(j)
               !green_mesh%rake_p(cont+green_mesh%interp_i) = ang_vy(j) *angle_p(j)
               if (ind(j) .eq. 1) then
                green_mesh%rake_p(cont) = ang_vx(j)*angle_p(j)
                green_mesh%rake_p(cont+green_mesh%interp_i) = ang_vy(j) *angle_p(j)
                cost_ang = cost_ang + angle_p(j)**2.
               endif
               write(iunit,*) green_mesh%rake_p(cont), green_mesh%rake_p(cont+green_mesh%interp_i), angle(j)*180./pi, ind(j)
             !  write(iunit,*) ang_vx(j), green_mesh%rake_p(cont), ang_vy(j), green_mesh%rake_p(cont+57)
               j = j + 1
               cont = cont + 1
            enddo
        enddo

        close(iunit)

        cost_ang = (cost_ang / 2.)

        print *, cost_ang, 'cost angle'
        print *, green_mesh%costa
        green_mesh%costa = green_mesh%costa + 0.001 * cost_ang
        print *, green_mesh%costa
!stop

        deallocate(ind)
        deallocate(ang_vx,ang_vy)
        deallocate(angle,angle_p)
      endsubroutine rake_angle


! --------------------------------------------------------------------
! SUBROUTINE  MeanVariance():
!    This subroutine computes the mean, variance and standard
! deviation.
! --------------------------------------------------------------------

   SUBROUTINE  MeanVariance(Data, SIZE, n, Mean, Variance, StdDev)
      IMPLICIT  NONE
      INTEGER, INTENT(IN)                 :: SIZE
      INTEGER, INTENT(IN)                 :: n
      REAL, DIMENSION(1:SIZE), INTENT(IN) :: Data
      REAL, INTENT(OUT)                   :: Mean, Variance, StdDev
      INTEGER                             :: i

      Mean = 0.0
      DO i = 1, n
         Mean = Mean + Data(i)
      END DO
      Mean = Mean / n

      Variance = 0.0
      DO i = 1, n
         Variance = Variance + (Data(i) - Mean)**2
      END DO
      Variance = Variance / n
      StdDev   = SQRT(Variance)
   END SUBROUTINE  MeanVariance



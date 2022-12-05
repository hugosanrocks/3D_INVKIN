

     ! Regularization terms controling the rake angle
     ! using unitary vectors.

      subroutine rake_uni(green_mesh)

      implicit none
      INCLUDE 'green.h'
!     optim type strcuture
      TYPE (mesh) :: green_mesh

      ! Arrays to save rake angle and rake angle prior
      real, dimension(:), allocatable :: angle_r, angle_p, angle_q, ang_vx, ang_vy
      integer, dimension(:), allocatable :: ind
      real :: stk, dip, rak, modu, vector(2)
      integer i, j, k, cont, n, ini, fin
      integer :: iunit, ii, cont2
      real,dimension(:),allocatable :: serie, serie2, serie3
      real pi2, cost_ang, trigg, cof, maxtrace, coef, val
      real term_1, term_2, term_3, term_4, x, z, x2, z2, x3, z3, grad(green_mesh%modelsize2)
      real pi/3.1415926535897932/
      pi2 = pi * 2.

      cof = 0.01!0.01
      print *, ' Measuring the rake angle'
 
      !samples, subfaults
      !Arrays to store the rake angles
      allocate(angle_r(green_mesh%interp_i*2*green_mesh%msub))
      allocate(angle_q(green_mesh%interp_i*2*green_mesh%msub))
      allocate(angle_p(green_mesh%interp_i*2*green_mesh%msub))
      allocate(ang_vx(green_mesh%interp_i*green_mesh%msub))
      allocate(ang_vy(green_mesh%interp_i*green_mesh%msub))
      allocate(ind(green_mesh%interp_i*green_mesh%msub))
      allocate(serie(green_mesh%interp_i))
      allocate(serie2(green_mesh%msub))
      allocate(serie3(green_mesh%interp_i*green_mesh%msub))

      n = green_mesh%interp_i
      angle_r(:) = 0.
      angle_p(:) = 0.
      angle_q(:) = 0.
      ang_vx(:) = 0.

      ini = 1
      fin = green_mesh%interp_i

      iunit=32
      open(iunit,file=green_mesh%out//'angle.out',status='unknown')
      open(71,file=green_mesh%out//'anglep.out',status='unknown')

      !From 1D scalar to a 2D vector
      if (green_mesh%rake_opt .eq. 1) then
      call model_d(green_mesh%model,green_mesh%model2, &
  &                green_mesh%interp_i,green_mesh%msub,&
  &                green_mesh%slipm,green_mesh%dir_i,&
  &                green_mesh%vnorm_i)
      endif

      !Global array for regularization
      green_mesh%rake_p(:) = 0.

      serie(:) = 0.
      serie2(:) = 0.
      serie3(:) = 0.

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
        serie3(j) = modu
        serie(k) = modu
        !2D rake angle of current solution
        !Avoiding NaN when modulus = 0.
        !When modu = 0. set the angle to the prior information
        if (modu .eq. 0.) then
         angle_r(cont) = green_mesh%vslip2(green_mesh%vnorm_i(i),1)
         angle_r(cont+green_mesh%interp_i) = &
  &                      green_mesh%vslip2(green_mesh%vnorm_i(i),2)
        else
         angle_r(cont) = stk / modu
         angle_r(cont+green_mesh%interp_i) = dip / modu
        endif
        !2D rake vector of angle from prior knowledge
        angle_p(cont) = green_mesh%vslip2(green_mesh%vnorm_i(i),1)
        angle_p(cont+green_mesh%interp_i) = & 
  &                     green_mesh%vslip2(green_mesh%vnorm_i(i),2)
        j = j + 1
        cont = cont + 1
       enddo
       serie2(i) = maxval(serie)
      enddo

      j = 1
      do i=1,green_mesh%msub
       !trigg = maxval(green_mesh%model2)*0.2         !serie2(i)*cof
       trigg = serie2(i)*cof
       !Select one subfault at each time
       cont=1+(i-1)*green_mesh%interp_i*2
       do k=1,green_mesh%interp_i
        if (serie3(j) .ge. (trigg)) then
         ind(j) = 1
         angle_q(cont) = angle_r(cont) - angle_p(cont)
         angle_q(cont+green_mesh%interp_i) = &
  &        angle_r(cont+green_mesh%interp_i) - &
  &        angle_p(cont+green_mesh%interp_i)
        write(71,*) angle_p(cont), angle_r(cont), angle_q(cont), 'vx'
        write(71,*) angle_p(cont+green_mesh%interp_i), angle_r(cont+green_mesh%interp_i), angle_q(cont+green_mesh%interp_i), 'vz'
         !Interprete the angle
         if ( angle_r(cont+green_mesh%interp_i) .gt. 0.) then
          if ( angle_r(cont) .gt. 0.) then
           ang_vx(j) = (pi - atan(angle_r(cont+green_mesh%interp_i)/angle_r(cont)))
          elseif ( angle_r(cont) .lt. 0.) then
           ang_vx(j) =  (atan(-1.*angle_r(cont+green_mesh%interp_i)/angle_r(cont)))
          endif
         elseif ( angle_r(cont+green_mesh%interp_i) .lt. 0.) then
          if ( angle_r(cont) .gt. 0.) then
           ang_vx(j) = (pi + atan(-1.*angle_r(cont+green_mesh%interp_i)/angle_r(cont)))
          elseif ( angle_r(cont) .lt. 0.) then
           ang_vx(j) = (pi*2. - atan(angle_r(cont+green_mesh%interp_i)/angle_r(cont)))
          endif
         elseif (( angle_r(cont+green_mesh%interp_i) .eq. 0.) .and. &
  &              (angle_r(cont) .eq. 0.)) then
         if (j .eq. 1) then
          ang_vx(j) = 0.
         else
           ang_vx(j) = ang_vx(j-1)
         endif
         endif
        else
         ind(j) = 0
        endif
        write(iunit,*) j, serie3(j), ind(j), ang_vx(j)*180./pi
        val = ang_vx(j)*180./pi
        !Hard boundaries at allowed rake limits
        !Discontinuity case, around 360 - 0
        if ( ind(j) .eq. 1 ) then
         if ( green_mesh%rak_case .eq. 2) then
          if ( (val .ge. green_mesh%rake_lim(1,2)) .and. &
  &            (val .le. green_mesh%rake_lim(1,3)) ) then
             !Inside this limit do nothing
          !Lower limit
          elseif (( val .lt. green_mesh%rake_lim(1,2)) .and. & 
  &               ( val .gt. green_mesh%rake_lim(1,3)) ) then
            dip = green_mesh%dip_i(green_mesh%vnorm_i(i))
            stk = green_mesh%stk_i(green_mesh%vnorm_i(i))
            rak = green_mesh%rake_lim(green_mesh%vnorm_i(i),2)
            call vector_transformation(stk,dip,rak,vector)
            green_mesh%model2(cont) = serie3(j)*vector(1)
            green_mesh%model2(cont+green_mesh%interp_i) = serie3(j)*vector(2)
            !Can be also set to zero
            !green_mesh%model2(cont) = 0.
            !green_mesh%model2(cont+green_mesh%interp_i) = 0.
          !Upper limit
          elseif (( val .gt. green_mesh%rake_lim(1,3)) .and. &
  &               ( val .lt. green_mesh%rake_lim(1,2)) ) then
            dip = green_mesh%dip_i(green_mesh%vnorm_i(i))
            stk = green_mesh%stk_i(green_mesh%vnorm_i(i))
            rak = green_mesh%rake_lim(green_mesh%vnorm_i(i),3)
            call vector_transformation(stk,dip,rak,vector)
            green_mesh%model2(cont) = serie3(j)*vector(1)
            green_mesh%model2(cont+green_mesh%interp_i) = serie3(j)*vector(2)
            !Can be also set to zero
            !green_mesh%model2(cont) = 0.
            !green_mesh%model2(cont+green_mesh%interp_i) = 0.
          endif
         elseif ( green_mesh%rak_case .eq. 1) then
          !Hard boundaries at allowed rake limits
          if ( (val .le. green_mesh%rake_lim(1,2)) .and. &
  &          (val .ge. green_mesh%rake_lim(1,3)) ) then
             !Inside this limit do nothing
          !Upper limit NORMAL CASE
          elseif ( val .gt. green_mesh%rake_lim(1,2)) then
            !Lower limit
            dip = green_mesh%dip_i(green_mesh%vnorm_i(i))
            stk = green_mesh%stk_i(green_mesh%vnorm_i(i))
            rak = green_mesh%rake_lim(green_mesh%vnorm_i(i),2)
            call vector_transformation(stk,dip,rak,vector)
            green_mesh%model2(cont) = serie3(j)*vector(1)
            green_mesh%model2(cont+green_mesh%interp_i) = serie3(j)*vector(2)
          elseif ( val .lt. green_mesh%rake_lim(1,3)) then
          !Upper limit
            dip = green_mesh%dip_i(green_mesh%vnorm_i(i))
            stk = green_mesh%stk_i(green_mesh%vnorm_i(i))
            rak = green_mesh%rake_lim(green_mesh%vnorm_i(i),3)
            call vector_transformation(stk,dip,rak,vector)
            green_mesh%model2(cont) = serie3(j)*vector(1)
            green_mesh%model2(cont+green_mesh%interp_i) = serie3(j)*vector(2)
          endif
         endif
        else
        endif
        j = j + 1
        cont = cont + 1
       enddo
      enddo
      close(iunit)
      close(71)


      !Estimate gradient
      grad(:) = 0.
      cost_ang = 0.
      j = 1
      do i=1,green_mesh%msub
       trigg = serie2(i)*cof
       !Select one subfault at each time
       cont = 1+(i-1)*green_mesh%interp_i*2
       do k=1,green_mesh%interp_i
!        if (( angle_r(cont+green_mesh%interp_i) .gt. 0. ) .and. &
!   &        ( angle_r(cont) .gt. 0.)) then
         if ( serie3(j) .ge. trigg ) then
         x = angle_r(cont)
         z = angle_r(cont+green_mesh%interp_I) 
         x2 = 1 / ( (z**2. / x) + x  )
         z2 = -1.*z / ( z**2. + x**2.)
         x3 = atan(angle_p(cont+green_mesh%interp_i)/angle_p(cont))
         cost_ang = cost_ang + (atan(z/x) - x3)**2.
         grad(cont) = (atan(z/x) - x3  )*x2
         grad(cont+green_mesh%interp_i) = (atan(z/x)-x3 )*z2
        endif 
        j = j + 1
        cont = cont + 1
       enddo
      enddo

      green_mesh%gradad(:) = grad(:)
      green_mesh%costm = cost_ang / 2.

      deallocate(serie,serie2,serie3)
      deallocate(ind)
      deallocate(ang_vx,ang_vy)
      deallocate(angle_r,angle_p,angle_q)
      endsubroutine rake_uni


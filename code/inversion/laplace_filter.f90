       subroutine init_laplace_filter(green_mesh)

       IMPLICIT NONE
       INCLUDE 'green.h'
       TYPE(mesh) :: green_mesh

       integer :: lenx, lenz, lent, i, j, k, p, q
       real, dimension(:), allocatable :: xc, zc, tc
       real, dimension(:,:,:), allocatable :: filter1, filter2
       real :: dist, factor, alpha, beta, lx, lz, lt, h
       real pi/3.1415926535897932/
       real :: theta1, theta2, theta, phi, varphi, x2, z2, t2
       real :: a, b, c, d, e, f, g, hh, ii, matrix(3,3), vecin(3), vecout(3)
       real :: a2, b2, c2, d2, e2, f2, g2, hh2, ii2, x3, z3, t3, theta3
       integer iunit, contf, contf2, cut1, cut2, cut3, cut4
       real spike(21,9,50), y(41,17,99)
       integer lenx2, lenz2, lent2, lenxy, lenzy, lenty
       real resul(21,9,50), v(3)
       integer i2, j2, k2, cont
       integer ii1, iii2, ii3, ii4, ini1, ini2

       v(1) = -1./sqrt(2.)
       v(2) = 1./sqrt(2.)
       v(3) = 0.

       iunit = 77
       open(iunit,file=green_mesh%dat//'laplace.dat',status='old',&
  &         action='read')
       read(iunit,*) lx, lz, lt
       read(iunit,*) lenx, lenz, lent
       read(iunit,*) theta, phi, varphi
       close(iunit)

       !left side
       theta1 = -1.*theta*pi/180.
       !right side
       theta2 = theta*pi/180.
       phi    = phi*pi/180.
       varphi = varphi*pi/180.

       theta3 = 0.
       !rotation on X-Z strike-dip plane
       a2=  cos(phi)*cos(varphi) + sin(phi)*sin(theta3)*sin(varphi)
       b2=  cos(phi)*sin(theta3)*sin(varphi) -1.*cos(varphi)*sin(phi)
       c2=  cos(theta3)*sin(varphi)
       d2=  cos(theta3)*sin(phi)
       e2=  cos(phi)*cos(theta3)
       f2= -1.*sin(theta3)
       g2=  cos(varphi)*sin(phi)*sin(theta3) -1.*cos(phi)*sin(varphi)
       hh2= sin(phi)*sin(varphi) + cos(phi)*cos(varphi)*sin(theta3)
       ii2= cos(theta3)*cos(varphi)

! Arbitrary rotation on time axis
!       a = cos(theta) + v(1)**2*(1.-cos(theta))
!       b = v(1)*v(2)*(1.-cos(theta))-1.*v(3)*sin(theta)
!       c = v(1)*v(3)*(1.-cos(theta)) + v(2)*sin(theta)
!       d = v(1)*v(2)*(1.-cos(theta)) + v(3)*sin(theta)
!       e = cos(theta) + v(2)**2*(1.-cos(theta))
!       f = v(2)*v(3)*(1.-cos(theta)) -1.*v(1)*sin(theta)
!       g = v(3)*v(1)*(1.-cos(theta)) -1.*v(2)*sin(theta)
!       hh = v(3)*v(2)*(1.-cos(theta))+v(1)*sin(theta)
!       ii = cos(theta) + v(3)**2*(1.-cos(theta))

       h = 1.!green_mesh%stk_s
       alpha = lx/lz
       beta  = lx/lt
       factor = (alpha*beta)/((8.*pi)*(lx/h)**3)


       allocate(xc(lenx),zc(lenz),tc(lent))
       allocate(filter1(lenx,lenz,lent))
       allocate(filter2(lenx,lenz,lent))

       do i=1,lenx
         xc(i) = real(-1*10+(i-1))*h
       enddo
       ! print *, xc
       do i=1,lenz
         zc(i) = real(-1*4+(i-1))*h
       enddo
       ! print *, zc
       do i=1,lent
         tc(i) = real(-1*25+(i-1))*h
       enddo
       ! print *, tc

       lenx2 = 21
       lenz2 = 9
       lent2 = 50
       spike(:,:,:) = 0.
       spike(7,3,1) = 1.
!       spike(6,2,6) = 1.
!       spike(6,2,7) = 0.5
!       spike(6,2,8) = 0.25
!       spike(6,2,9) = 0.12
!       spike(6,2,10) = 0.06
!       spike(6,2,11) = 0.03

!       spike(6,5,5) = 2.
!       spike(6,5,6) = 1.
!       spike(6,5,7) = 0.5
!       spike(6,5,8) = 0.25
!       spike(6,5,9) = 0.12
!       spike(6,5,10) = 0.06
!       spike(6,5,11) = 0.03



       lenxy = lenx + lenx2 - 1
       lenzy = lenz + lenz2 - 1
       lenty = lent + lent2 - 1

       filter1(:,:,:) = 0.
       filter2(:,:,:) = 0.
       do contf = 1, 1
        if (contf .eq. 1) then
          theta = theta1
        else
          theta = theta2
        endif
        do k=1,lent
         do j=1,lenz
          do i=1,lenx
           a = cos(theta) + v(1)**2*(1.-cos(theta))
           b = v(1)*v(2)*(1.-cos(theta))-1.*v(3)*sin(theta)
           c = v(1)*v(3)*(1.-cos(theta)) + v(2)*sin(theta)
           d = v(1)*v(2)*(1.-cos(theta)) + v(3)*sin(theta)
           e = cos(theta) + v(2)**2*(1.-cos(theta))
           f = v(2)*v(3)*(1.-cos(theta)) -1.*v(1)*sin(theta)
           g = v(3)*v(1)*(1.-cos(theta)) -1.*v(2)*sin(theta)
           hh = v(3)*v(2)*(1.-cos(theta))+v(1)*sin(theta)
           ii = cos(theta) + v(3)**2*(1.-cos(theta))
           x2 = xc(i)*a + zc(j)*b + tc(k)*c
           z2 = xc(i)*d + zc(j)*e + tc(k)*f
           t2 = xc(i)*g + zc(j)*hh + tc(k)*ii
           x3 = x2*a2 + z2*b2 + t2*c2
           z3 = x2*d2 + z2*e2 + t2*f2
           t3 = x2*g2 + z2*hh2+ t2*ii2
           dist = sqrt((x3)**2 + (alpha**2*(z3**2))+ &
  &                  (beta**2*(t3**2)))
           if (contf .eq. 1) then
            filter1(i,j,k) = factor * exp(-1.*(dist/lx))
           else
            filter2(i,j,k) = factor * exp(-1.*(dist/lx))
           endif
          enddo
         enddo
        enddo
        if (contf .eq. 1) then
         call conv3d(spike,lenx2,lenz2,lent2,filter1,lenx,lenz,lent,y,lenxy,lenzy,lenty)
        else
         call conv3d(spike,lenx2,lenz2,lent2,filter2,lenx,lenz,lent,y,lenxy,lenzy,lenty)
        endif
        if (contf .eq. 1) then
         k2=1
         do ii1=26,26+49
          j2=1
          do iii2=5,5+8 !8
           i2=1
           do ii3=11,31
             resul(i2,j2,k2) = y(ii3,iii2,ii1)
             i2 = i2+1
           enddo
           j2 = j2+1
          enddo
          k2 = k2+1
         enddo
        else
         k2=1
         do ii1=26,26+49
          j2=1
          do iii2=5,5+8
           i2=12
           do ii3=22,31
             resul(i2,j2,k2) = y(ii3,iii2,ii1)
             i2 = i2+1
           enddo
           j2 = j2+1
          enddo
          k2 = k2+1
         enddo
        endif

      enddo

        do ii1=1,50
        do iii2=1,9
          do ii3=1,21
             write(29,*) resul(ii3,iii2,ii1)
          enddo
        enddo
       enddo

      ! do k=1,lent
      ! do j=1,lenz
      !    do i=1,lenx
      !       write(26,*) filter(i,j,k)
      !    enddo
      !  enddo
      ! enddo

       !do k=1,lent2
       !do j=1,lenz2
       !   do i=1,lenx2
       !      write(27,*) spike(i,j,k)
       !   enddo
       ! enddo
       !enddo

       !do k=1,lenty
       !do j=1,lenzy
       !   do i=1,lenxy
       !      write(28,*) y(i,j,k)
       !   enddo
       ! enddo
       !enddo

       deallocate(xc,zc,tc)
       deallocate(filter1,filter2)
       endsubroutine init_laplace_filter

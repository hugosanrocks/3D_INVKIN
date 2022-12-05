       subroutine laplace_matrix(green_mesh)

       IMPLICIT NONE
       INCLUDE 'green.h'
       TYPE(mesh) :: green_mesh

       integer :: lenx, lenz, lent, i, j, k, ii, jj, kk, iii
       real, dimension(:), allocatable :: xc, zc, tc, lx, lz, lt
       real, dimension(:), allocatable :: theta, phi
       real, dimension(:), allocatable :: filter, vector, result, result2, res
       real :: dist, factor, alpha, beta, lxm, lzm, ltm, h
       real :: xm, zm, tm
       real :: a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, vec1(3), vec2(3), vec(3)
       integer :: lxs, lzs, lts, cont1, cont2
       real :: x, z, t, val
       integer rot_opt, pp
       real pi/3.1415926535897932/

       rot_opt = 1

       lxs = 21
       lzs = 21
       lts = 50

       allocate(xc(lxs),zc(lzs),tc(lts))
       allocate(lx(lxs*lzs*lts),lz(lxs*lzs*lts),lt(lxs*lzs*lts))
       allocate(filter(lxs*lzs*lts))
       allocate(vector(lxs*lzs*lts))
       allocate(result(lxs*lzs*lts))
       allocate(result2(lxs*lzs*lts))
       allocate(res(lxs*lzs*lts))
       allocate(theta(lxs*lzs*lts))
       allocate(phi(lxs*lzs*lts))


       do i=1, lxs
        xc(i) = -10. + real(i-1)
       enddo
       do i=1, lzs
        zc(i) = -10. + real(i-1)
       enddo
       do i=1, lts
        tc(i) = -25. + real(i-1)
       enddo
       lx(:) = 1.
       lz(:) = 1.
       lt(:) = 1.
       !Angles to rotate in X-Z plane
       theta(:) = 0.
       !Angle to rotate in Time
       phi(:) = 0.
       !Reference vector to rotate in time
       vec1(:) = [-1./sqrt(2.),  1./sqrt(2.), 0.]
       vec2(:) = [ 1./sqrt(2.), -1./sqrt(2.), 0.]

       h = 1.

       result(:) = 0.
       vector(:) = 0.
       i=198 !121
       vector(i) = 1.
       lx(i) = 4.
       lz(i) = 4.
       lt(i) = 1.
       theta(i) = 0.
       i=198 !341
       vector(i) = 1.
       lx(i) = 4.
       lz(i) = 4.
       lt(i) = 1.
       theta(i) = 0.
       i=198 !121
       vector(i) = 1.
       lx(i) = 4.
       lz(i) = 4.
       lt(i) = 1.
       theta(i) = 0.
       i=198 !341
       vector(i) = 1.
       lx(i) = 4.
       lz(i) = 4.
       lt(i) = 1.
       theta(i) = 0.

!       vector(204) = 0.
!       vector(244) = 1.
!       vector(245) = 1.
!       lx(204) = 1.
!       lz(204) = 1.
!       lx(244) = 4.
!       lz(244) = 4.
!       lx(245) = 4.
!       lz(245) = 4.
!       lx(i) = 4.
!       lz(i) = 4.
!       lt(i) = 1.

!      lx(244) = 4.
!      lz(244) = 4.
!      lt(244) = 1.
!       theta(67) = 0.*pi/180.
!       theta(47) = 0.*pi/180.

       cont1 = 1
       do k=1,lts
         do j=1,lzs
           do i=1,lxs
            if ( (i .lt. 11) ) then  !bien
             phi(cont1) = 0.*pi/180.
            elseif  (i .ge. 11 ) then !bien
             phi(cont1) = 0.*pi/180.
            endif
            cont1 = cont1 + 1
           enddo
         enddo
       enddo

phi(:) = 0.

       cont2 = 1
       do kk = 1, lts
        tm = tc(kk)
        do jj = 1, lzs
         zm = zc(jj)
         do ii = 1, lxs
           xm = xc(ii)
           cont1 = 1
           do k = 1, lts
            do j = 1, lzs
             do i = 1, lxs
              vec(:) = vec1(:)
              if ( (ii .lt. 11) .and. (jj .lt. 11) ) then
               theta(cont1) = 0.
               phi(cont1) = -70.*pi/180.
              elseif ( (ii .ge. 11) .and. (jj .lt. 11) ) then
               theta(cont1) = 0.5*pi
               !j1 = cos(2.*theta(cont1))
               !k1 = -1.*sin(2.*theta(cont1))
               phi(cont1) = -70.*pi/180.
              elseif ( (ii .lt. 11) .and. (jj .ge. 11) ) then
               theta(cont1) = 0.5*pi
               !j1 = cos(2.*theta(cont1))
               !k1 = -1.*sin(2.*theta(cont1))
               phi(cont1) = 70.*pi/180.
              elseif ( (ii .ge. 11) .and. (jj .ge. 11) ) then
               theta(cont1) = 0.
               !j1 = cos(theta(cont1))
               !k1 = sin(theta(cont1))
               !vec(:) = vec1(:)
               phi(cont1) = 70.*pi/180.
              endif
               a1 = cos(phi(cont1)) + vec(1)**2*(1.-cos(phi(cont1)))
               b1 = vec(1)*vec(2)*(1.-cos(phi(cont1)))-1.*vec(3)*sin(phi(cont1))
               c1 = vec(1)*vec(3)*(1.-cos(phi(cont1))) + vec(2)*sin(phi(cont1))
               d1 = vec(1)*vec(2)*(1.-cos(phi(cont1))) + vec(3)*sin(phi(cont1))
               e1 = cos(phi(cont1)) + vec(2)**2*(1.-cos(phi(cont1)))
               f1 = vec(2)*vec(3)*(1.-cos(phi(cont1))) -1.*vec(1)*sin(phi(cont1))
               g1 = vec(3)*vec(1)*(1.-cos(phi(cont1))) -1.*vec(2)*sin(phi(cont1))
               h1 = vec(3)*vec(2)*(1.-cos(phi(cont1)))+vec(1)*sin(phi(cont1))
               i1 = cos(phi(cont1)) + vec(3)**2*(1.-cos(phi(cont1)))
              if (( (i .lt. 11) .and. (j .lt. 11) ) .or. &
  &               ( (i .ge. 11) .and. (j .ge. 11) )) then
               j1 = cos(theta(cont1))
               k1 = sin(theta(cont1))
               x = (xc(i)*j1 -1.*zc(j)*k1)*a1 + (xc(i)*k1 + zc(j)*j1)*b1 & 
  &                + tc(k)*c1
               z = (xc(i)*j1 -1.*zc(j)*k1)*d1 + (xc(i)*k1 + zc(j)*j1)*e1 & 
  &                + tc(k)*f1
               t = (xc(i)*j1 -1.*zc(j)*k1)*g1 + (xc(i)*k1 + zc(j)*j1)*h1 & 
  &                + tc(k)*i1
               xm =  (xc(ii)*j1 -1.*zc(jj)*k1)*a1 + (xc(ii)*k1 + zc(jj)*j1)*b1 &
  &                + tc(kk)*c1
               zm =  (xc(ii)*j1 -1.*zc(jj)*k1)*d1 + (xc(ii)*k1 + zc(jj)*j1)*e1 &
  &                + tc(kk)*f1
               tm =  (xc(ii)*j1 -1.*zc(jj)*k1)*g1 + (xc(ii)*k1 + zc(jj)*j1)*h1 &
  &                + tc(kk)*i1
               elseif (( (i .ge. 11) .and. (j .lt. 11) ) .or. &
  &                    ( (i .lt. 11) .and. (j .ge. 11) )) then
               theta(cont1) = 0.5*pi
               j1 = cos(2.*theta(cont1))
               k1 = -1.*sin(2.*theta(cont1))
               x = (xc(i)*j1 -1.*zc(j)*k1)*a1 + (-1.*xc(i)*k1 -1.*zc(j)*j1)*b1 &
  &                + tc(k)*c1
               z = (xc(i)*j1 -1.*zc(j)*k1)*d1 + (-1.*xc(i)*k1 -1.*zc(j)*j1)*e1 &
  &                + tc(k)*f1
               t = (xc(i)*j1 -1.*zc(j)*k1)*g1 + (-1.*xc(i)*k1 -1.*zc(j)*j1)*h1 &
  &                + tc(k)*i1
               xm =  (xc(ii)*j1 -1.*zc(jj)*k1)*a1 + (-1.*xc(ii)*k1 -1.*zc(jj)*j1)*b1 &
  &                + tc(kk)*c1
               zm =  (xc(ii)*j1 -1.*zc(jj)*k1)*d1 + (-1.*xc(ii)*k1 -1.*zc(jj)*j1)*e1 &
  &                + tc(kk)*f1
               tm =  (xc(ii)*j1 -1.*zc(jj)*k1)*g1 + (-1.*xc(ii)*k1 -1.*zc(jj)*j1)*h1 &
  &                + tc(kk)*i1
               endif
               lxm = lx(cont1)
               lzm = lz(cont1)
               ltm = lt(cont1)
               alpha = lxm / lzm
               beta = lxm / ltm
               factor = (alpha*beta)/(8.*pi*(lxm/h)**3)
               dist = sqrt((x-xm)**2. &
  &                       +alpha**2*(z-zm)**2 &
  &                       +beta**2*(t-tm)**2 )
               filter(cont1) = factor*exp(-1.*dist/lxm)
               cont1 = cont1 + 1
             enddo
            enddo
           enddo
          val = 0.
          do iii=1,lxs*lzs*lts
            result(cont2) = result(cont2) + filter(iii)*vector(iii)
           enddo
         cont2 = cont2 + 1 
         enddo
        enddo
       enddo

       do i=1,lxs*lzs*lts
        write(30,*) result(i)
        write(311,*) vector(i)
       enddo

       deallocate(theta,phi)
       deallocate(filter,vector,result,result2,res)
       deallocate(xc,zc,tc,lx,lz,lt)
       endsubroutine laplace_matrix

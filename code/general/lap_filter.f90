       subroutine lap_filter(green_mesh)

       IMPLICIT NONE
       INCLUDE 'green.h'
       TYPE(mesh) :: green_mesh

       integer :: i, j, k, ii, jj, kk, iii
       real, dimension(:), allocatable :: xc, zc, tc, lx, lz, lt
       real, dimension(:), allocatable :: theta, phi
       real, dimension(:), allocatable :: filter, result, vector
       real :: dist, factor, alpha, beta, lxm, lzm, ltm, h
       real :: xm, zm, tm
       real :: a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, vec1(3), vec2(3), vec(3)
       integer :: lxs, lzs, lts, cont1, cont2
       real :: x, z, t
       Integer :: rot_opt
       real pi/3.1415926535897932/

       rot_opt = 1

       lxs = 35
       lzs = 22
       lts = 39


       allocate(xc(lxs),zc(lzs),tc(lts))
       allocate(lx(lxs*lzs*lts),lz(lxs*lzs*lts),lt(lxs*lzs*lts))
       allocate(filter(lxs*lzs*lts))
       allocate(vector(lxs*lzs*lts))
       allocate(result(lxs*lzs*lts))
       allocate(theta(lxs*lzs*lts))
       allocate(phi(lxs*lzs*lts))

       !coordinates
       do i=1, lxs
        xc(i) = 1. + real(i-1)
       enddo
       do i=1, lzs
        zc(i) = 1. + real(i-1)
       enddo
       do i=1, lts
        tc(i) = 0. + real(i-1)*0.25
       enddo

       !correlation lengths       
       lx(:) = 3.
       lz(:) = 2.
       lt(:) = 0.25

       !Angles to rotate in X-Z plane
       theta(:) = 10.
       !Angle to rotate in Time
       phi(:) = 0.
       !Reference vector to rotate in time
       vec1(:) = [-1./sqrt(2.),  1./sqrt(2.), 0.]
       vec2(:) = [ 1./sqrt(2.), -1./sqrt(2.), 0.]

       h = 1.

       !output and input vectors
       result(:) = 0.

       ii=648
       k=1
       do i=1,648
        do j=1,41
         jj = (j-1)*ii + i
         vector(jj) = green_mesh%grad1(k)
         k=k+1
        enddo
       enddo

       !i=342+10368 !121
       !lx(i) = 4.
       !lz(i) = 2.
       !lt(i) = 1.
       !print *, 'theta: phi:'
       !theta(i) = 45.
       !read(*,*) theta(i), phi(i)
       theta(:) = theta(:)*(pi/180.)
       phi(:) = phi(:)*(pi/180.)

       !Reference center (tm,zm,xm)
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
               !generalized rotation with respect to a vector (vec)
               a1 = cos(phi(cont1)) + vec(1)**2*(1.-cos(phi(cont1)))
               b1 = vec(1)*vec(2)*(1.-cos(phi(cont1)))-1.*vec(3)*sin(phi(cont1))
               c1 = vec(1)*vec(3)*(1.-cos(phi(cont1))) + vec(2)*sin(phi(cont1))
               d1 = vec(1)*vec(2)*(1.-cos(phi(cont1))) + vec(3)*sin(phi(cont1))
               e1 = cos(phi(cont1)) + vec(2)**2*(1.-cos(phi(cont1)))
               f1 = vec(2)*vec(3)*(1.-cos(phi(cont1))) -1.*vec(1)*sin(phi(cont1))
               g1 = vec(3)*vec(1)*(1.-cos(phi(cont1))) -1.*vec(2)*sin(phi(cont1))
               h1 = vec(3)*vec(2)*(1.-cos(phi(cont1)))+vec(1)*sin(phi(cont1))
               i1 = cos(phi(cont1)) + vec(3)**2*(1.-cos(phi(cont1)))
               !rotation around x-z
               j1 = cos(theta(cont1))
               k1 = sin(theta(cont1))
               !rotate coordinates and reference coordinates
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
               !correlation distances
               lxm = lx(cont1)
               lzm = lz(cont1)
               ltm = lt(cont1)
               alpha = lxm / lzm
               beta = lxm / ltm
               !normalization factor
               factor = (alpha*beta)/(8.*pi*(lxm/h)**3)
               !distance
               dist = sqrt((x-xm)**2. &
  &                       +alpha**2*(z-zm)**2 &
  &                       +beta**2*(t-tm)**2 )
               !filter
               filter(cont1) = factor*exp(-1.*dist/lxm)
               cont1 = cont1 + 1
             enddo
            enddo
           enddo
          !cumulate multiplication
          do iii=1,lxs*lzs*lts
            result(cont2) = result(cont2) + filter(iii)*vector(iii)
          enddo
         cont2 = cont2 + 1 
         enddo
        enddo
       enddo

       ii=648
       k=1
       do i=1,648
        do j=1,41
         jj = (j-1)*ii + i
         green_mesh%grad1(k) = result(jj)
         k=k+1
        enddo
       enddo


       deallocate(theta,phi)
       deallocate(filter,result,vector)
       deallocate(xc,zc,tc,lx,lz,lt)
       endsubroutine lap_filter


       !subroutine to filter only 2D slices
       !of a 3D volume

       subroutine lap_filter_2d(green_mesh)

       IMPLICIT NONE
       INCLUDE 'green.h'
       TYPE(mesh) :: green_mesh

       integer :: i, j, k, ii, jj, kk, iii, jjj, comp
       real, dimension(:), allocatable :: xc, zc, lx, lz
       real, dimension(:), allocatable :: theta
       real, dimension(:), allocatable :: filter, result, vector
       real :: dist, factor, alpha, lxm, lzm, h
       real :: xm, zm
       real :: j1, k1
       integer :: lxs, lzs, lts, cont1, cont2, col1, col2
       real :: x, z
       integer :: xhyp, zhyp, idhyp
       real pi/3.1415926535897932/


       !dimensions
       lxs = green_mesh%nsstk
       lzs = green_mesh%nsdip
       lts = green_mesh%interp_i

       !allocate memory
       allocate(xc(lxs),zc(lzs))
       allocate(lx(lxs*lzs),lz(lxs*lzs))
       allocate(filter(lxs*lzs))
       allocate(vector(lxs*lzs))
       allocate(result(lxs*lzs))
       allocate(theta(lxs*lzs))

       !coordinates
       do i=1, lxs
        xc(i) = 1. + real(i-1)
       enddo
       do i=1, lzs
        zc(i) = 1. + real(i-1)
       enddo

       !correlation lengths       
       lx(:) = real(green_mesh%r2)
       lz(:) = real(green_mesh%r1)

       !step of grid
       h = green_mesh%stk_s

      !node with the hypocenter
      !general definition of hypo
      xhyp = green_mesh%xhyp
      zhyp = green_mesh%zhyp
      idhyp = lxs*(green_mesh%zhyp-1)+ xhyp
      print *, xhyp, zhyp, idhyp, 'hypo'
      !toy
      !xhyp = 1
      !zhyp = 1
      !idhyp = 1
      !SIV1
      !xhyp = 26
      !zhyp = 11
      !idhyp=422
      !422
      !JAP
      !xhyp = 7
      !zhyp = 7
      !idhyp= 181
      !181
      !Concentric smoothing approach near hypocenter
      !Angles to rotate in X-Z plane
      k=1 
      do i=1,green_mesh%nsdip
        do j=1,green_mesh%nsstk
        !siv1 0.6
        !jap 0.1
         if ( ( i .le. zhyp ) .and. ( j .le. xhyp ) ) then
          theta(k) = (1.-green_mesh%distcorr(k,idhyp)/(maxval(green_mesh%distcorr(:,idhyp)))) * (40.)
         elseif ( ( i .gt. zhyp ) .and. ( j .le. xhyp ) ) then
          theta(k) = (1.-green_mesh%distcorr(k,idhyp)/(maxval(green_mesh%distcorr(:,idhyp)))) * (-40.)
         elseif ( ( i .le. zhyp ) .and. ( j .gt. xhyp ) ) then
          theta(k) = (1.-green_mesh%distcorr(k,idhyp)/(maxval(green_mesh%distcorr(:,idhyp)))) * (-40.)
         elseif ( ( i .gt. zhyp ) .and. ( j .gt. xhyp ) ) then
          theta(k) = (1.- green_mesh%distcorr(k,idhyp)/(maxval(green_mesh%distcorr(:,idhyp)))) * (40.)
         endif
         k=k+1
        enddo
       enddo

       !theta(:) = 85.
       theta(:) = theta(:)*(pi/180.)

       !open(60,file='grad_nofil.ascii',status='unknown')
       !open(61,file='grad_fil.ascii',status='unknown')

      if ( green_mesh%rake_opt .eq. 1 ) then
       comp = 1
      else
       comp = 2
      endif

      do jjj=1,comp                !2D vector if necessary
      
      do kk=1,green_mesh%interp_i

       !output vectors
       result(:) = 0.

       !input vector
       !2D slices for each time sample
       col1 = 1 + (jjj-1)*green_mesh%msub
       col2 = col1+green_mesh%msub-1
       vector(:) = green_mesh%slipr2(kk,col1:col2)

       !only to prove the response of the filter
       !vector(:) = 0.
       !vector(216+400) = 1.
       !vector(220+400) = 1.
       !vector(224+400) = 1.
       !vector(230+400) = 1.
       !vector(235+400) = 1.

       !vector(132) = 1.
       !vector(83) = 1.
       !vector(213) = 1.
       !vector(208) = 1.
       !theta(177) = -30.*pi/180.

       !Reference center (tm,zm,xm)
       cont2 = 1
        do jj = 1, lzs
         do ii = 1, lxs
           cont1 = 1
            do j = 1, lzs
             do i = 1, lxs
               !rotation around x-z
               j1 = cos(theta(cont1))
               k1 = sin(theta(cont1))
               !rotate coordinates and reference coordinates
               x = xc(i)*j1 + zc(j)*k1
               z = -1.*xc(i)*k1 + zc(j)*j1
               xm= xc(ii)*j1 + zc(jj)*k1
               zm= -1.*xc(ii)*k1 + zc(jj)*j1
               !correlation distances
               lxm = lx(cont1)
               lzm = lz(cont1)
               alpha = lxm / lzm
               !normalization factor
               factor = (alpha)/(2.*pi*(lxm/h)**2.)
               !distance
               dist = sqrt( (x-xm)**2. + alpha**2*(z-zm)**2 )
               !filter
               filter(cont1) = factor*exp(-1.*dist/lxm)
               cont1 = cont1 + 1
             enddo
            enddo
          !cumulate multiplication
          do iii=1,lxs*lzs
            result(cont2) = result(cont2) + filter(iii)*vector(iii)
          enddo
         cont2 = cont2 + 1 
         enddo
        enddo
        !assign output
        green_mesh%slipr2(kk,col1:col2) = result(:)
       enddo
       !only to check the effect of the filter
       !do i=1,green_mesh%msub
       ! write(60,*) vector(i)
       ! write(61,*) result(i)
       !enddo
       !stop

       enddo

       !close(60)
       !close(61)

       deallocate(theta)
       deallocate(filter,result,vector)
       deallocate(xc,zc,lx,lz)
       endsubroutine lap_filter_2d

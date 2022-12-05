

     !This program estimates the position of nodes
     !discretizing a finite fault plane.

     !INPUT
     !strike        angle
     !dip           angle
     !length_stike  km
     !lenght_dip    km

     !OUTPUT
     !coor          positions of nodes (x,y,z) km
                    !hypocenter is taken as (0,0,0)

     program ffault_discretization

     
      implicit none
      real         hypo(3)
      real strike1, dip1, angle1, strike2, dip2, angle2, dip
      real len_strike1, len_dip1, len_strike2, len_dip2
      real len_strike3, len_dip3, len_strike4, len_dip4
      real         dx, dy, dz, d_str, d_dip, area, di
      real         dist, x, y, stkr, dipr
      real         corners(6,3), vstk(3), vdip(3), vnorm(3)
      real         disty
      integer i, j, k, nstk1, ndip1, nstk2, ndip2
      real,dimension(:,:),allocatable :: coor, coor2
      real pi/3.1415926535897932/

     hypo = [0., 0., 12.45]

     strike1 = 205.
     dip1 = 72.
     strike2 = 235.
     dip2 = 65.

     !sizes to the right and left of hypocenter

     !WEST                        EAST      UP
     !------------------!/   s        
     !       dip1       !/    t      /dip3
     !  stk1  *   stk2  !/     k    /
     !       dip2       !/      e  /
     !------------------!/        /dip4
     !                            /         DOWN

     !1st plane
     len_strike1 = 6.*2.+1.   !Distance west of hypo
     len_dip1 = 5.*2.+1.      !Distance above hypo
     len_strike2 = 2.*2.+1.   !Distance east of hypo
     len_dip2 = 3.*2.+1.      !Distance below hypo
     !2nd plane
     len_strike3 = 20.*2.+1
     len_dip3 = 9.*2.+1

     nstk1 = 9
     ndip1 = 9
     nstk2 = 20
     ndip2 = 9

     d_str = 2.               !Node separation
     d_dip = 2.               !Node separation

     angle1 = strike1
     call distaz2xy(len_strike1,angle1,x,y)
     print *, x, y, 'cor1'
     corners(1,:) = [x, y, hypo(3)]
     angle1 = strike1 - 180.
     call distaz2xy(len_strike2,angle1,x,y)
     !print *, x, y
     corners(2,:) = [x, y, hypo(3)]

     stkr = atan(y/x)
     dipr = dip1*pi/180.
     print *, stkr*180./pi, 'stkr', dipr*180./pi, 'dipr'
     print *, 'Remember cartesian and geog coordinate system'

     vstk(1) = -1.*cos(stkr)
     vstk(2) = -1.*sin(stkr)
     vstk(3) = 0.
     print *, vstk, 'vstk'

     vnorm(1)= -1.*sin(dipr)*sin(stkr)  !*-1 original
     vnorm(2)= sin(dipr)*cos(stkr)
     vnorm(3)= -1.*cos(dipr)
     print *, vnorm, 'vnorm'

     vdip(1)= vstk(2)*vnorm(3)-vnorm(2)*vstk(3)
     vdip(2)= -1.*(vstk(1)*vnorm(3)-vnorm(1)*vstk(3))
     vdip(3)= vstk(1)*vnorm(2)-vnorm(1)*vstk(2)
     print *, vdip, 'vdip'

     corners(3,:) = corners(1,:) +(len_dip1*vdip(:))
     corners(4,:) = corners(2,:) +(len_dip1*vdip(:))
     corners(5,:) = corners(1,:) -(len_dip2*vdip(:))
     corners(6,:) = corners(2,:) -(len_dip2*vdip(:))

      !Write to check corners of the 1st plane
      do i=1,6
       print *, corners(i,:), i
      enddo

     di = d_str / 2.

     allocate(coor((nstk1+nstk2),3),coor2((nstk1+nstk2)*ndip1,3))

     coor(1,:) = corners(3,:) -1.*(di*vstk(:)) -1.*(di*vdip(:))
       do i=2,nstk1
        coor(i,:) = coor(i-1,:) - (d_str*vstk(:))
       ! print *, coor(i,:)
       enddo
      ! print *, vdip, 'dip'
       k=1
       do j=1,ndip1
       do i=1,nstk1
        coor2(k,:) = coor(i,:) -1.*(((real(j)-1.)*d_dip)*vdip(:))
     !   print *, coor2(k,:), i, j
        k = k+1
       enddo
       enddo

      open(22,file='fault.dat',status='unknown')
      do i=1,nstk1*ndip1
        write(22,*) coor2(i,1:3)
      enddo
      close(22)


!       print *, corners(3,:) - coor2(1,:)
!       print *, corners(4,:) - coor2(7,:)
!       print *, corners(5,:) - coor2(7*8+1,:)
!       print *, corners(6,:) - coor2(7*9,:)


     !Now the origin is at the last node along strike
     !of the first plane and at hypo depth
     angle2 = strike2
     call distaz2xy(len_strike3,angle2,x,y)

     stkr = atan(y/x)
     dipr = dip2*pi/180.
!     print *, stkr*180./pi, 'stk2', dipr*180./pi, 'dip2'

     vstk(1) = -1.*cos(stkr)
     vstk(2) = -1.*sin(stkr)
     vstk(3) = 0.
!     print *, vstk, 'vstk'

     vnorm(1)= -1.*sin(dipr)*sin(stkr)  !*-1 original
     vnorm(2)= sin(dipr)*cos(stkr)
     vnorm(3)= -1.*cos(dipr)
!     print *, vnorm, 'vnorm'

     vdip(1)= vstk(2)*vnorm(3)-vnorm(2)*vstk(3)
     vdip(2)= -1.*(vstk(1)*vnorm(3)-vnorm(1)*vstk(3))
     vdip(3)= vstk(1)*vnorm(2)-vnorm(1)*vstk(2)
!     print *, vdip, 'vdip'

     !Take the shallowest right corner from 1st plane and
     !connect it to the second plane
     corners(3,:) = corners(4,:)
     corners(4,:) = corners(3,:) + [-1.*x, -1.*y, 0.]
     corners(5,:) = corners(3,:) -1.*(len_dip3*vdip)
     corners(6,:) = corners(4,:) -1.*(len_dip3*vdip)
      !Write to check corners of the 1st plane
      do i=1,6
!       print *, corners(i,:), i
      enddo

     coor(1,:) = corners(3,:) -1.*(di*vstk(:)) -1.*(di*vdip(:))
       do i=2,nstk2
        coor(i,:) = coor(i-1,:) - (d_str*vstk(:))
       ! print *, coor(i,:)
       enddo
       ! print *, vdip, 'dip'
       do j=1,ndip2
       do i=1,nstk2
        coor2(k,:) = coor(i,:) -1.*(((real(j)-1.)*d_dip)*vdip(:))
     !   print *, coor2(k,:), i, j
        k = k+1
       enddo
       enddo

      open(22,file='fault.dat2',status='unknown')
      do i=1,(nstk1+nstk2)*ndip1
        write(22,*) coor2(i,1:3)
      enddo
      close(22)



     deallocate(coor,coor2)
     endprogram ffault_discretization



      !This subroutine converts the polar coordinates (r,az) of 
      !a station, which are given with respect to an epicenter,
      ! location to Cartesian coordinates (x,y).

      subroutine distaz2xy(dist,az,x,y)

      IMPLICIT NONE
      integer i, j
      real,intent(inout) :: dist, az, x, y
      real ang
      real pi/3.1415926535897932/

      x=0.
      y=0.

       !1st case, N-E quadrant
       if ( (az .ge. 0.).and. (az .le. 90.) ) then
         ang = az*pi/180.
         x = dist*sin(ang)
         y = dist*cos(ang)
       !2nd case, S-E quadrant
       elseif ( (az .gt. 90.) .and. (az .le. 180.) ) then
         ang = (az-90.)*pi/180.
         x = dist*cos(ang)
         y = (-1.*dist*sin(ang))
       elseif ( (az .gt. 180.) .and. (az .le. 270.) ) then
         ang = (az-180.)*pi/180.
         x = (-1.*dist*sin(ang))
         y = (-1.*dist*cos(ang))
       elseif ( (az .gt. 270.) .and. (az .le. 360.)  ) then
         ang = (az-270.)*pi/180.
         x = (-1.*dist*cos(ang))
         y = dist*sin(ang)
       endif


      endsubroutine distaz2xy




      subroutine check_angle(az)

      IMPLICIT NONE
      integer i, j
      real,intent(inout) :: az
      real dist, x, y, ang
      real pi/3.1415926535897932/

      x=0.
      y=0.

       !1st case, N-E quadrant
       if ( (az .ge. 0.).and. (az .le. 90.) ) then
         ang = az*pi/180.
         x = dist*sin(ang)
         y = dist*cos(ang)
       !2nd case, S-E quadrant
       elseif ( (az .gt. 90.) .and. (az .le. 180.) ) then
         ang = (az-90.)*pi/180.
         x = dist*cos(ang)
         y = (-1.*dist*sin(ang))
       elseif ( (az .gt. 180.) .and. (az .le. 270.) ) then
         ang = (az-180.)*pi/180.
         x = (-1.*dist*sin(ang))
         y = (-1.*dist*cos(ang))
       elseif ( (az .gt. 270.) .and. (az .le. 360.)  ) then
         ang = (az-270.)*pi/180.
         x = (-1.*dist*cos(ang))
         y = dist*sin(ang)
       endif


      endsubroutine check_angle


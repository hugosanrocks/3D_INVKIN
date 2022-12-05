

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
      real         strike, dip, len_strike, len_dip, angle
      real         dx, dy, dz, d_str, d_dip, area, di
      real         dist, x, y, strike2, stkr, dipr
      real         corners(6,3), vstk(3), vdip(3), vnorm(3)
      integer i, j, k
      real,dimension(:,:),allocatable :: coor, coor2
      real pi/3.1415926535897932/

     hypo = [0., 0., 20.]

     strike = 289.
     dip = 12.
     len_strike = 19.*4.5
     len_dip = 19.*4.5
     d_str = 4.5
     d_dip = 4.5
     area = len_strike * len_dip
     dist = len_strike / 2.

     call distaz2xy(dist,strike,x,y)
     corners(1,:) = [x, y, hypo(3)]
     corners(2,:) = -1.*corners(1,:)
     corners(2,3) = hypo(3)

     stkr = asin(y/x)
     dipr = dip*pi/180.
     
     vstk(1) = cos(stkr)
     vstk(2) = sin(stkr)
     vstk(3) = 0.
     
     vnorm(1)= -1.*sin(dipr)*sin(stkr)
     vnorm(2)= sin(dipr)*cos(stkr)
     vnorm(3)= -1.*cos(dipr)

     vdip(1)= vstk(2)*vnorm(3)-vnorm(2)*vstk(3)
     vdip(2)= -1.*(vstk(1)*vnorm(3)-vnorm(1)*vstk(3))
     vdip(3)= vstk(1)*vnorm(2)-vnorm(1)*vstk(2)

     corners(3,:) = corners(1,:) +(dist*vdip(:))
     corners(4,:) = corners(2,:) +(dist*vdip(:))
     corners(5,:) = corners(1,:) -(dist*vdip(:))
     corners(6,:) = corners(2,:) -(dist*vdip(:))
      do i=1,6
       print *, corners(i,:)
      enddo

     vstk(1) = x/sqrt(x**2+y**2)
     vstk(2) = y/sqrt(x**2+y**2)
     vstk(3) = 0.

     dx = dist*2. / 18.
     di = dx/2.

     allocate(coor(18,3),coor2(18*18,3))

     coor(1,:) = corners(5,:) - (di*vstk(:)) + (di*vdip(:))
     print *, dx, di
       do i=2,18
        coor(i,:) = coor(i-1,:) - (dx*vstk(:))
       ! print *, coor(i,:)
       enddo
       print *, vdip, 'dip'
       k=1
       do j=1,18
       do i=1,18
        coor2(k,:) = coor(i,:) + (((real(j)-1)*dx)*vdip(:))
        print *, coor2(k,:), i, j
        k = k+1
       enddo
       enddo

      open(22,file='fault.dat',status='unknown')
      do i=1,18*18
        write(22,*) coor2(i,1:3)
      enddo
      close(22)


       print *, corners(5,:) - coor2(1,:)
       print *, corners(6,:) - coor2(18,:)
       print *, corners(3,:) - coor2(18*17+1,:)
       print *, corners(4,:) - coor2(18*18,:)

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

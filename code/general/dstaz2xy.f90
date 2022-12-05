
      !This program converts the polar coordinates (r,az) of 
      !a station, which are given with respect to an epicenter,
      !location to Cartesian coordinates (x,y).

      program distaz2xy

      IMPLICIT NONE
      integer nsta, i, j, iunitin, iunitout
      real dist, az, ang
      real x, y
      character*4 sta
      real pi/3.1415926535897932/

      iunitout = 33

      open(iunitout,file='stations.xy',status='unknown',&
  &        action='write')


      iunitin=22
      open(iunitin,file='stations.dstaz',status='old',&
  &        action='read')
      read(iunitin,*) nsta

      do i=1,nsta
       read(iunitin,*) sta, dist, az
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
       write(iunitout,*) sta, x, y
      enddo

      close(iunitin)
      close(iunitout)

      endprogram distaz2xy

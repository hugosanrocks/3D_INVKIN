      module triangles

       contains

       !=============================================!
       !function triangle_area(x1,x2,x3,y1,y2,y3)    !
       !function returning the area of any triangle  !
       !based on the coordinates of the three corners!
       !=============================================!

       function triangle_area(x_in,y_in)

       implicit none
       real, dimension(3) :: x_in, y_in
       real :: triangle_area

        x_in(:) = x_in(:)*1000.
        y_in(:) = y_in(:)*1000.
        triangle_area=abs(x_in(1)*(y_in(2)-y_in(3)) + &
  &                       x_in(2)*(y_in(3)-y_in(1)) + &
                          x_in(3)*(y_in(1)-y_in(2))) / 2.
       
       return
       endfunction triangle_area


       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       !function triangle_transform
       !Transforms the coordinates of an array of points coordinates
       !(xw(siz,1:2) inside any type of triangular surface 
       !(x_in(3),y_in(3)) into the equivalent coordinates 
       !inside a standard scalene triangle.
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

       function triangle_transform(x_in,y_in,xw,siz)

       implicit none
       integer :: siz
       real :: x_in(3), y_in(3)
       real, dimension(siz,3) :: xw
       real, dimension(siz,2) :: triangle_transform
       real :: x(siz), y(siz)
       integer :: i

       !transform all the points
       do i=1,siz
        x(i) = x_in(1)*(1-xw(i,1)-xw(i,2))+ &
  &            x_in(2)*xw(i,1)+x_in(3)*xw(i,2)
        y(i) = y_in(1)*(1-xw(i,1)-xw(i,2))+ &
  &            y_in(2)*xw(i,1)+y_in(3)*xw(i,2)
        triangle_transform(i,1:2) = [x(i), y(i)]
       enddo

       return
       endfunction triangle_transform


       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       !function integrate_triangle
       !Computes the 2D surface integral of a function 
       !that lives inside a standard scalene triangle.
       !input: val       values of the function at the nodes
       !                 of the triangle element.
       !       x_in y_in coordinates of irregualr triangle element
       !       xw        coordinates of gauss points and weights
       !       siz       number of desired gauss points
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

       function integrate_triangle(pint,x_in,y_in,xw,val,siz)

       implicit none
       integer :: siz
       real, dimension(siz,2) :: pint
       real, dimension(siz) :: val
       real, dimension(siz,3) :: xw
       real, dimension(3) :: x_in, y_in
       real :: area, val_int
       real :: integrate_triangle
       integer :: i

       !weighted sum
       integrate_triangle = 0.
       do i=1,siz
        val_int = triangle_interpol(pint(i,1),pint(i,2),x_in,y_in,val)
        integrate_triangle = integrate_triangle + val_int*xw(i,3)
       enddo

       !normalization to original triangle area
       area = triangle_area(x_in,y_in)
       integrate_triangle = area * integrate_triangle

       return
       endfunction integrate_triangle


       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       !Function triangle_interpol
       !Interpolates the value of a function inside a standard 
       !scalene triangle given the values of the function at 
       !the three corners of the scalene triangle.
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

       function triangle_interpol(x_in,y_in,x,y,val)

       implicit none
       real, dimension(3) :: val
       real, dimension(3) :: x, y
       real :: x_in, y_in
       real :: w1, w2, w3
       real :: triangle_interpol

       !estimate weights
       w1 = ((y(2) - y(3))*(x_in-x(3))+(x(3)-x(2))*(y_in-y(3))) / &
  &         ((y(2) - y(3))*(x(1)-x(3))+(x(3)-x(2))*(y(1)-y(3)))
       w2 = ((y(3) - y(1))*(x_in-x(3))+(x(1)-x(3))*(y_in-y(3))) / &
  &         ((y(2) - y(3))*(x(1)-x(3))+(x(3)-x(2))*(y(1)-y(3)))
       w3 = 1. - w1 - w2

       !get the value of interpolation
       triangle_interpol = w1*val(1) + w2*val(2) + w3*val(3)

       return
       endfunction triangle_interpol


       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       !   Function triangle_gauss_points provides the Gaussian     !
       !   points and weights for the Gaussian quadrature of        !
       !   order n for the standard triangles.                      !
       !                                                            !
       !   Input: n   - the order of the Gaussian quadrature (n<=7) !
       !                                                            !
       !   Output: xw - a n by 3, matrix:                            !
       !       1st column gives the x-coordinates of points         !
       !       2nd column gives the y-coordinates of points         !
       !       3rd column gives the weights                         !
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

       function triangle_gauss_points(n,siz)

       implicit none
       integer :: n, siz
       real, dimension(siz,3) :: triangle_gauss_points
       real, dimension(siz,3) :: xw

       xw(:,3) = 0.
       if (n .eq. 1) then
         xw(1,:)=[0.33333333333333,    0.33333333333333,    1.00000000000000]
       elseif (n .eq. 2) then
         xw(1,:)=[   0.,    1./2.,    1./3.]
         xw(2,:)=[1./2.,       0.,    1./3.]
         xw(3,:)=[1./2.,    1./2.,    1./3.]
       elseif (n .eq. 3) then
         xw(1,:)=[0.33333333333333,    0.33333333333333,   -0.56250000000000]
         xw(2,:)=[0.20000000000000,    0.20000000000000,    0.52083333333333]
         xw(3,:)=[0.20000000000000,    0.60000000000000,    0.52083333333333]
         xw(4,:)=[0.60000000000000,    0.20000000000000,    0.52083333333333]
       elseif (n .eq. 4) then
         xw(1,:)=[0.44594849091597,    0.44594849091597,    0.22338158967801]
         xw(2,:)=[0.44594849091597,    0.10810301816807,    0.22338158967801]
         xw(3,:)=[0.10810301816807,    0.44594849091597,    0.22338158967801]
         xw(4,:)=[0.09157621350977,    0.09157621350977,    0.10995174365532]
         xw(5,:)=[0.09157621350977,    0.81684757298046,    0.10995174365532]
         xw(6,:)=[0.81684757298046,    0.09157621350977,    0.10995174365532]
       elseif (n .eq. 5) then
         xw(1,:)=[0.33333333333333,    0.33333333333333,    0.22500000000000]
         xw(2,:)=[0.47014206410511,    0.47014206410511,    0.13239415278851]
         xw(3,:)=[0.47014206410511,    0.05971587178977,    0.13239415278851]
         xw(4,:)=[0.05971587178977,    0.47014206410511,    0.13239415278851]
         xw(5,:)=[0.10128650732346,    0.10128650732346,    0.12593918054483]
         xw(6,:)=[0.10128650732346,    0.79742698535309,    0.12593918054483]
         xw(7,:)=[0.79742698535309,    0.10128650732346,    0.12593918054483]
       elseif (n .eq. 6) then
         xw(1,:)=[0.24928674517091,    0.24928674517091,    0.11678627572638]
         xw(2,:)=[0.24928674517091,    0.50142650965818,    0.11678627572638]
         xw(3,:)=[0.50142650965818,    0.24928674517091,    0.11678627572638]
         xw(4,:)=[0.06308901449150,    0.06308901449150,    0.05084490637021]
         xw(5,:)=[0.06308901449150,    0.87382197101700,    0.05084490637021]
         xw(6,:)=[0.87382197101700,    0.06308901449150,    0.05084490637021]
         xw(7,:)=[0.31035245103378,    0.63650249912140,    0.08285107561837]
         xw(8,:)=[0.63650249912140,    0.05314504984482,    0.08285107561837]
         xw(9,:)=[0.05314504984482,    0.31035245103378,    0.08285107561837]
         xw(10,:)=[0.63650249912140,    0.31035245103378,    0.08285107561837]
         xw(11,:)=[0.31035245103378,    0.05314504984482,    0.08285107561837]
         xw(12,:)=[0.05314504984482,    0.63650249912140,    0.08285107561837]
       elseif (n .eq. 7) then
         xw(1,:)=[0.33333333333333,    0.33333333333333,   -0.14957004446768]
         xw(2,:)=[0.26034596607904,    0.26034596607904,    0.17561525743321]
         xw(3,:)=[0.26034596607904,    0.47930806784192,    0.17561525743321]
         xw(4,:)=[0.47930806784192,    0.26034596607904,    0.17561525743321]
         xw(5,:)=[0.06513010290222,    0.06513010290222,    0.05334723560884]
         xw(6,:)=[0.06513010290222,    0.86973979419557,    0.05334723560884]
         xw(7,:)=[0.86973979419557,    0.06513010290222,    0.05334723560884]
         xw(8,:)=[0.31286549600487,    0.63844418856981,    0.07711376089026]
         xw(9,:)=[0.63844418856981,    0.04869031542532,    0.07711376089026]
         xw(10,:)=[0.04869031542532,    0.31286549600487,    0.07711376089026]
         xw(11,:)=[0.63844418856981,    0.31286549600487,    0.07711376089026]
         xw(12,:)=[0.31286549600487,    0.04869031542532,    0.07711376089026]
         xw(13,:)=[0.04869031542532,    0.63844418856981,    0.07711376089026]
       elseif (n .eq. 8) then
         xw(1,:)=[0.33333333333333,    0.33333333333333,    0.14431560767779]
         xw(2,:)=[0.45929258829272,    0.45929258829272,    0.09509163426728]
         xw(3,:)=[0.45929258829272,    0.08141482341455,    0.09509163426728]
         xw(4,:)=[0.08141482341455,    0.45929258829272,    0.09509163426728]
         xw(5,:)=[0.17056930775176,    0.17056930775176,    0.10321737053472]
         xw(6,:)=[0.17056930775176,    0.65886138449648,    0.10321737053472]
         xw(7,:)=[0.65886138449648,    0.17056930775176,    0.10321737053472]
         xw(8,:)=[0.05054722831703,    0.05054722831703,    0.03245849762320]
         xw(9,:)=[0.05054722831703,    0.89890554336594,    0.03245849762320]
         xw(10,:)=[0.89890554336594,    0.05054722831703,    0.03245849762320]
         xw(11,:)=[0.26311282963464,    0.72849239295540,    0.02723031417443]
         xw(12,:)=[0.72849239295540,    0.00839477740996,    0.02723031417443]
         xw(13,:)=[0.00839477740996,    0.26311282963464,    0.02723031417443]
         xw(14,:)=[0.72849239295540,    0.26311282963464,    0.02723031417443]
         xw(15,:)=[0.26311282963464,    0.00839477740996,    0.02723031417443]
         xw(16,:)=[0.00839477740996,    0.72849239295540,    0.02723031417443]
       else
        write(*,*) 'Bad input n for triangle_gauss_points'
       endif

       !return the Gauss points and weights
       triangle_gauss_points = xw
       return
       endfunction triangle_gauss_points

       !=======================================================!
       ! function point_in_triangle:                           !
       !                                                       !
       ! Function to determine if a point                      !
       ! (px,py) lives inside of a given 2D triangle           !
       ! which node coordinates are given by                   !
       ! p0x,p0y,p1x,p1y,p2x,p2y -> ordered conterclockwise    !
       !                                                       !
       ! if the point is inside point_in_triang = 1
       ! if not point_in_triangle = 0
       !=======================================================!

       function point_in_triang(px,py,p0x,p0y,p1x,p1y,p2x,p2y)

        implicit none
        integer :: point_in_triang
        real :: s, t, px, py
        real :: p0x, p0y, p1x, p1y, p2x, p2y
        real :: area, d1, d2, d3, det1, det2, det3
        real :: dp1, dp2, dp3

        ! Evaluate s, t and 1-s-t. 
        ! The point (px,py) is inside the triangle if and
        ! only if they are all positive.

        !signed area of triangle
        Area = 0.5*(-1.*p1y*p2x + p0y*(-1.*p1x + p2x) + &
  &            p0x*(p1y - p2y) + p1x*p2y)

        !variable to account for inside or not test
        s = 1./(2.*Area)*(p0y*p2x - p0x*p2y + &
  &         (p2y - p0y)*px + (p0x - p2x)*py)
        t = 1./(2.*Area)*(p0x*p1y - p0y*p1x + &
  &         (p0y - p1y)*px + (p1x - p0x)*py)

        !edge distances used to know if point is at edges
        dp1 = sqrt((p0x-p1x)**2.+(p0y-p1y)**2.)
        dp2 = sqrt((p1x-p2x)**2.+(p1y-p2y)**2.)
        dp3 = sqrt((p0x-p2x)**2.+(p0y-p2y)**2.)

        !distance from (px,py) to triangle nodes
        d1 = sqrt((px-p0x)**2.+(py-p0y)**2.)
        d2 = sqrt((px-p1x)**2.+(py-p1y)**2.)
        d3 = sqrt((px-p2x)**2.+(py-p2y)**2.)

        !determinants used for point at edge test
        det2 = px*(p1y-p2y)-py*(p1x-p2x)+(p1x*p2y-p2x*p1y)
        det1 = px*(p0y-p1y)-py*(p0x-p1x)+(p0x*p1y-p1x*p0y)
        det3 = px*(p0y-p2y)-py*(p0x-p2x)+(p0x*p2y-p2x*p0y)

         !test to pass
         !inside or not test
         if ( ( s .gt. 0. ) .and. ( t .gt. 0. ) &
            .and. ( (1.-s-t ) .gt. 0. ) ) then
          point_in_triang = 1
         else
          !if point is outside
          point_in_triang = 0
          !if the point is at the node (vertex) positions
          if ( (d1 .eq. 0.) .or. (d2 .eq. 0.) .or. (d3 .eq. 0.) ) then
           point_in_triang = 1
          !if point is at any of three edges
          elseif ( ((det1 .eq. 0.).and. ((d1 .lt. dp1) .and. (d2 .lt. dp1))) .or. &
  &                ((det2 .eq. 0.).and. ((d2 .lt. dp2) .and. (d2 .lt. dp3))) .or. &
  &                ((det3 .eq. 0.).and. ((d3 .lt. dp3) .and. (d1 .lt. dp3))) ) then
           point_in_triang = 1
          endif
        endif

       endfunction


      end module triangles

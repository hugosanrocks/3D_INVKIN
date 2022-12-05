      program test_triangles_module

      !Program used to test "triangles.f90" module
      !built to compute the 2D integral of a function
      !discretized by random triangles

      use triangles
      implicit none
      real :: x1, x2, x3, y1, y2, y3, x, y
      !to test area
      real :: area
      !to test gauss points
      real,dimension(:,:),allocatable :: xw
      integer :: n, siz
      !to test traingle transformation
      real :: x_in(3), y_in(3)
      real,dimension(:,:),allocatable :: val_out
      !to test bilinear interpolation
      real :: val_int, val_nint(3)
      !to test surface integration
      real :: val(3)
      integer :: i

      !coordinates of random triangle
      x1 = 0.
      x2 = 1.
      x3 = 0.
      y1 = 0.
      y2 = 0.
      y3 = 1.
      x_in = [x1, x2, x3]
      y_in = [y1, y2, y3]

      !test function of area
      area = triangle_area(x_in,y_in)
      !print *, 'Area:', area

      !test function of gauss points
      n = 2
      siz = 3
      allocate(xw(siz,3))
      xw = triangle_gauss_points(n,siz)
      !print *, 'Gauss points:'
      do i=1,siz
      ! print *, xw(i,:)
      enddo

      !test transformation to unitary standard triangle
      allocate(val_out(siz,2))
      val_out = triangle_transform(x_in,y_in,xw,siz)
      !print *, 'Points in standard triangle:'
      do i=1,siz
      ! print *, val_out(i,:)
      enddo

      !interpolation of function at that point
      !values of function at the nodes
      val_nint = [1., 2., 0.]
      !nodes of triangle
      xw(1,:) = [0.5, 1., 1.]
      xw(2,:) = [0., 0.,  1.]
      xw(3,:) = [1., 0., 1.]
      print *, '                       x           y            val'
      do i=1,3
       print *,'                ', xw(i,1), xw(i,2), val_nint(i)
      enddo
      val_int = triangle_interpol(xw(1,1),xw(1,2),xw(:,1),xw(:,2),val_nint)
      print *, 'Interpolated:', xw(1,1),xw(1,2),val_int
      val_int = triangle_interpol(xw(2,1),xw(2,2),xw(:,1),xw(:,2),val_nint)
      print *, 'Interpolated:', xw(2,1),xw(2,2),val_int
      val_int = triangle_interpol(xw(3,1),xw(3,2),xw(:,1),xw(:,2),val_nint)
      print *, 'Interpolated:', xw(3,1),xw(3,2),val_int
      x=0.5
      y=0.31
      val_int = triangle_interpol(x,y,xw(:,1),xw(:,2),val_nint)
      print *, 'Interpolated:', x,y,val_int


      !integration of the function at that triangle
      !val(:) = [1., 1., 1.]
      !do i=1,100*500
       
      ! area = integrate_triangle(x_in,y_in,xw,val,siz)
      !enddo
      !print *, 'Intregal inside surface:', area

      !deallocate(xw,val_out)
      endprogram test_triangles_module

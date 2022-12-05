       program test

       implicit none
       real array(5), x(5), value, xwant
       real array2d(5,5), y(5), value2, ywant
       integer nx, ny, i, j
       real xint(9), yint(9), arrayint(9,9), delta
       nx = 5
       ny = 5
       do i= 1,5
          x(i) = real(i)
          y(i) = real(i)
       enddo

       xint(1) = 1.
       yint(1) = 1.
       delta = 0.5
       do i= 2,9
          xint(i) = xint(i-1) + delta
          yint(i) = yint(i-1) + delta
       enddo

       do j=1,nx
       do i=1,nx
         array2d(i,j) = x(i)**2+y(j)**2
       enddo
       enddo

       do i= 1,9
       do j= 1,9
       xwant = xint(j)
       ywant = yint(i)
       call FIND2D(Xwant,Ywant,Value,X,Y,Array2d,nx,ny)
       write(44,*) value
       enddo
       enddo

       


       endprogram test

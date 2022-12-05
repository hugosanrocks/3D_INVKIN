        program test

        implicit none
        real x(5), y(5), xint(12), yint(12)
        integer nx, ny, nxint, nyint
        real mat(5,5), matint(12,12), xw, yw, res
        integer i, j, k

        nx = 5
        ny = 5

        nxint=12
        nyint=12

        do i=1,nx
         x(i) = real(i)
         do j=1,ny
           mat(i,j) = real(i*j)
           y(j) = real(j)
         enddo
         print *, mat(i,:)
        enddo
        print *,'x',x
        print *,'y',y

        do i=1,nxint
         xint(i) = real(i)/2.
         yint(i) = real(i)/2.
        enddo
        print *, 'xint', xint
        print *, 'yint', yint

        do i=1,nxint
         do j=1,nyint
          xw=xint(i)
          yw=yint(j)
          call FIND2D(xw,yw,res,x,y,mat,nx,ny)
          matint(i,j) = res
          print *, xw, yw, res
         enddo
        enddo




        endprogram test

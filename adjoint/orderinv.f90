       subroutine order_inv(x,lx)

       implicit none
       integer, intent(inout) :: lx
       real, intent(inout) :: x(lx)
       real :: y(lx)
       integer i,j


       !X is the dummy input array to invert in 1D (mirror)

       !Inversion in 1D
       j=lx
       DO i=1,lx
         y(i)=x(j)
       j=j-1
       ENDDO

       !Assign the inverted vector to the output array
       x(:)=y(:)

       end subroutine order_inv

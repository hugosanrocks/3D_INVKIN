       subroutine res_order(x,h,y,lx)

       implicit none
       integer, intent(inout) :: lx
       real, intent(inout) :: x(lx), h(lx), y(lx)
       integer i,j


       j=lx
       DO i=1,LX
         Y(i)=x(j)-h(i)
       j=j-1
       ENDDO

       end subroutine res_order

       subroutine conv(x,lx,h,lh,y,ly)

       !1D CONVOLUTION SUBROUTINE

       implicit none
       !Dimension of vectors
       integer, intent(inout) :: LX, LH, LY
       !Vectors to be convolved
       real, intent(inout) :: x(lx), h(lh), y(ly)
       !Counters
       integer i,j

       y(:) = 0.

       DO i=1,LX
        DO j=1,Lh
         Y(i+j-1)=Y(i+j-1)+X(i)*h(j)
        ENDDO
       ENDDO


       end subroutine conv

       subroutine conv2d(x,mx,lx,h,mh,lh,y,my,ly)

       IMPLICIT NONE
       !Dimension of vectors
       integer, intent(inout) :: lx, mx,lh, mh, ly, my
       !Vectors to be convolved
       real, intent(inout) :: x(mx,lx), h(mh,lh), y(my,ly)
       !Counters
       integer i, j, k, p

       y(:,:) = 0.
       do k = 1, mx
        do p = 1, mh
         do i = 1, lx
          do j = 1, lh
           y(k+p-1,i+j-1) = y(k+p-1,i+j-1)+ &
  &                         x(k,i)*h(p,j)
          enddo
         enddo
        enddo
       enddo

       endsubroutine conv2d


       subroutine conv3d(x,mx,lx,nx,h,mh,lh,nh,y,my,ly,ny)

       IMPLICIT NONE
       !Dimension of vectors
       integer, intent(inout) :: mx, lx, nx, mh, lh, nh, my, ly, ny
       !Vectors to be convolved
       real, intent(inout) :: x(mx,lx,nx), h(mh,lh,nh), y(my,ly,ny)
       !Counters
       integer i, j, k, p, q, r

       y(:,:,:) = 0.
       do r = 1, nx
        do q = 1, nh
         do k = 1, mx
          do p = 1, mh
           do i = 1, lx
            do j = 1, lh
             y(k+p-1,i+j-1,r+q-1) = y(k+p-1,i+j-1,r+q-1)+ &
  &                           x(k,i,r)*h(p,j,q)
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo

       endsubroutine conv3d


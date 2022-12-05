      subroutine xcorr1d(x,lx,h,lh,y,ly)

      !1D cross correlation of two vectors
      
      !INPUT
      !x, h     vectors to correlate
      !lx, lh   length of vectors
      !ly       length of the output vector
      !         ly = lx + lh - 1

      !OUTPUT
      !ly       vector to store correlation

      integer,intent(in) :: lx, lh, ly
      real,intent(inout) :: x(lx), h(lh), y(ly)
      integer :: i, j
      real,dimension(:),allocatable :: l

      allocate(l(lh))

      !Flip flop signal for correlation
      j = lh
      do i=1,lh
       l(i)=h(j)
       j = j - 1
      enddo

      !1D cross correlation
      do i=1,lx
       do j=0,lh-1
        y(i+j)=y(i+j)+x(i)*l(j+1)
       enddo
      enddo

      deallocate(l)
      endsubroutine xcorr1d

      subroutine fcost(r,ly,lx,coss)

      !COST function subroutine (L2)

      implicit none

      !Variables needed only by this program
      INTEGER, INTENT(INOUT) :: ly, lx
      REAL, INTENT(INOUT) :: r(ly,lx)
      REAL, INTENT(INOUT) :: coss
      REAL cost(lx,lx), normv(lx), y(ly,lx), norm(5)
      integer i, j, k

      cost(:,:) = 0.

      !NORMALIZE THE CONTRIBUTION OF EACH STATION
      !normv = maxval(r,dim=1)
      !j=1
      !do i=1,5
       !norm(i) = maxval(normv(j:j+2))
       !do k=0,2
       !y(:,j+k) = r(:,j+k) / norm(i)
       !enddo
      !j=j+3
      !enddo

      y(:,:) = r(:,:)

      !Stack the sqaure of the residuals
      do i=1,lx
        do j=1,lx
         do k=1,ly
          cost(i,j)=cost(i,j)+ y(k,j)*y(k,i)
         enddo
        enddo
      enddo

      coss=0.
      do i=1,lx
       coss = coss + cost(i,i)
      enddo


      end subroutine fcost


!          program trial

!          real x(4), y(4), cov
!          integer n

!          n=3
!          x(1)=22
!          x(2)=33
!          x(3)=10

!          y(1)=20
!          y(2)=30
!          y(3)=1

!          call covar(x,y,n,cov)
!          print *, cov

!          endprogram trial

          subroutine covar(x,y,n,cov)

          implicit none
          integer, intent(inout) :: n
          real, intent(inout) :: x(n), y(n)
          real, intent(out) :: cov
          real matrix(n,2), m
          integer i
  
          cov = 0.

          !Estimate and remove the mean of each vector
          call mean(x,n,m)
          x(:) = x(:) - m
          call mean(y,n,m)
          y(:) = y(:) - m

          !Order the two vectors in a matrix
          matrix(:,1) = x(:)
          matrix(:,2) = y(:)

          call prod_sum(matrix,n,cov)


          endsubroutine covar


           subroutine mean(dat,n,val)

           implicit none
           integer, intent(inout) :: n
           real, intent(inout) :: dat(n)
           real, intent(out) :: val

           integer i

           !Estimate the mean of a vector
           val = 0.
           do i=1,n
             val = val + (dat(i) / real(n))
           enddo

           end subroutine mean


            subroutine prod_sum(ma,n,cov)

            implicit none
            integer, intent(inout) :: n
            real,intent(inout) :: cov, ma(n,2)
            integer i,j,k


              !Estimate the covariance between the two vectors
              do k=1,n
                   cov = cov + ma(k,1)*ma(k,2)
              enddo

          !  do i=1,2
          !   do j=1,i
          !    cov(i,j) = cov(j,i)
          !   enddo
          !  enddo

           cov = cov / real(n-1)

            endsubroutine prod_sum

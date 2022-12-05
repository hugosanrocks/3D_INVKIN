!*********************************************!
!*    SEISCOPE OPTIMIZATION TOOLBOX          *!
!*********************************************!
! The routine normL2 returns the Euclidian    !
! norm of a vector x of size n                !
!---------------------------------------------!
! INPUT  : integer n                          !
!          real,dimension(n) x                !
! OUTPUT : real norm_x                        !
!---------------------------------------------!
subroutine normL2(n,x,norm_x)
  
  implicit none

  !IN
  integer :: n
  real,dimension(n) :: x
  !IN/OUT
  real :: norm_x
  !Local variables
  integer :: i
  
  norm_x=0.
  do i=1,n
     norm_x=norm_x+x(i)**2
  enddo
  norm_x=sqrt(norm_x)
  
end subroutine normL2

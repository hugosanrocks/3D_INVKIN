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

  integer, intent(in) :: n
  real,dimension(n), intent(in) :: x

  real, intent(out) :: norm_x

  integer :: i
  
  norm_x=0.0
  do i=1,n
     norm_x=norm_x+x(i)**2
  end do
  norm_x=sqrt(norm_x)
  
end subroutine normL2

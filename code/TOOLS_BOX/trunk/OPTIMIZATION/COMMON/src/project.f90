!*********************************************!
!*    SEISCOPE OPTIMIZATION TOOLBOX          *!
!*********************************************!
! The routine project compute the projection  !
! of the vector x into the box defined by     !
! the bound constraints optim%lb and optim%ub !
!---------------------------------------------!
! INPUT  : integer n                          !
!          optim_type optim                   !
! IN/OUT : real,dimension(n) x                !
!---------------------------------------------!
subroutine project(n,optim,x)

  implicit none
  
  include 'optim_type.h'
  
  
  !IN
  integer :: n
  type(optim_type) :: optim
  !IN/OUT
  real,dimension(n) :: x
  !Local variables
  integer :: i
  
  do i=1,n
     if(x(i).gt.optim%ub(i)) then
        x(i)=optim%ub(i)-optim%threshold
     endif
     if(x(i).lt.optim%lb(i)) then
        x(i)=optim%lb(i)+optim%threshold
     endif
  enddo
  
end subroutine project

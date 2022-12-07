!*********************************************!
!*    SEISCOPE OPTIMIZATION TOOLBOX          *!
!*********************************************!
! The routine scalL2 returns the scalar       !
! product between two vectors x and y         ! 
! of size n                                   !
!---------------------------------------------!
! INPUT:  integer n                           !
!         real,dimension(n) x                 !
!         real,dimension(n) y                 !
! OUTPUT: real scal_xy                        !
!---------------------------------------------!
subroutine scalL2(n,x,y,scal_xy)
  implicit none 

  integer, intent(in) :: n
  real,dimension(n), intent(in) :: x,y
  real, intent(out) :: scal_xy
  
  integer :: i

  scal_xy=0.
  do i=1,n
     scal_xy=scal_xy+x(i)*y(i)
  enddo
  
end subroutine scalL2

!*********************************************!
!*    SEISCOPE OPTIMIZATION TOOLBOX          *!
!*********************************************!
! This routine is used to deallocate all the  !
! arrays that have been used during the       !
! PSTD optimization                           !
!---------------------------------------------!
! INPUT/OUT:  optim_type optim                !
!---------------------------------------------!
subroutine finalize_LBFGS(optim)
  
  implicit none
  include 'optim_type.h'  
  
  !IN/OUT
  type(optim_type) :: optim !data structure   
  
  deallocate(optim%sk)
  deallocate(optim%yk)
  deallocate(optim%xk)
  deallocate(optim%grad)
  deallocate(optim%descent)
  
end subroutine finalize_LBFGS

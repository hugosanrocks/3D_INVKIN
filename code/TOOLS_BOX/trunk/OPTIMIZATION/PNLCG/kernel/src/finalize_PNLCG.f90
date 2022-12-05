!*********************************************!
!*    SEISCOPE OPTIMIZATION TOOLBOX          *!
!*********************************************!
! This routine is used to deallocate all the  !
! arrays that have been used during the       !
! PSTD optimization                           !
!---------------------------------------------!
! INPUT/OUT:  optim_type optim                !
!---------------------------------------------!
subroutine finalize_PNLCG(optim)
  
  implicit none
  include 'optim_type.h'  
  
  !IN/OUT
  type(optim_type) :: optim !data structure   
  
  deallocate(optim%xk)
  deallocate(optim%grad)
  deallocate(optim%grad_prev)
  deallocate(optim%descent)
  deallocate(optim%descent_prev)
  
end subroutine finalize_PNLCG

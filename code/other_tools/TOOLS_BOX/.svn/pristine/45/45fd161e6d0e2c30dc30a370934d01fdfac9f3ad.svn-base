!*********************************************!
!*    SEISCOPE OPTIMIZATION TOOLBOX          *!
!*********************************************!
! This routine is used to deallocate all the  !
! arrays that have been used during the       !
! PSTD optimization                           !
!---------------------------------------------!
! INPUT/OUT:  optim_type optim                !
!---------------------------------------------!
subroutine finalize_TRN(optim)
  
  implicit none
  include 'optim_type.h'  
  
  !IN/OUT
  type(optim_type) :: optim !data structure   
  
  deallocate(optim%xk)
  deallocate(optim%grad)
  deallocate(optim%descent)
  deallocate(optim%descent_prev)
  deallocate(optim%residual)
  deallocate(optim%d)
  deallocate(optim%Hd)
  deallocate(optim%eisenvect) 
  
end subroutine finalize_TRN

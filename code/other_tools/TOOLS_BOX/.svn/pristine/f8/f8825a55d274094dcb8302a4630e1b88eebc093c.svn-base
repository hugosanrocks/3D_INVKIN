!*********************************************!
!*    SEISCOPE OPTIMIZATION TOOLBOX          *!
!*********************************************!
! This routine implements a simple            !
! convergence test based on the relative      !
! cost  decrease.                             !
! If the current relative cost is lower than  !
! a value set by the user in optim%conv       !
! then test_conv is returned equal to .true.  !
! otherwise, it is returned equal to .false.  ! 
!---------------------------------------------!
! INPUT:  optim_type optim                    !
! OUTPUT: logical test_conv                   !
!---------------------------------------------!
subroutine std_test_conv(optim,fcost,test_conv)
  
  implicit none
  include 'optim_type.h'  

  !IN
  real :: fcost
  type(optim_type) :: optim !data structure   
  !INT/OUT
  logical :: test_conv
    
  test_conv=((fcost/optim%f0<optim%conv).or.&
       (optim%cpt_iter.ge.optim%niter_max))
  
end subroutine std_test_conv

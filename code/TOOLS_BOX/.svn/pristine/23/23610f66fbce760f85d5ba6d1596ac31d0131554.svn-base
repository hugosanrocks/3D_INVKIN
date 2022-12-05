!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This file contains the routine for the computation  !
! of forcing term following the first formula of      !
! Eisenstat and Walker, see                           !
! Choosing the forcing term in an Inexact Newton      !
! method, Eisenstat and Walker, 1994,                 !
! SIAM Journal on Scientific Computing 17 (1), 16-32  !
!------------------------------------------------------
! INPUT : integer n dimension                         !
!         real,dimension(n) grad current gradient     !
!         character*4 FLAG communication flag         !
! INPUT/OUTPUT : optim_type optim (data structure)    !
!-----------------------------------------------------!
subroutine forcing_term_TRN(n,grad,optim)
  
  implicit none
  include 'optim_type.h'  
  !IN
  integer :: n                               !dimension of the problem
  real,dimension(n) :: grad                  !gradient at the current point
  !IN/OUT 
  type(optim_type) :: optim                  !data structure   
  !Local variable
  real :: eta_save,eta_save_power,norm_eisenvect
  
  !-----------------------------------------------------!
  ! Computation of the forcing term optim%eta following !
  ! the formula                                         !
  !-----------------------------------------------------!
  eta_save=optim%eta  
  optim%eisenvect(:)=grad(:)-optim%residual(:)
  call normL2(n,optim%eisenvect,norm_eisenvect)
  optim%eta=norm_eisenvect/optim%norm_grad_m1
  
  !-----------------------------------------------------!
  ! Additional safeguard if optim%eta is too large      !       
  !-----------------------------------------------------!
  eta_save_power=eta_save**((1.+sqrt(5.))/2.)
  if(eta_save_power>0.1) then
     optim%eta=max(optim%eta,eta_save_power)
  endif
  if(optim%eta>1.) then
     optim%eta=0.9     
  endif
  
end subroutine forcing_term_TRN

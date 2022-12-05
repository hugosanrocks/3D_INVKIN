!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This file contains the routines for the computation !
! of the Newton descent direction using a             !
! preconditioned matrix-free conjugate gradient solver!
! The algorithm is taken from                         !
! Numerical Optimization, Nocedal, 2nd edition, 2006  !
! Algorithm 5.3 p.119                                 !
! The stopping criterion is taken from                !
! Choosing the forcing term in an Inexact Newton      !
! method, Eisenstat and Walker, 1994,                 !
! SIAM Journal on Scientific Computing 17 (1), 16-32  !
! The routine is implemented in a reverse             !
! communication format. The user is requested         !
! - to compute  Hessian-vector products at each       !
!   iteration of the process                          !
! - to apply its preconditioner to the residuals      !
!   at each iteration of the process                  !
! Because of the preconditioning, the routine has been!
! split into three routines, compared to the TRN      !
! version                                             !
! The first routine is descent_PTRN0 for              !
! initialization                                      !   
! The second routine is descent_PTRN1 for the first   !
! part of the conjugate gradient iteration, before    !
! preconditioning.                                    !
! The third routine is descent_PTRN2 for the second   !
! part of the conjugate gradient iteration, after     !
! preconditioning.                                    !

!-----------------------------------------------------!
! descent_PTRN0
!-----------------------------------------------------!
! INPUT : integer n dimension                         !
!         real,dimension(n) grad current gradient     !
!         real,dimension(n) grad_preco preconditioned !
!                           current gradient          !
!         character*4 FLAG communication flag         !
! INPUT/OUTPUT : optim_type optim (data structure)    !
!-----------------------------------------------------!
subroutine descent_PTRN0(n,grad,grad_preco,optim,FLAG)
  
  implicit none
  include 'optim_type.h'
  
  !IN
  character*4 :: FLAG                                      !communication FLAG
  integer :: n                                             !dimension of the problem
  real,dimension(n) :: grad                                !gradient at current point
  real,dimension(n) :: grad_preco                          !preconditioned gradient at current point
  !IN/OUT
  type(optim_type) :: optim
  
  !-----------------------------------------------------!
  ! Initialization of the conjugate gradient process    !
  !-----------------------------------------------------!
  optim%residual(:)=grad(:)
  optim%residual_preco(:)=grad_preco(:)
  optim%d(:)=-1.*optim%residual_preco(:)
  optim%Hd(:)=0.
  optim%descent(:)=0.     
  optim%qk_CG=0.
  optim%hessian_term=0.     
  call normL2(n,optim%residual,optim%norm_residual)
  optim%conv_CG=.false.    
  optim%cpt_iter_CG=0
  call print_info_PTRN(optim,0e0,FLAG)     

end subroutine descent_PTRN0

!-----------------------------------------------------!
! descent_PTRN1
!-----------------------------------------------------!
! INPUT : integer n dimension                         !
!         real,dimension(n) grad current gradient     !
!         character*4 FLAG communication flag         !
! INPUT/OUTPUT : optim_type optim (data structure)    !
!-----------------------------------------------------!
subroutine descent_PTRN1(n,grad,optim,FLAG)
  
  implicit none
  include 'optim_type.h'
  
  !IN
  character*4 :: FLAG                        !communication FLAG
  integer :: n                               !dimension of the problem
  real,dimension(n) :: grad                  !gradient at current point
  !IN/OUT
  type(optim_type) :: optim                  !data structure
  !Local variables
  real :: grad_term
  real,dimension(:),allocatable :: mgrad
  
  !-----------------------------------------------------!
  ! start one conjugate gradient iteration       !
  !-----------------------------------------------------!    
  call scalL2(n,optim%d,optim%Hd,optim%dHd)        
  if(optim%dHd<0.) then
     !-----------------------------------------------------!
     ! if dHd < 0, detection of a negative eigenvalue of   !
     ! the Hessian operator, stop the process              !        
     !-----------------------------------------------------!  
     optim%conv_CG=.true.
     write(21,*) 'Negative curvature'
     if(optim%cpt_iter_CG.eq.0) then
        !-----------------------------------------------------!
        ! if this is the first iteration, thenreturn the      !
        ! opposite of the preconditioned gradient as descent  !
        ! direction: preconditioned steepest descent direction!
        !-----------------------------------------------------!   
        optim%descent(:)=optim%d(:) 
        !-----------------------------------------------------!
        ! if the debug option is activated, compute the       !
        ! quadratic function minimized during the conjugate   !
        ! gradient process (check is this value decresae      !
        ! throughout the CG iterations )                      !
        !-----------------------------------------------------!
        if(optim%debug) then
           call scalL2(n,optim%residual,optim%residual_preco,optim%res_scal_respreco)
           optim%alpha_CG=optim%res_scal_respreco/optim%dHd         
           optim%qkm1_CG=optim%qk_CG
           allocate(mgrad(n))
           mgrad(:)=-1.*grad(:)           
           call scalL2(n,optim%descent,mgrad,grad_term)           
           optim%hessian_term=optim%hessian_term+(optim%alpha_CG**2)*optim%dHd
           optim%qk_CG=-grad_term+0.5*optim%hessian_term  
           deallocate(mgrad)
        endif
     endif
  else         
     !-----------------------------------------------------!
     ! if dHd > 0, then start one conjugate gradient       !
     ! iteration                                           !
     !-----------------------------------------------------!   
     !Update descent direction
     call scalL2(n,optim%residual,optim%residual_preco,optim%res_scal_respreco)
     optim%alpha_CG=optim%res_scal_respreco/optim%dHd
     optim%descent_prev(:)=optim%descent(:)
     optim%descent(:)=optim%descent(:)+optim%alpha_CG*optim%d(:)
     optim%residual(:)=optim%residual(:)+optim%alpha_CG*optim%Hd(:)        
     !STOP HERE and wait for preconditioning
  end if
  
end subroutine descent_PTRN1

!-----------------------------------------------------!
! descent_PTRN2
!-----------------------------------------------------!
! INPUT : integer n dimension                         !
!         real,dimension(n) grad current gradient     !
!         character*4 FLAG communication flag         !
! INPUT/OUTPUT : optim_type optim (data structure)    !
!-----------------------------------------------------!
subroutine descent_PTRN2(n,grad,optim,FLAG)
  
  implicit none
  include 'optim_type.h'
  
  !IN
  character*4 :: FLAG                                 !communication FLAG
  integer :: n                                        !dimension of the problem
  real,dimension(n) :: grad                           !gradient at current point
  !IN/OUT
  type(optim_type) :: optim                           !data structure
  !Local variables
  real :: res_scal_respreco_prev,&
       beta,grad_term,descent_scal_Hd
  real,dimension(:),allocatable :: mgrad
  
  !-----------------------------------------------------!
  ! continue the current conjugate gradient iteration   !
  !-----------------------------------------------------! 
  !Update CG direction
  res_scal_respreco_prev=optim%res_scal_respreco
  call scalL2(n,optim%residual,optim%residual_preco,&
       optim%res_scal_respreco)
  beta=(optim%res_scal_respreco)/(res_scal_respreco_prev)
  optim%d(:)=-1.*optim%residual_preco(:)+beta*optim%d(:)                  
  !Update iteration counter 
  optim%cpt_iter_CG=optim%cpt_iter_CG+1  
  !-----------------------------------------------------!
  ! if the debug option is activated, compute the       !
  ! quadratic function minimized during the conjugate   !
  ! gradient process (check is this value decresae      !
  ! throughout the CG iterations )                      !
  !-----------------------------------------------------!  
  if(optim%debug) then
     optim%qkm1_CG=optim%qk_CG
     allocate(mgrad(n))
     mgrad(:)=-1.*grad(:)
     call scalL2(n,optim%descent,mgrad,grad_term)
     call scalL2(n,optim%descent_prev,optim%Hd,descent_scal_Hd)
     optim%hessian_term=optim%hessian_term+(optim%alpha_CG**2)*optim%dHd+&
          2.*optim%alpha_CG*descent_scal_Hd
     optim%qk_CG=-grad_term+0.5*optim%hessian_term
     deallocate(mgrad)
  endif
  !-----------------------------------------------------!
  ! Check if the Eisenstat stopping critertion is       !
  ! satisfied                                           !
  !-----------------------------------------------------! 
  call normL2(n,optim%residual,optim%norm_residual)
  optim%conv_CG=(&
       (optim%norm_residual.le.(optim%eta*optim%norm_grad)).or.&
       (optim%cpt_iter_CG.ge.optim%niter_max_CG))        
  !-----------------------------------------------------!
  ! Print information on the current CG iteration       !
  !-----------------------------------------------------!
  call print_info_PTRN(optim,0e0,FLAG)     
  
end subroutine descent_PTRN2


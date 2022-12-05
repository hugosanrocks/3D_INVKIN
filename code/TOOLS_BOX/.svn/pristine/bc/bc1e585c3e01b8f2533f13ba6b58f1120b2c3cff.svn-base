!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This routine is the reverse communication mode      !
! preconditioned nonlinear conjugate gradient         !
! algorithm of Dai and Yuan                           !
! Y. DAI AND Y. YUAN, A nonlinear conjugate gradient  !
! method with a strong global convergence property,   !
! SIAM Journal on Optimization, 10 (1999), pp. 177–182!
!                                                     !
! See also Nocedal, Numerical optimization,           !
! 2nd edition p.132                                   !
!                                                     !
! This routine performs an iterative                  !
! minimization of a function f following the          !
! recurrence                                          !
!                                                     !
! x_{0}=x0                                            !
! x_{k+1}=x_k+\alpha_k d_k                            !
!                                                     !
! where the descent direction d_k is                  !
!                                                     !
! d_k=-Q_k \nabla f_k + \beta_k d_{k-1}               !
!                                                     !
! where Q_k       : preconditioner at iteration k     !
!      \nabla f_k : gradient of f in x_k              !
!      \beta_k    : scalar computed through the       !
!                   Dai and Yuan formula              !
!                                                     !
! and alpha_k is the steplength computed through the  !
! common linesearch algorithm of the TOOLBOX          !
!                                                     !
! The first call to the algorithm must be done with   !
! FLAG='INIT'. For this first call, the initial point !
! x0 is given through the variable x. The input       !
! variables fcost and grad, grad_preco must correspond!
! respectively to the misfit, gradient and            ! 
! preconditioned gradient at x0.                      !
!                                                     !
! The reverse communication with the user is          !
! performed through the variable FLAG. This           !
! variable indicates to the user on return what action! 
! he has to do, or the state of the algorithm.        !
! Possible values are                                 !
! - FLAG='GRAD' => the user must compute the cost,    !
!                  gradient and preconditioned        !
!                  gradient at current point x        ! 
! - FLAG='CONV' => a minimizer has been found         !
! - FLAG='NSTE' => a new step is performed            !
! - FLAG='FAIL' => the linesearch has failed          !
!-----------------------------------------------------!
! INPUT  : integer :: n (dimension)                   ! 
!          real fcost (current cost)                  !
!          real,dimension(n) grad                     !
!          real,dimension(n) grad_preco               !
! INPUT/OUTPUT : real,dimension(n) x                  !
!                optim_typ optim (data structure)     !
!                character*4 FLAG (communication)     !
!-----------------------------------------------------!
subroutine PNLCG(n,x,fcost,grad,grad_preco,optim,FLAG)
  
  implicit none
  include 'optim_type.h'  
  
  !IN
  integer  :: n                        !dimension of the problem
  real :: fcost                        !cost associated with x
  real,dimension(n) :: grad,grad_preco !gradient and preconditioned gradient at x 
  !IN/OUT  
  character*4 :: FLAG  
  real,dimension(n) :: x               !current point
  type(optim_type) :: optim            !data structure   
  !Local variable
  logical :: test_conv
  real :: norm_grad
  
  if(FLAG.eq.'INIT') then
     !-----------------------------------------------------!
     ! if FLAG is INIT, call the dedicated initialization  !
     ! subroutine to allocate data structure optim and     !
     ! initialize the linesearch process                   !
     !-----------------------------------------------------!
     call init_PNLCG(n,x,fcost,grad,grad_preco,optim)
     call std_linesearch(n,x,fcost,grad,optim)
     call print_info(n,'CG',optim,fcost,FLAG)
     FLAG='GRAD'     
     optim%nfwd_pb=optim%nfwd_pb+1
     !Store current gradient before the user compute the new one
     optim%grad_prev(:)=grad(:)      
  else
     !-----------------------------------------------------!
     ! else call the linesearch process                    !  
     !-----------------------------------------------------!
     call std_linesearch(n,x,fcost,grad,optim)     
     if(optim%task.eq.'NEW_STEP') then 
        optim%cpt_iter=optim%cpt_iter+1
        !-----------------------------------------------------!
        ! test for convergence                                !
        !-----------------------------------------------------!
        call std_test_conv(optim,fcost,test_conv)
        if(test_conv) then
           FLAG='CONV'
           !print info on current iteration        
           optim%grad(:)=grad(:)
           call print_info(n,'CG',optim,fcost,FLAG)
           call finalize_PNLCG(optim)
        else
           FLAG='NSTE'
           !-----------------------------------------------------!
           ! if a NEW_STEP is taken, compute a new descent       !
           ! direction using current descent, gradient and       !
           ! preconditioned gradient                             !
           !-----------------------------------------------------!
           call descent_PNLCG(n,grad,grad_preco,optim)                    
           !print info on current iteration        
           optim%grad(:)=grad(:)
           call print_info(n,'CG',optim,fcost,FLAG)
        endif
     elseif(optim%task.eq.'NEW_GRAD') then
        !-----------------------------------------------------!
        ! if the linesearch needs a new gradient then ask the !  
        ! user to provide it
        !-----------------------------------------------------!
        FLAG='GRAD'         
        optim%nfwd_pb=optim%nfwd_pb+1
        !Store current gradient before the user compute the new one
        optim%grad_prev(:)=grad(:)        
     elseif(optim%task.eq.'FAILURE!') then        
        !-----------------------------------------------------!
        ! if the linesearch has failed, inform the user       !
        !-----------------------------------------------------!
        FLAG='FAIL'
        !print info on current iteration
        optim%grad(:)=grad(:)
        call print_info(n,'CG',optim,fcost,FLAG)
        call finalize_PNLCG(optim)
     endif
  endif
  
end subroutine PNLCG




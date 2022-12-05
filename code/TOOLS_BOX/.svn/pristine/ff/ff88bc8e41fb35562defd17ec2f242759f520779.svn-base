!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This routine is the reverse communication mode      !
! l-BFGS algorithm. This algorithm is described in    !
! Numerical Optimization, Nocedal, 2nd edition, 2006  !
! Algorithm 7.4 p. 178, Algorithm 7.5 p. 179          !
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
! d_k=-Q_k \nabla f_k                                 !
!                                                     !
! with Q_k        : l-BFGS approximation of the       !
!                   inverse Hessian at iteration k    !
!\nabla f_k       : gradient of f in x_k              !
!                                                     !
! and alpha_k is the steplength computed through the  !
! common linesearch algorithm of the TOOLBOX          !
!                                                     !
! The first call to the algorithm must be done with   !
! FLAG='INIT'. For this first call, the initial point !
! x0 is given through the variable x, and the input   !
! variable fcost and grad must correspond respectively!
! to the misfit and gradient at x0.                   !
!                                                     !
! The reverse communication with the user is          !
! performed through the variable FLAG. This           !
! variable indicates to the user on return what action! 
! he has to do, or the state of the algorithm.        !
! Possible values are                                 !
! - FLAG='GRAD' => the user must compute the cost and !
!                  (preconditioned) gradient at       !
!                  current point x                    !
! - FLAG='CONV' => a minimizer has been found         !
! - FLAG='NSTE' => a new step is performed            !
! - FLAG='FAIL' => the linesearch has failed          !
!-----------------------------------------------------!
! INPUT  : integer :: n (dimension)                   ! 
!          real fcost (current cost)                  !
!          real,dimension(n) grad                     !
! INPUT/OUTPUT : real,dimension(n) x                  !
!                optim_type optim (data structure)    !
!                character*4 FLAG (communication)     !
!-----------------------------------------------------!
subroutine LBFGS(n,x,fcost,grad,optim,FLAG)
  
  implicit none
  include 'optim_type.h'  
  
  !IN
  integer  :: n                         !dimension of the problem
  integer  :: l                         !number of stored pairs for the 
                                        !l-BFGS approximation of the 
                                        !inverse Hessian
  real :: fcost                         !cost associated with x
  real,dimension(n) :: grad             !associated gradient 
  !IN/OUT  
  character*4 :: FLAG  
  real,dimension(n) :: x                !current point
  type(optim_type) :: optim             !data structure   
  !Local variables
  logical :: test_conv
  real :: norm_grad

  if(FLAG.eq.'INIT') then
     !-----------------------------------------------------!
     ! if FLAG is INIT, call the dedicated initialization  !
     ! subroutine to allocate data structure optim and     !
     ! initialize the linesearch process                   !
     !-----------------------------------------------------!
     call init_LBFGS(n,optim%l,x,fcost,grad,optim)
     call std_linesearch(n,x,fcost,grad,optim)
     call print_info(n,'LB',optim,fcost,FLAG)
     FLAG='GRAD'     
     optim%nfwd_pb=optim%nfwd_pb+1
  else
     !-----------------------------------------------------!
     ! else call the linesearch process                    !  
     !-----------------------------------------------------!
     call std_linesearch(n,x,fcost,grad,optim)     
     if(optim%task.eq.'NEW_STEP') then
        !-----------------------------------------------------!
        ! test for convergence                                !
        !-----------------------------------------------------!
        optim%cpt_iter=optim%cpt_iter+1                   
        call std_test_conv(optim,fcost,test_conv)
        if(test_conv) then
           FLAG='CONV'           
           !print info on current iteration
           optim%grad(:)=grad(:)
           call print_info(n,'LB',optim,fcost,FLAG)        
           call finalize_LBFGS(optim)
        else
           FLAG='NSTE'
           !-----------------------------------------------------!
           ! if a NEW_STEP is taken, compute a new descent       !
           ! direction using current gradient and l-BFGS         !
           ! approximation of the inverse Hessian                !
           ! preconditioned gradient                             !
           !-----------------------------------------------------!
           !LBFGS update        
           call update_LBFGS(n,optim%cpt_lbfgs,optim%l,x,grad,optim%sk,optim%yk)
           !Computation of the new descent direction
           call descent_LBFGS(n,grad,optim%sk,optim%yk,optim%cpt_lbfgs,&
                optim%l,optim%descent)         
           !LBFGS store
           call save_LBFGS(n,optim%cpt_lbfgs,optim%l,x,grad,optim%sk,optim%yk)           
           !print info on current iteration
           optim%grad(:)=grad(:)
           call print_info(n,'LB',optim,fcost,FLAG)        
        endif
     elseif(optim%task.eq.'NEW_GRAD') then 
        !-----------------------------------------------------!
        ! if the linesearch needs a new gradient then ask the !  
        ! user to provide it
        !-----------------------------------------------------!
        FLAG='GRAD'         
        optim%nfwd_pb=optim%nfwd_pb+1        
     elseif(optim%task.eq.'FAILURE!') then        
        !-----------------------------------------------------!
        ! if the linesearch has failed, inform the user       !
        !-----------------------------------------------------!
        FLAG='FAIL'
        !print info on current iteration
        optim%grad(:)=grad(:)
        call print_info(n,'LB',optim,fcost,FLAG)
        call finalize_LBFGS(optim)
     endif
  endif
  
end subroutine LBFGS




!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This routine is used for the initialization of the  !
! reverse communication mode preconditioned nonlinear !
! conjugate gradient algorithm from the TOOLBOX.      !
!                                                     !
! The parameters for the linesearch are set here, as  !
! well as the data structure allocations required for !
! the use of the algorithm                            !
!-----------------------------------------------------!
! INPUT : integer n, dimension of the problem         !
!       : real fcost                                  !
!       : real,dimension(n) x,grad,grad_preco         !
!                           first iterate, gradient   !
!                           preconditioned gradient   !
! INPUT/OUTPUT : optim_type optim data structure      ! 
!-----------------------------------------------------!
subroutine init_PNLCG(n,x,fcost,grad,grad_preco,optim)
  
  implicit none
  include 'optim_type.h'  
  !IN
  integer :: n
  real :: fcost
  real,dimension(n) :: x,grad,grad_preco
  !IN/OUT
  type(optim_type) :: optim !data structure   

  !---------------------------------------!
  ! set counters                          !
  !---------------------------------------!
  optim%cpt_iter=0
  optim%f0=fcost
  optim%nfwd_pb=0

  !---------------------------------------!
  ! initialize linesearch parameters      !
  ! by default, the max number of         !
  ! linesearch iteration is set to 20     !
  ! and the initial steplength is set to 1!
  !---------------------------------------! 
  optim%m1=1e-4 ! Wolfe conditions parameter 1 (Nocedal value)
  optim%m2=0.9  ! Wolfe conditions parameter 2 (Nocedal value)
  optim%mult_factor=10 ! Bracketting parameter (Gilbert value)
  optim%fk=fcost
  optim%nls_max=20 ! max number of linesearch
  optim%cpt_ls=0
  optim%first_ls=.true.
  optim%alpha=1. ! first value for the linesearch steplength 
  
  !---------------------------------------!
  ! memory allocations                    !
  !---------------------------------------!
  allocate(optim%grad_prev(n))
  allocate(optim%descent_prev(n))
  allocate(optim%xk(n))
  optim%xk(:)=x(:)
  allocate(optim%grad(n))
  optim%grad(:)=grad(:)
  allocate(optim%descent(n))
  
  !---------------------------------------!
  ! first descent direction               !
  !---------------------------------------!
  optim%descent(:)=-1.*grad_preco(:)
  

end subroutine init_PNLCG

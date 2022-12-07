!*****************************************************!
!*          SEISCOPE OPTIMIZATION TOOLBOX            *!
!*****************************************************!
! This file contains the routine for the              !
! initialization of the truncated Newton algorithm    !
!------------------------------------------------------
! INPUT : integer n dimension                         !
!         real,dimension(n) x current point           !
!         real,dimension(n) grad current gradient     !
!         character*4 FLAG communication flag         !
! INPUT/OUTPUT : optim_type optim (data structure)    !
!-----------------------------------------------------!
subroutine init_PTRN(n,x,fcost,grad,optim)
  
  implicit none
  include 'optim_type.h'  
  !IN
  integer :: n 
  real :: fcost
  real,dimension(n) :: x,grad
  !IN/OUT
  type(optim_type) :: optim !data structure   
    
  !---------------------------------------!
  ! set counters                          !
  !---------------------------------------!
  optim%cpt_iter=0
  optim%f0=fcost
  optim%nfwd_pb=0
  optim%nhess=0
  optim%eta=0.9  
  !optim%eta=0.1
  
  !---------------------------------------!
  ! initialize linesearch parameters      !
  ! by default, the max number of         !
  ! linesearch iteration is set to 20     !
  ! and the initial steplength is set to 1!
  !---------------------------------------! 
  optim%m1=1e-4  ! Wolfe conditions parameter 1 (Nocedal value)
  optim%m2=0.9   ! Wolfe conditions parameter 2 (Nocedal value)
  optim%mult_factor=10 !Bracketting parameter (Gilbert value)
  optim%fk=fcost
  optim%nls_max=20 ! max number of linesearch
  optim%cpt_ls=0
  optim%first_ls=.true.
  optim%alpha=1.   ! first value for the linesearch steplength

  !---------------------------------------!
  ! memory allocations                    !
  !---------------------------------------!
  allocate(optim%xk(n))
  optim%xk(:)=x(:)     
  allocate(optim%grad(n))
  optim%grad(:)=grad(:)
  allocate(optim%descent(n))
  allocate(optim%descent_prev(n))
  allocate(optim%residual(n))
  allocate(optim%residual_preco(n))
  allocate(optim%d(n))
  allocate(optim%Hd(n))
  allocate(optim%eisenvect(n))

  !---------------------------------------!
  ! norm of the first gradient            !
  !---------------------------------------!
  call normL2(n,grad,optim%norm_grad)
   
  
end subroutine init_PTRN
